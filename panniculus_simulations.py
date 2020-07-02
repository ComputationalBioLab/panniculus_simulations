#!/usr/bin/env python3

# This code runs a quasi-static hyeperelastic (neo-Hookean) simulation of a 
# heel under compression from a reaction force.
# 
# To run this code, FENICS must be installed. 
# It was written with FENICS 2019.1 installed on Ubuntu Linux 18.04 OS.
#
# IMPORTANT NOTE AND CREDITS:
#
# This code has been adapted from various FENICS forum posts, demos and 
# tutorials. Two, in particular, are worth mentioning and are hereby credited:
# 1. Hyperelasticity FENICS demo:
# https://fenicsproject.org/docs/dolfin/latest/python/demos/hyperelasticity/demo_hyperelasticity.py.html
# 2. Forum post by Marc(marchirschvogel):
# https://fenicsproject.discourse.group/t/nonlinear-hyperelasticity-in-dolfin-x-minimal-example-and-a-few-issues/3098
# 3. Hyperelasticity tutorial
# https://fenics-handson.readthedocs.io/en/latest/elasticity/doc.html#reference-solution
#
# Adaptation of the above-mentioned resources for this particular problem
# was written by Alberto Corrias - June 2020

import time
import sys, os, subprocess, time
import math

from dolfin import *
import numpy as np


mesh_dir = 'meshes/'
results_dir = 'output_files/'
parameters["form_compiler"]["quadrature_degree"] = 4

def RunSimulation(mesh_name,simulation_has_PC):
    start = time.time()

    mesh_filename = mesh_dir + mesh_name + '.xml'

    if (simulation_has_PC == True):
        results_filename = results_dir + mesh_name + '12N.xdmf';#For 12 N results
    else:
        results_filename = results_dir + mesh_name + '-NOPC-12N.xdmf';#For 12 N results

    # in which interval to write results
    write_results_every = 1
    # mpi communicator (for parallel runs)
    comm = MPI.comm_world
    
    # read in mesh (materials are tagged within the mesh)
    mesh = Mesh(mesh_filename)
    materials = MeshFunction('size_t', mesh, 3, mesh.domains())
    all_domains = materials.array();

    # Elasticity parameters - See material_params_calculations.py for more info
    # where E is the Young's modulus and nu is the Poisson ratio
    # All values of E are in Pa
    E_soft_tisse  = 2990.0; 
    nu_soft_tissue =0.495;
    E_tendon = 0.1945*1e6 
    nu_tendon = 0.495
    E_skin = 0.970853*1e6 
    nu_skin = 0.495
    E_PC = 0.206894*1e6 
    nu_PC = 0.457 
    E_bone= 7e9 
    nu_bone= 0.3
    
    #Derived from E and nu, mu is shear modulus and K is bulk modulus
    mu_soft_tissue = E_soft_tisse/(2.*(1.+nu_soft_tissue))
    mu_skin = E_skin/(2.*(1.+nu_skin))
    mu_tendon = E_tendon/(2.*(1.+nu_tendon)) 
    mu_PC= E_PC/(2.*(1.+nu_PC))
    mu_bone= E_bone/(2.*(1.+nu_bone))
    kappa_soft_tissue =E_soft_tisse/(3.*(1.-2.*nu_soft_tissue))
    kappa_skin = E_skin/(3.*(1-2*nu_skin))
    kappa_tendon= E_tendon/(3.*(1-2*nu_tendon))
    kappa_PC=E_PC/(3.*(1-2*nu_PC))
    kappa_bone= E_bone/(3.*(1-2*nu_bone))
    
    class neo_hook_param(UserExpression):
        def __init__(self, materials, p_soft_tissue, p_bone, p_AT, p_skin, p_PC, **kwargs):
            super().__init__(kwargs)
            self.materials = materials
            self.p_soft_tissue = p_soft_tissue
            self.p_bone = p_bone
            self.p_AT = p_AT
            self.p_skin = p_skin
            self.p_PC = p_PC
            
        #Domain codes for PC meshes (7 domains)
        #0: PC
        #1: Bone
        #2: Achilles tendon
        #3: Soft tissue in front
        #4: Skin
        #5: Skin (bottom part)
        #6: Soft tissue        
        def eval_cell(self, values, x, cell):
                if self.materials[cell.index] == 0:
                    if (simulation_has_PC == True):
                        values[0] = self.p_PC
                    else:
                        values[0] = self.p_soft_tissue
                if self.materials[cell.index] == 1:
                        values[0] = self.p_bone
                if self.materials[cell.index] == 2:
                        values[0] = self.p_AT
                if self.materials[cell.index] == 3:
                        values[0] = self.p_soft_tissue
                if self.materials[cell.index] == 4:
                        values[0] = self.p_skin
                if self.materials[cell.index] == 5:
                        values[0] = self.p_skin
                if self.materials[cell.index] == 6:
                        values[0] = self.p_soft_tissue
                        
        def value_shape(self):
            return ();
        
    #Vector, fucntions and function spaces
    V_u = VectorFunctionSpace(mesh, "CG", 1)
    
    du    = TrialFunction(V_u)            # Incremental displacement
    var_u = TestFunction(V_u)             # Test function
    u     = Function(V_u, name='Displacement')
    del_u = Function(V_u) # disp increment for Newton solver
    
    
    #For post-processing
    V_strain     = TensorFunctionSpace(mesh, "DG", 0)#Used for GL strain
    V_stress     = TensorFunctionSpace(mesh, "DG", 0)#Cauchy stress
    V_von_mises  = FunctionSpace(mesh,"DG",0)
    V_equiv_strain = FunctionSpace(mesh,"DG",0)
    
    # print infos on problem size (this is # of nodes * 3 for 3D problems)
    if comm.rank == 0:
        print("Number of degrees of freedom: %i" % (V_u.dim()))
        sys.stdout.flush()
    
    ##############################################################################
    #BOUNDARY CONDITIONS
    #First, figure out lower y coordinate. This anatomy is oriented downward and
    #we clamp the top of the foot (lowest y coordinate) and max z coordinate
    min_y_coord = 1e9;
    max_y_coord = -1e9
    min_z_coord = 1e9
    max_z_coord = -1e9
    for n in mesh.coordinates():
        y = n[1]
        z = n[2]
        if (y<min_y_coord):
            min_y_coord = y
        if (y>max_y_coord):
            max_y_coord = y
        if (z<min_z_coord):
            min_z_coord = z
        if (z>max_z_coord):
            max_z_coord = z
    
    def clamped_boundary(x, on_boundary):
        return (x[1] < (min_y_coord + 1e-4)) and on_boundary
    
    def clamped_boundary_2(x, on_boundary):
        return (x[2]> (max_z_coord - 1e-4)) and on_boundary
    
    # Fixed boundary condition for displacement
    fixed = Constant((0.0, 0.0, 0.0))
    bc0 = DirichletBC(V_u, fixed, clamped_boundary)
    bc1 = DirichletBC(V_u, fixed, clamped_boundary_2) 
    dir_bcs = [bc0,bc1]
    ##############################################################################
    
    ##############################################################################
    #APPLIED REACTION FORCE TO HE BOTTOM OF THE HEEL
    # Load stepping. Change only nstep, if needed. 
    # Magnitude is in the load function
    T       = 1.0
    nstep   = 10
    dt      = T/nstep
    
    area_of_contact = 257.59*10**-6 #from Solidworks, 224.23 mm^2 ... 257.59 if cutoff at 23 mm
    force = 12 #change to 15 for 15 N simulations
    max_surf_traction = force/area_of_contact;
    def applied_traction(t):
        return Constant(max_surf_traction*t/T)
    
    class SurfaceBoundary(SubDomain):
        def inside(self, x, on_boundary):
            cutoff_for_applying_force = 0.023+min_y_coord
            return on_boundary and x[1] > cutoff_for_applying_force
            
    # Mark facets where apply surface traction with 1
    boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    surface_boundary = SurfaceBoundary()
    surface_boundary.mark(boundary_parts, 1)
    ds4 = ds(subdomain_data=boundary_parts, subdomain_id=1) # Area where to apply
    
    # direction of application (energy function will dot to this)
    n0 = Constant((0.0,-1.0, 0.0));
    ##############################################################################
    
    
    I = Identity(len(u))  # Identity tensor, useful for many functions below
    
    #2nd Piola-Kirchhoff stress, r is displacement
    def PK2_stress(r):
        #S = 2.*diff(psi,C) This is also possible
        
        F = I + grad(r) # deformation gradient
        C = variable(F.T*F)   # right Cauchy-Green tensor
        
        #mu (or G) is the Shear modulus
        mu = neo_hook_param(materials,mu_soft_tissue,mu_bone, mu_tendon, mu_skin, mu_PC, degree=0)
        #kappa = E/3*(1.-2.*nu) is the bulk modulus
        kappa =  neo_hook_param(materials,kappa_soft_tissue, kappa_bone, kappa_tendon, kappa_skin,kappa_PC, degree=0)
        # Cauchy-Green invariants
        #Ic   = tr(C) Unused
        #IIc  = 0.5*(tr(C)**2. - tr(C*C)) Unused
        IIIc = det(C)
        
        #Formulation from Korelc reference of the PK2 stress, page 142 - isochoric volumetric split
        Jf = sqrt(IIIc)
        S = mu*Jf**(-2./3.)*(I-(1.0/3.0)*tr(C)*inv(C))+(kappa/2.0)*(Jf**2 - 1)*inv(C)
        return S
    
    # internal virtual work
    def dW_total(r):
        F = I + grad(r) # deformation gradient
        J = det(F)            # determinant of deformation gradient
        C = variable(F.T*F)   # right Cauchy-Green tensor
        
        dW_internal =  0.5*inner(PK2_stress(r),derivative(C, r, var_u))*dx #Note: the 0.5 is because PK2=2*dW/dC
        
        # Traction at boundary ds4 
        bF = J*dot(inv(F).T, applied_traction(t)*n0)
        dW_external = inner(bF, var_u)*ds4
        
        return (dW_internal - dW_external)
   
    #Given r (displacement), compute and returns stresses and deformations
    def post_process(r):
        F = I + grad(r)
        J = det(F)
        PK2 = PK2_stress(r)
        cauchy_stress = (1./J) * F * PK2 * F.T
        
        #Von Mises
        deviatoric_stress = cauchy_stress - (1./3.)*tr(cauchy_stress)*I
        von_mises_stress = sqrt((3./2)*inner(deviatoric_stress, deviatoric_stress))
        
        #Green Lagrange strain
        E = 0.5*(F.T*F-I)
        
        #equivalent strain
        deviatoric_strain = E - (1./3.)*tr(E)*I        
        equiv_strain =  (sqrt((2./3)*inner(deviatoric_strain, deviatoric_strain)))
        
        return cauchy_stress, von_mises_stress, E, equiv_strain
    
    #output file options    
    results_xdmf_file = XDMFFile(results_filename)
    results_xdmf_file.parameters["flush_output"] = True
    results_xdmf_file.parameters["functions_share_mesh"] = True
    results_xdmf_file.parameters["rewrite_function_mesh"] = False
    
        
    # load/time stepping
    interval = np.linspace(0, T, nstep+1)
    
    for (N, dt) in enumerate(np.diff(interval)):
        
        t = interval[N+1]
    
        # nonlinear variational form: internal minus external virtual work
        varform = dW_total(u)
    
        # Jacobian of nonlinear variational form
        jac = derivative(varform, u, du)
    
        # Newton loop
        struct_res_norm = 1.0
        struct_inc_norm = 1.0
        tol_struct_res = 1.0e-8
        tol_struct_inc = 1.0e-8
        it = 0
        maxiter = 25
        
        # assemble system
        K, r = assemble_system(jac, -varform, dir_bcs)
    
        # get initial norms
        struct_disp_norm_0 = norm(u.vector(),'l2')
        struct_res_norm_0 = norm(r,'l2')
        
        if comm.rank == 0:
            print("Predictor (disp 2-norm %.4e) yields residual 2-norm:\n      %.4e" % (struct_disp_norm_0,struct_res_norm_0))
            print("iter  struct res 2-norm  struct incr 2-norm")
            sys.stdout.flush()
        
        #newton loop
        while it < maxiter:
            
            it += 1
            
            # solve linearized system
            solve(K, del_u.vector(), r, 'mumps')# superlu also works
            
            # update solution
            u.vector().axpy(1.0, del_u.vector())
            
            # assemble system
            K, r = assemble_system(jac, -varform, dir_bcs)
    
            struct_res_norm = norm(r,'l2')
            struct_inc_norm = norm(del_u.vector(),'l2')
    
            if comm.rank == 0:
                print("%i     %.4e         %.4e" % (it, struct_res_norm, struct_inc_norm))
                sys.stdout.flush()
            
            # check if converged
            if struct_res_norm <= tol_struct_res and struct_inc_norm <= tol_struct_inc:
                break
        
        else:#This part is executed only if no break happens during the while loop    
            if comm.rank == 0:
                print("Newton did not converge after %i iterations!" % (it))
                sys.stdout.flush()
            sys.exit()#Newton iterations failed to converge, exit
            
        # write results every write_results_every steps
        if (N+1) % write_results_every == 0:
            
            # Save solution to XDMF format
            results_xdmf_file.write(u, t)
            
            #post-process
            cauchy_stress, von_mises, gl_strain, equiv_strain = post_process(u);
    
            vm = project(von_mises, V_von_mises) 
            vm.rename("VonMises", "VonMises")
            results_xdmf_file.write(vm, t)
            
            eq_s = project(equiv_strain,V_equiv_strain)
            eq_s.rename("EquivStrain","EquivStrain")
            results_xdmf_file.write(eq_s,t)
            
            gl_strn = project(gl_strain,V_strain)
            gl_strn.rename("GreenLagrangeStrain","GreenLagrangeStrain")
            results_xdmf_file.write(gl_strn, t)
            
            cauchy_strs = project(cauchy_stress, V_stress)
            cauchy_strs.rename("CauchyStress","CauchyStress")
            results_xdmf_file.write(cauchy_strs, t)
    
        if comm.rank == 0: # only proc 0 should print this
            print("### LOAD STEP %i / %i successfully completed, LOAD: %.4f" % (N+1,nstep,applied_traction(t)))
            print("-----------------------------------------------------------")
            sys.stdout.flush()
    
    
    if comm.rank == 0: # only proc 0 should print this
        print('Time for computation: %.1f min' % ( (time.time()-start)/60. ))
        sys.stdout.flush()
        
#FOR CONVERGENCE TESTS
mesh_names = ['MESH-PC-COARSE.geo','MESH-PC-FINE.geo', 'MESH-1p2-PC.geo', 'MESH099-PC.geo', 'MESH099-PC.geo']
has_pcs = [True,True,True,True,False] #last mesh without PC, convergence with PC only

#For single simulations...
#mesh_names = ['MESH-PC-COARSE.geo']
#has_pcs = [True]
for i in range(0,len(mesh_names)):
        RunSimulation(mesh_names[i],has_pcs[i])
