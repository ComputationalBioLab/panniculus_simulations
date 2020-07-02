# Here, mu (often, G) is the shear moddulus and kappa is the bulk modulus
# mu = E/(2*(1+nu))
# kappa = E/(3*(1-2*nu))
# E = 2*G*(1+nu)
    
#soft_tissue
print('--------SOFT TISSUE (MPa where applicable)---------')
#Parameters from Isvilanonda et al (Ogden material)
mu_i = 9.685*1e-9*1e6 #Pa
alpha_i = 24.326 #dimensionless
mu_i_2 = 1e-3*1e6 #Pa 
alpha_i_2 = 2.0 #dimensionless
#Table parameters (published)
G_st = 0.5*(mu_i*alpha_i+mu_i_2*alpha_i_2)
nu_ST = 0.495
E_st = 2*G_st*(1+nu_ST)
K_st =E_st/(3.*(1.-2.*nu_ST))
print('G is ', G_st*1e-6, 'nu is ',nu_ST, 'E is ',E_st*1e-6, 'K is ', K_st*1e-6)
print('---------------------------------------------------')
print('----------TENDON (MPa where applicable)------------')
#Friedman reference. Table 7.1, page 126
E_tendon = 0.1945*1e6#In Pa
nu_tendon = 0.495
mu_tendon = E_tendon/(2.*(1.+nu_tendon)) 
kappa_tendon= E_tendon/(3.*(1-2*nu_tendon))
print('G is ', mu_tendon*1e-6, 'nu is ',nu_tendon, 'E is ',E_tendon*1e-6, 'K is ', kappa_tendon*1e-6)
print('---------------------------------------------------')
print('----------PC (MPa where applicable)----------------')
#From Bryan's book chapter citing Koo et al
G_PC = 71000 #This one is taken from the experiment on subject in Figure 3 -> G0=71 KPa
nu_PC = 0.457#From Hansang et al, 0.457 mentioned just befor ethe discussion
kappa_PC = 2*G_PC*(1+nu_PC)/(3.*(1.-2.*nu_PC))
E_PC = 2*G_PC*(1+nu_PC)
print('G is ', G_PC*1e-6, 'nu is ',nu_PC, 'E is ',E_PC*1e-6, 'K is ', kappa_PC*1e-6)
print('---------------------------------------------------')

print('----------SKIN (MPa where applicable)--------------')
E_skin = 0.970853*1e6 #0.970853 MPa Table 7.1 of Friedman et al
nu_skin = 0.495 #Table 7.1 of Friedman et al
mu_skin = E_skin/(2.*(1.+nu_skin)) 
kappa_skin = E_skin/(3.*(1-2*nu_skin))
print('G is ', mu_skin*1e-6, 'nu is ',nu_skin, 'E is ',E_skin*1e-6, 'K is ', kappa_skin*1e-6)
print('---------------------------------------------------')
print('----------BONE (MPa where applicable)--------------')
E_bone = 7*1e9 #Table 7.1 of Friedman et al
nu_bone = 0.3 #Table 7.1 of Friedman et al
mu_bone = E_bone/(2.*(1.+nu_bone)) 
kappa_bone = E_bone/(3.*(1-2*nu_bone))
print('G is ', mu_bone*1e-6, 'nu is ',nu_bone, 'E is ',E_bone*1e-6, 'K is ', kappa_bone*1e-6)
print('---------------------------------------------------')


