

from fenics import *
from dolfin import *
import math
import matplotlib.pyplot as plt

def compute_centroid(v1,v2,v3,v4):
    ret = [0.0,0.0,0.0]
    ret[0] = 0.25* (v1[0] + v2[0] + v3[0] + v4[0])
    ret[1] = 0.25* (v1[1] + v2[1] + v3[1] + v4[1])
    ret[2] = 0.25* (v1[2] + v2[2] + v3[2] + v4[2])
    return ret

def compute_distance(v1,v2):
    return sqrt((v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + + (v1[2]-v2[2])**2)
        
mesh_dir = 'meshes/';
mesh_names = ['MESH_1.geo','MESH_2.geo', 'MESH_3.geo', 'MESH_4.geo', 'MESH_4.geo']
csv_filenames = ['MESH_1.geo','MESH_2.geo', 'MESH_3.geo', 'MESH_4.geo', 'MESH_4-NOPC.geo']
csv_dir = 'post_process/'

#point (centroid) of max vequiv strain
target_for_convergence = [0.016989667,0.0296650675,0.051473093500000004]
#For plotting
equivs_at_targets = [];
num_els = [];
try:    
    output_file = open('post_process_data.txt', "w")
except:
    print('ERROR: Unable to open  output file')
    exit(); 
    
for mesh_index in range(0,len(mesh_names)):
    mesh_name = mesh_names[mesh_index]
    filename = mesh_dir + mesh_name + '.xml'

    csv_point_name = csv_dir + csv_filenames[mesh_index] + '.csv'
    csv_cell_name = csv_dir + csv_filenames[mesh_index] + '-CELL.csv'

    mesh = Mesh(filename)
    mesh.init()
    materials = MeshFunction('size_t', mesh, 3, mesh.domains())
    #Visualize domains only
    #domains_file =  File(mesh_dir + mesh_name + '.pvd') 
    #domains_file << materials
    
    #Some info in array format
    all_cells = mesh.cells();
    all_domains = materials.array();
    all_coords= mesh.coordinates();
    num_els.append(len(all_cells))
    
    target_cell=-1
    cell_counter = 0
    min_dist_from_target = 1e9
    for cell in all_cells:
        vv1 = all_coords[ cell[0] ]
        vv2 = all_coords[ cell[1] ]
        vv3 = all_coords[ cell[2] ]
        vv4 = all_coords[ cell[3] ]
        distance = compute_distance(target_for_convergence, compute_centroid(vv1,vv2,vv3,vv4))
        if distance<min_dist_from_target:
            target_cell = cell_counter
            min_dist_from_target = distance
        cell_counter = cell_counter + 1
        
    try:
        csv_point_file = open(csv_point_name, "r")
    except:
        print('ERROR: Unable to open input file ' + csv_point_name)
        exit();
    try:
        csv_cell_file = open(csv_cell_name, "r")
    except:
        print('ERROR: Unable to open input file ' + csv_cell_name)
        exit();    
    
    line_counter = 0;
    cauchy_indices = []
    green_inidices = []
    equiv_strain_idx = -1
    von_mises_idx  = -1
    vtk_index = -1
    v1_index = -1
    v2_index = -1
    v3_index = -1
    v4_index = -1
    max_cauchy = []
    min_cauchy = []
    max_green= []
    min_green= []
    max_equiv_strain = -1e9
    min_equiv_strain = 1e9
    max_vonmises = -1e18
    min_vonmises = 1e18
    vtk_indices = []
    for i in range(0,9):
        cauchy_indices.append(-10)
        green_inidices.append(-10)
        max_cauchy.append(-1000000.6)
        min_cauchy.append(1e18)
        max_green.append(-1e9)
        min_green.append(1e9)
    
    threshold_equiv_strain = 2.0
    threshold_von_mises = 1e5
    high_von_mises_counter = 0
    high_strain_counter = 0
    for line in csv_cell_file:
        split_line = line.split(',')
        
        if (line_counter==0):#Figure out where things are in the file
            for index in range(0,len(split_line)):
                if split_line[index].find('Cauchy') != -1 :
                    if (split_line[index].find('0') != -1): cauchy_indices[0] = index
                    if (split_line[index].find('1') != -1): cauchy_indices[1] = index
                    if (split_line[index].find('2') != -1): cauchy_indices[2] = index
                    if (split_line[index].find('3') != -1): cauchy_indices[3] = index
                    if (split_line[index].find('4') != -1): cauchy_indices[4] = index
                    if (split_line[index].find('5') != -1): cauchy_indices[5] = index
                    if (split_line[index].find('6') != -1): cauchy_indices[6] = index
                    if (split_line[index].find('7') != -1): cauchy_indices[7] = index
                    if (split_line[index].find('8') != -1): cauchy_indices[8] = index
                if split_line[index].find('Green') != -1 :
                    if (split_line[index].find('0') != -1): green_inidices[0] = index
                    if (split_line[index].find('1') != -1): green_inidices[1] = index
                    if (split_line[index].find('2') != -1): green_inidices[2] = index
                    if (split_line[index].find('3') != -1): green_inidices[3] = index
                    if (split_line[index].find('4') != -1): green_inidices[4] = index
                    if (split_line[index].find('5') != -1): green_inidices[5] = index
                    if (split_line[index].find('6') != -1): green_inidices[6] = index
                    if (split_line[index].find('7') != -1): green_inidices[7] = index
                    if (split_line[index].find('8') != -1): green_inidices[8] = index          
                if split_line[index].find('Equiv') != -1 : equiv_strain_idx = index
                if split_line[index].find('Von') != -1 : von_mises_idx = index
                if split_line[index].find('vtkOriginalIndices') != -1 : vtk_index = index
                if split_line[index].find('Point Index 0') != -1 : v1_index = index
                if split_line[index].find('Point Index 1') != -1 : v2_index = index
                if split_line[index].find('Point Index 2') != -1 : v3_index = index
                if split_line[index].find('Point Index 3') != -1 : v4_index = index
        else:
            for n in range(0,9):
                if (float(split_line[cauchy_indices[n]])>max_cauchy[n]): max_cauchy[n] = float(split_line[cauchy_indices[n]])
                if (float(split_line[green_inidices[n]])>max_green[n]): max_green[n] = float(split_line[green_inidices[n]])
                if (float(split_line[cauchy_indices[n]])<min_cauchy[n]): min_cauchy[n] = float(split_line[cauchy_indices[n]])
                if (float(split_line[green_inidices[n]])<min_green[n]): min_green[n] = float(split_line[green_inidices[n]])
            
            if (float(split_line[equiv_strain_idx])>max_equiv_strain): 
                max_equiv_strain = float(split_line[equiv_strain_idx])
#                vx1 = all_coords[int(split_line[v1_index])]
#                vx2 = all_coords[int(split_line[v2_index])]
#                vx3 = all_coords[int(split_line[v3_index])]
#                vx4 = all_coords[int(split_line[v4_index])]
#                centr = compute_centroid(vx1,vx2,vx3,vx4)
#                print(line_counter)
            if (float(split_line[equiv_strain_idx])<min_equiv_strain): min_equiv_strain = float(split_line[equiv_strain_idx])
            if (float(split_line[equiv_strain_idx])>threshold_equiv_strain): high_strain_counter = high_strain_counter + 1
                
            if (float(split_line[von_mises_idx])>max_vonmises): max_vonmises = float(split_line[von_mises_idx])
            if (float(split_line[von_mises_idx])<min_vonmises): min_vonmises = float(split_line[von_mises_idx])
            if (float(split_line[von_mises_idx])>threshold_von_mises): high_von_mises_counter = high_von_mises_counter + 1
            
            vtk_indices.append(split_line[vtk_index])
            if (line_counter == (target_cell +1)):
                equivs_at_targets.append(float(split_line[equiv_strain_idx]))
        line_counter = line_counter + 1
    low_y_coord = 1e9;
    max_y_coord = -1e9
    for n in mesh.coordinates():
        y = n[1]
        if (y<low_y_coord):
            lower_clamped_y_coord = y
        if (y>max_y_coord):
            max_y_coord = y
    
    csv_point_file.close()
    csv_cell_file.close()
    print('-----------------------------------------')
    print('Mesh ',mesh_name, ' has ', len(mesh.coordinates()), ' nodes, ', len(all_cells), ' elements, ',max(all_domains)+1, ' domains.')
    print('CSV file has ',vtk_indices[len(vtk_indices)-1].rstrip("\n"),' cells and ', line_counter-1, ' rows' )
    print('Von mises range: ',min_vonmises*1e-6, ' - ', max_vonmises*1e-6, ' MPa' )
    print('Equivalent strain range: ',min_equiv_strain, ' - ', max_equiv_strain, )
    print('Number of Von mises above threshold: ',high_von_mises_counter,' which is ',100*high_von_mises_counter/len(all_cells),' percent')
    print('Number of equivalent strain above threshold: ',high_strain_counter,' which is ',100*high_strain_counter/len(all_cells),' percent')
    print('-----------------------------------------')
    
    output_string = 'Analysing ' + csv_cell_name + '\n'
    output_string =  output_string + 'Mesh ' + str(mesh_name) + ' has ' + str(len(mesh.coordinates()))+ ' nodes, ' + str(len(all_cells)) + ' elements, ' + str(max(all_domains)+1) + ' domains. \n'
    output_string =  output_string + 'CSV file has ' + str(vtk_indices[len(vtk_indices)-1].rstrip("\n")) + ' cells and '+ str(line_counter-1) + ' rows \n'
    output_string =  output_string + 'Von mises range: ' + str(min_vonmises*1e-6) + ' - '+ str(max_vonmises*1e-6) + ' MPa\n'
    output_string =  output_string + 'Equivalent strain range: ' + str(min_equiv_strain)+ ' - '+ str(max_equiv_strain) + '\n'
    output_string =  output_string + 'Number of Von mises above ' + str(threshold_von_mises) + ' Pa threshold: ' + str(high_von_mises_counter) + ' which is ' + str(100*high_von_mises_counter/len(all_cells))+' percent \n'
    output_string =  output_string + 'Number of equivalent strain above ' + str(threshold_equiv_strain) +' threshold: ' + str(high_strain_counter) + ' which is ' + str(100*high_strain_counter/len(all_cells)) + ' percent \n'
    output_file.write('--------------------------------------------------------------------- \n')
    output_file.write(output_string)
    output_file.write('--------------------------------------------------------------------- \n')
output_file.close()
plt.plot(num_els,equivs_at_targets,'r-x')
