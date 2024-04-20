import numpy as np
from utility_functions_bifurcation import *
from stylianou_geometry_equations import R_i


iota_b_test = 35*np.pi/180
iota_gamma_test = 39.83 * np.pi/180
delta_alpha_test = 3.6*np.pi/180


#Create 5 generation tree from 15 to 18

gen_15 = [420, 0, 335, 2010, 1675, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
gen_16 = [335, 0, 270, 1620, 1350, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
gen_17 = [270, 0, 215, 1290, 1075, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
gen_18 = [215, 0, 175, 1050,  875, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
tree_parameters = [gen_15, gen_16, gen_17, gen_18]
axial_rotations = [0, np.pi/2, 0, np.pi/2]

#get position and normal vectors of each bifurcation unit
r_pos = np.zeros((4,3))
r_neg = np.zeros((4,3))
n_pos = np.zeros((4,3))
n_neg = np.zeros((4,3))

for ind in range(4):
    r_pos[ind], r_neg[ind] = calculate_outlet_positions(tree_parameters[ind])
    n_pos[ind], n_neg[ind] = calculate_outlet_normals(tree_parameters[ind])

#get transformed position and normal vectors for the full tree
tree_positions_pos = np.zeros((4,3))
tree_positions_neg = np.zeros((4,3))
tree_normals_pos = np.zeros((4,3))
tree_normals_neg = np.zeros((4,3))

tree_positions_pos[0] = r_pos[0]
tree_positions_neg[0] = r_neg[0]
tree_normals_pos[0] = n_pos[0]
tree_normals_neg[0] = n_neg[0]

for j in range(1,4):
    R1 = unit_vector_match_rotation(np.array([1,0,0]), tree_normals_pos[j-1])
    R2 = rotation_about_vector_u(tree_normals_pos[j-1], axial_rotations[j])
    
    tree_positions_pos[j] = (R2 @ (R1 @ r_pos[j])) + tree_positions_pos[j-1]
    tree_positions_neg[j] = (R2 @ (R1 @ r_neg[j])) + tree_positions_pos[j-1]
    tree_normals_pos[j] =(R2 @ (R1 @ n_pos[j]))
    tree_normals_neg[j] =(R2 @ (R1 @ n_neg[j]))


#Get visualization planes
visualization_plane_normals = np.zeros((4,3))
visualization_plane_positions = np.zeros((4,3))

visualization_plane_normals[0] = np.cross(n_pos[0], n_neg[0])

for j in range(1,4):
    visualization_plane_normals[j] = np.cross(tree_normals_pos[j], tree_normals_neg[j])
    visualization_plane_positions[j] = tree_positions_pos[j-1]

#convert position vectors from micrometers to meters
tree_positions_pos = tree_positions_pos*1E-6
tree_positions_neg = tree_positions_neg*1E-6
visualization_plane_positions = visualization_plane_positions*1E-6

#Write it to a text file for easy copy/paste (csv would be easier but star doesn't make anything easy to copy and paste...)
output = open("Four_Generation_Plane_Values.txt", 'w')
for k in range(4):
    output.write(f"Generation {k+15} Plane Values\n\n\n")
    output.write( "Flow Plane Position = [\n")
    output.write(f"    {tree_positions_pos[k,0]}\n")
    output.write(f"    {tree_positions_pos[k,1]}\n")
    output.write(f"    {tree_positions_pos[k,2]}\n")
    output.write("]\n\n")
    output.write( "Flow Plane Normal = [\n")
    output.write(f"    {tree_normals_pos[k,0]}\n")
    output.write(f"    {tree_normals_pos[k,1]}\n")
    output.write(f"    {tree_normals_pos[k,2]}\n")
    output.write("]\n\n\n")
    output.write( "Visualization Plane Position = [\n")
    output.write(f"    {visualization_plane_positions[k,0]}\n")
    output.write(f"    {visualization_plane_positions[k,1]}\n")
    output.write(f"    {visualization_plane_positions[k,2]}\n")
    output.write("]\n\n")
    output.write( "Visualization Plane Normal = [\n")
    output.write(f"    {visualization_plane_normals[k,0]}\n")
    output.write(f"    {visualization_plane_normals[k,1]}\n")
    output.write(f"    {visualization_plane_normals[k,2]}\n")
    output.write("]\n\n")
    output.write("-----------------------------------------------------------------------------\n\n\n")

output.close()
