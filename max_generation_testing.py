"""
Testing how many generations I can create from a .csv file
"""
import csv
import numpy as np
import open3d as o3d
from bifurcation_tree_modules.single_pathway_bifurcating_model_class import single_pathway_bifurcating_model



with open("lung_dimensions_mm.csv") as fp:
    reader = csv.reader(fp, delimiter=",", quotechar='"')
    next(reader, None)  # skip the headers
    data_read = [row for row in reader]

lung_dimensions = np.array(data_read).astype(np.float64)
lung_dimensions[:,1:] = lung_dimensions[:,1:]
n_gens = np.shape(lung_dimensions)[0]


branching_angle_test = 35*np.pi/180
carina_angle_test = 3.6*np.pi/180

lung_parameters = [0] * (n_gens - 1)
axial_rotations = [0] * (n_gens - 1)
# axial_rotations = np.random.rand(n_gens - 1) * np.pi * 2
# 

for j in range(1, n_gens):
    parent_radius_j = lung_dimensions[j-1, 1] * 1000
    daughter_radius_j = lung_dimensions[j,1] * 1000
    daughter_length_j = lung_dimensions[j,2] * 1000
    outer_radius_j = daughter_radius_j * 5
    gen_j_parameters = [parent_radius_j, 0, daughter_radius_j, daughter_length_j, outer_radius_j ,branching_angle_test, carina_angle_test, 1]
    lung_parameters[j-1] = gen_j_parameters
    axial_rotations[j-1] = np.pi * np.mod(j-1, 2 ) - np.pi/2 + (np.random.rand()*np.pi/2 - (np.pi/4))

axial_rotations[0] = 0
n_min = 0

bifurcating_tree = single_pathway_bifurcating_model(lung_parameters[n_min:], axial_rotations[n_min:])
combined_mesh = bifurcating_tree.tree_mesh

outlet_positions = bifurcating_tree.tree_outlet_positions
negative_positions = bifurcating_tree.negative_outlet_positions
outlet_normals = bifurcating_tree.tree_outlet_normals
negative_normals = bifurcating_tree.negative_outlet_normals

visualization = [combined_mesh]
combined_mesh.compute_triangle_normals()

o3d.visualization.draw_geometries(visualization, window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)

o3d.io.write_triangle_mesh("Bifurcating_Tree_Geometry.stl", combined_mesh)

bifurcation_output = open("bifurcation_output.txt", "w")
bifurcation_output.write(f"Generated Unilateral Bifurcating Pathway with {n_gens} Generations\n\n")
bifurcation_output.write("#############################################################################\n\n\n")
bifurcation_output.write("AXIAL ROTATIONS = \n\t[")
for i in range(n_gens - 1):
    bifurcation_output.write(f"{axial_rotations[i]}, ")
bifurcation_output.write("\b\b]\n\n\n")
bifurcation_output.write("CONTINUATION POSITIONS = \n")
for ind, pos in enumerate(outlet_positions):
    bifurcation_output.write(f"[{pos[0]/1000000}, {pos[1]/1000000}, {pos[2]/1000000}]\n")
bifurcation_output.write("\n\nCONTINUATION NORMALS = \n")
for ind, norm in enumerate(outlet_normals):
    bifurcation_output.write(f"[{norm[0]/1000000}, {norm[1]/1000000}, {norm[2]/1000000}]\n")
bifurcation_output.write("\n\n#####################################################################################\n\n\n")
bifurcation_output.write("OUTLET POSITIONS = \n")
for ind, pos in enumerate(negative_positions):
    bifurcation_output.write(f"[{pos[0]/1000000}, {pos[1]/1000000}, {pos[2]/1000000}]\n")
bifurcation_output.write("\n\nCONTINUATION NORMALS = \n")
for ind, norm in enumerate(negative_normals):
    bifurcation_output.write(f"[{norm[0]/1000000}, {norm[1]/1000000}, {norm[2]/1000000}]\n")
bifurcation_output.write("\n\n#####################################################################################\n\n\n")
bifurcation_output.write("VISUALIZATION PLANE NORMALS = \n")
for j in range(n_gens - 1):
    vis_norm = np.cross(outlet_normals[j], negative_normals[j])
    bifurcation_output.write(f"[{vis_norm[0]/1000000}, {vis_norm[1]/1000000}, {vis_norm[2]/1000000}]\n")
