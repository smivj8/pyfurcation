"""
Testing how many generations I can create from a .csv file
"""
import csv
import numpy as np
import open3d as o3d
from single_pathway_bifurcating_model_class import single_pathway_bifurcating_model



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
    parent_radius_j = lung_dimensions[j-1, 1]
    daughter_radius_j = lung_dimensions[j,1]
    daughter_length_j = lung_dimensions[j,2]
    outer_radius_j = daughter_radius_j * 5
    gen_j_parameters = [parent_radius_j, 0, daughter_radius_j, daughter_length_j, outer_radius_j ,branching_angle_test, carina_angle_test, 1]
    lung_parameters[j-1] = gen_j_parameters
    axial_rotations[j-1] = np.pi * np.mod(j-1, 2 ) + np.pi/2 + (np.random.rand()*np.pi/2 - (np.pi/4))

axial_rotations[0] = 0
n_min = 0

bifurcating_tree = single_pathway_bifurcating_model(lung_parameters[n_min:], axial_rotations[n_min:])
combined_mesh = bifurcating_tree.tree_mesh

visualization = [combined_mesh]
o3d.visualization.draw_geometries(visualization, window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)