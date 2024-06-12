"""
Testing how many generations I can create from a .csv file
"""
import csv
import numpy as np
import open3d as o3d
from stylianou_geometry_equations_legacy import R_i
from single_pathway_bifurcating_model_class import single_pathway_bifurcating_model


with open("lung_dimensions_mm.csv") as fp:
    reader = csv.reader(fp, delimiter=",", quotechar='"')
    next(reader, None)  # skip the headers
    data_read = [row for row in reader]

lung_dimensions = np.array(data_read).astype(np.float64)
lung_dimensions[:,1:] = lung_dimensions[:,1:] * 1E3
n_gens = np.shape(lung_dimensions)[0]

print(lung_dimensions)

iota_b_test = 35*np.pi/180
iota_gamma_test = 39.83 * np.pi/180
delta_alpha_test = 3.6*np.pi/180

lung_parameters = [0] * (n_gens - 1)
axial_rotations = [0] * (n_gens - 1)

for j in range(1, n_gens):
    r_j_parent = lung_dimensions[j-1, 1]
    r_j_daughter = lung_dimensions[j,1]
    l_j_daughter = lung_dimensions[j,2]
    r_outer = r_j_daughter * 5
    gen_j_parameters = [r_j_parent, 0, r_j_daughter, l_j_daughter, r_outer, R_i,iota_b_test, delta_alpha_test, iota_gamma_test, 1]
    lung_parameters[j-1] = gen_j_parameters
    axial_rotations[j-1] = np.pi/2 * np.mod(j+1, 2)

print(lung_parameters)
# bifurcating_tree = single_pathway_bifurcating_model(lung_parameters[17:], axial_rotations[17:])
# combined_mesh = bifurcating_tree.tree_mesh

# visualization = [combined_mesh]
# o3d.visualization.draw_geometries(visualization, window_name = "Visualization",
#                                    mesh_show_wireframe = True, mesh_show_back_face = True)