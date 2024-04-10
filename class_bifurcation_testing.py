import numpy as np
import open3d as o3d
from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_mesh_class import bifurcation_mesh
from bifurcation_free_vertices_class import bifurcation_free_vertices
from single_pathway_bifurcating_model_class import single_pathway_bifurcating_model
from bifurcation_unit_class import bifurcation_unit
from stylianou_geometry_equations import R_i
from utility_functions_bifurcation import *


D_0 = 100; R_0 = D_0/2; R_1 = D_0/2.48888; L_0 = 0; L_1 = 4.08888*R_1
R_o_0 = 5*R_1; iota_b_test = 35*np.pi/180; iota_gamma_test = 39.83 * np.pi/180
delta_alpha_test = 3.6*np.pi/180

D_1 = R_1*2; R_2 = D_1/2.48888; L_0 = 0; L_2 = 4.08888*R_2
R_o_1 = 5*R_2

test_parameters_0 = [R_0, 0, R_1, L_1, R_o_0, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
test_parameters_1 = [R_1, 0, R_2, L_2, R_o_1, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
test_tree_parameters = [test_parameters_0, test_parameters_1]
#Create 5 generation tree from 13 to 17
gen_13 = [475, 0, 410, 2700, 2050, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
gen_14 = [410, 0, 370, 1600, 1850, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
gen_15 = [370, 0, 250, 1330, 1250, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
gen_16 = [250, 0, 245, 1120, 1225, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
gen_17 = [245, 0, 200,  930, 1000, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
tree_parameters = [gen_13, gen_14, gen_15, gen_16, gen_17]
axial_rotations = [0, np.pi/2, 0, np.pi/2, 0]

bifurcating_tree = single_pathway_bifurcating_model(tree_parameters, axial_rotations)
combined_mesh = bifurcating_tree.tree_mesh

#Testing Rotation Functionality:
#Create position, normal, and location matrix
# vectors_0 = bifurcation_vectors(parameters_0)
# r_pos_0 = vectors_0.r_pos
# n_pos_0 = vectors_0.n_pos
# R1 = unit_vector_match_rotation(np.array([1,0,0]), n_pos_0)
# R2 = rotation_about_vector_u(n_pos_0, np.pi/2)

# bifurcation_unit_0 = bifurcation_unit(parameters_0)
# bifurcation_unit_1 = bifurcation_unit(parameters_1)
# bifurcation_unit_1.rotate_bifurcation_unit(R1)
# bifurcation_unit_1.rotate_bifurcation_unit(R2)
# bifurcation_unit_1.translate_bifurcation_unit(r_pos_0)

# mesh_0 = bifurcation_unit_0.bifurcation_unit_mesh
# mesh_1 = bifurcation_unit_1.bifurcation_unit_mesh

# inlet_0 = bifurcation_unit_0.inlet_free_vertices
# outlet_0 = bifurcation_unit_0.positive_outlet_free_vertices
# inlet_1 = bifurcation_unit_1.inlet_free_vertices
# outlet_1 = bifurcation_unit_1.positive_outlet_free_vertices

# free_vertices = np.vstack([inlet_0, outlet_0, inlet_1, outlet_1])

# free_vertices_pcd = o3d.geometry.PointCloud()
# free_vertices_pcd.points = o3d.utility.Vector3dVector(free_vertices)

# junction_mesh = create_junction_mesh(outlet_0, inlet_1, n_pos_0)


print("\nDONE GENERATING MESH\n")

print("\nDEBUGGING STOP IF NEEDED\n")

combined_mesh.remove_duplicated_vertices()
outlet_free_vertice = bifurcating_tree.cont_outlet_vertices
free_vertices_pcd = o3d.geometry.PointCloud()
free_vertices_pcd.points = o3d.utility.Vector3dVector(outlet_free_vertice)
visualization = [combined_mesh, free_vertices_pcd]


o3d.visualization.draw_geometries(visualization, window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)
