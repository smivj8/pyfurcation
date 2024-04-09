import numpy as np
from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_mesh_class import bifurcation_mesh
from bifurcation_free_vertices_class import bifurcation_free_vertices
from bifurcation_unit_class import bifurcation_unit
from stylianou_geometry_equations import R_i
from utility_functions_bifurcation import *
import open3d as o3d

D_0 = 100; R_0 = D_0/2; R_1 = D_0/2.48888; L_0 = 0; L_1 = 4.08888*R_1
R_o_0 = 5*R_1; iota_b_test = 35*np.pi/180; iota_gamma_test = 39.83 * np.pi/180
delta_alpha_test = 3.6*np.pi/180

D_1 = R_1*2; R_2 = D_1/2.48888; L_0 = 0; L_2 = 4.08888*R_2
R_o_1 = 5*R_2

parameters_0 = [R_0, 0, R_1, L_1, R_o_0, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
parameters_1 = [R_1, 0, R_2, L_2, R_o_1, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]


#Testing Rotation Functionality:
#Create position, normal, and location matrix
vectors_0 = bifurcation_vectors(parameters_0)
r_pos_0 = vectors_0.r_pos
n_pos_0 = vectors_0.n_pos
R1 = unit_vector_match_rotation(np.array([1,0,0]), n_pos_0)
R2 = rotation_about_vector_u(n_pos_0, np.pi/2)

bifurcation_unit_0 = bifurcation_unit(parameters_0)
bifurcation_unit_1 = bifurcation_unit(parameters_1)
bifurcation_unit_1.rotate_bifurcation_unit(R1)
bifurcation_unit_1.rotate_bifurcation_unit(R2)
bifurcation_unit_1.translate_bifurcation_unit(r_pos_0)

mesh_0 = bifurcation_unit_0.bifurcation_unit_mesh
mesh_1 = bifurcation_unit_1.bifurcation_unit_mesh

inlet_0 = bifurcation_unit_0.inlet_free_vertices
outlet_0 = bifurcation_unit_0.positive_outlet_free_vertices
inlet_1 = bifurcation_unit_1.inlet_free_vertices
outlet_1 = bifurcation_unit_1.positive_outlet_free_vertices

free_vertices = np.vstack([inlet_0, outlet_0, inlet_1, outlet_1])

free_vertices_pcd = o3d.geometry.PointCloud()
free_vertices_pcd.points = o3d.utility.Vector3dVector(free_vertices)


print("\n DONE GENERATING MESH\n")









print("\nDEBUGGING STOP IF NEEDED\n")



visualization = [mesh_0, mesh_1, free_vertices_pcd]

o3d.visualization.draw_geometries(visualization, window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)
