import numpy as np
from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_mesh_class import bifurcation_mesh
from bifurcation_free_vertices_class import bifurcation_free_vertices
from stylianou_geometry_equations import R_i
from utility_functions_bifurcation import rotation_matrix_about_y
import open3d as o3d

D_0 = 100; R_0 = D_0/2; R_1 = D_0/2.48888; L_0 = 0; L_1 = 4.08888*R_1
R_o_0 = 5*R_1; iota_b_test = 35*np.pi/180; iota_gamma_test = 39.83 * np.pi/180
delta_alpha_test = 3.6*np.pi/180

D_1 = R_1*2; R_2 = D_1/2.48888; L_0 = 0; L_2 = 4.08888*R_2
R_o_1 = 5*R_2

parameters_0 = [R_0, 0, R_1, L_1, R_o_0, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
parameters_1 = [R_1, 0, R_2, L_2, R_o_1, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]

bifurcation_mesh_0 = bifurcation_mesh(parameters_0)
bifurcation_mesh_1 = bifurcation_mesh(parameters_1)

mesh_0 = bifurcation_mesh_0.truncated_mesh
mesh_1 = bifurcation_mesh_1.truncated_mesh

print("\n DONE GENERATING MESH\n")

#Testing Rotation Functionality:
free_vertices_0 = bifurcation_free_vertices(parameters_0, mesh_0)
inlet_free_vertices_0 = free_vertices_0.inlet_free_vertices
pos_z_outlet_free_vertices_0 = free_vertices_0.positive_outlet_free_vertices
neg_z_outlet_free_vertices_0 = free_vertices_0.negative_outlet_free_vertices

free_vertices_1 = bifurcation_free_vertices(parameters_1, mesh_1)
inlet_free_vertices_1 = free_vertices_1.inlet_free_vertices
pos_z_outlet_free_vertices_1 = free_vertices_1.positive_outlet_free_vertices
neg_z_outlet_free_vertices_1 = free_vertices_1.negative_outlet_free_vertices



vectors_0 = bifurcation_vectors(parameters_0)
r_pos_0 = vectors_0.r_pos
r_neg_0 = vectors_0.r_neg
n_pos_0 = vectors_0.n_pos
n_neg_0 = vectors_0.n_neg


matrix = rotation_matrix_about_y(np.pi)
vectors_rotation = bifurcation_vectors(parameters_0)
vectors_rotation.rotate_vector(matrix)
r_pos_1 = vectors_rotation.r_pos
r_neg_1 = vectors_rotation.r_neg
n_pos_1 = vectors_rotation.n_pos
n_neg_1 = vectors_rotation.n_neg



print("\nDEBUGGING STOP IF NEEDED\n")

#CONVERTING numpy xyz to open3d pointcloud
inlet_pointcloud_0 = o3d.geometry.PointCloud()
inlet_pointcloud_0.points = o3d.utility.Vector3dVector(inlet_free_vertices_0)
inlet_pointcloud_1 = o3d.geometry.PointCloud()
inlet_pointcloud_1.points = o3d.utility.Vector3dVector(inlet_free_vertices_1)

outlet_pos_pointcloud_0 = o3d.geometry.PointCloud()
outlet_pos_pointcloud_0.points = o3d.utility.Vector3dVector(pos_z_outlet_free_vertices_0)
outlet_pos_pointcloud_1 = o3d.geometry.PointCloud()
outlet_pos_pointcloud_1.points = o3d.utility.Vector3dVector(pos_z_outlet_free_vertices_1)

outlet_neg_pointcloud_0 = o3d.geometry.PointCloud()
outlet_neg_pointcloud_0.points = o3d.utility.Vector3dVector(neg_z_outlet_free_vertices_0)
outlet_neg_pointcloud_1 = o3d.geometry.PointCloud()
outlet_neg_pointcloud_1.points = o3d.utility.Vector3dVector(neg_z_outlet_free_vertices_1)

visualization = [mesh_0, mesh_1, inlet_pointcloud_0, inlet_pointcloud_1, outlet_pos_pointcloud_0, 
                 outlet_pos_pointcloud_1, outlet_neg_pointcloud_0, outlet_neg_pointcloud_1]

o3d.visualization.draw_geometries(visualization, window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)
