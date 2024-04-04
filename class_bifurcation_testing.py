import numpy as np
from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_mesh_class import bifurcation_mesh
from stylianou_geometry_equations import R_i
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
mesh_0 = bifurcation_mesh_0.cropped_mesh
mesh_0 = bifurcation_mesh_1.cropped_mesh

mesh_0 = bifurcation_mesh_0.truncated_mesh

print("\n DONE GENERATING MESH\n")

inlet_free_vertices = bifurcation_mesh_0.get_free_edge_vertices(0)
pos_z_outlet_free_vertices = bifurcation_mesh_0.get_free_edge_vertices(1)
neg_z_outlet_free_vertices = bifurcation_mesh_0.get_free_edge_vertices(2)




#CONVERTING numpy xyz to open3d pointcloud
inlet_pointcloud = o3d.geometry.PointCloud()
inlet_pointcloud.points = o3d.utility.Vector3dVector(inlet_free_vertices)

outlet_pos_pointcloud = o3d.geometry.PointCloud()
outlet_pos_pointcloud.points = o3d.utility.Vector3dVector(pos_z_outlet_free_vertices)

outlet_neg_pointcloud = o3d.geometry.PointCloud()
outlet_neg_pointcloud.points = o3d.utility.Vector3dVector(neg_z_outlet_free_vertices)

o3d.visualization.draw_geometries([mesh_0, outlet_neg_pointcloud], window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)
