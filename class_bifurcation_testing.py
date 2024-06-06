import numpy as np
import open3d as o3d
from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_mesh_class import bifurcation_mesh
from bifurcation_free_vertices_class import bifurcation_free_vertices
from single_pathway_bifurcating_model_class import single_pathway_bifurcating_model
from bifurcation_unit_class import bifurcation_unit
from stylianou_geometry_equations import R_i
from utility_functions_bifurcation import *


iota_b_test = 35*np.pi/180
iota_gamma_test = 39.83 * np.pi/180
delta_alpha_test = 3.6*np.pi/180

#Create 5 generation tree from 13 to 17

gen_16 = [420, 0, 335, 2010, 1675, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
gen_17 = [335, 0, 270, 1620, 1350, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
gen_18 = [270, 0, 215, 1290, 1075, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
gen_19 = [215, 0, 175, 1050,  875, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
tree_parameters = [gen_16, gen_17, gen_18, gen_19]
axial_rotations = [0, np.pi/2, 0, np.pi/2]

bifurcating_tree = single_pathway_bifurcating_model(tree_parameters, axial_rotations)
combined_mesh = bifurcating_tree.tree_mesh
acinus_location = bifurcating_tree.tree_outlet_positions[-1]
acinus_orientation = bifurcating_tree.tree_outlet_normals[-1]

print("\nDONE GENERATING MESH\n")

print("\nDEBUGGING STOP IF NEEDED\n")

combined_mesh.remove_duplicated_vertices()
outlet_free_vertice = bifurcating_tree.cont_outlet_vertices
free_vertices_pcd = o3d.geometry.PointCloud()
free_vertices_pcd.points = o3d.utility.Vector3dVector(np.vstack([outlet_free_vertice]))
visualization = [combined_mesh, free_vertices_pcd]
combined_mesh.compute_triangle_normals()

o3d.visualization.draw_geometries(visualization, window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)

#o3d.io.write_triangle_mesh("Four_Generations_Mesh.stl", combined_mesh)
print(f"Acinus Location:    {acinus_location}\n")
print(f"Acinus Orientation: {acinus_orientation}\n")
