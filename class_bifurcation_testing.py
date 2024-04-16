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
tree_parameters = [gen_18, gen_19]#, gen_16, gen_17]
axial_rotations = [0, np.pi/2]

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
cont_outlet_vec = bifurcating_tree.cont_outlet_position
free_vertices_pcd = o3d.geometry.PointCloud()
free_vertices_pcd.points = o3d.utility.Vector3dVector(np.vstack([outlet_free_vertice, cont_outlet_vec]))
visualization = [combined_mesh, free_vertices_pcd]
combined_mesh.compute_triangle_normals()

o3d.visualization.draw_geometries(visualization, window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)

o3d.io.write_triangle_mesh("Second_Half_mesh.stl", combined_mesh)
