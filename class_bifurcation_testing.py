import numpy as np
import open3d as o3d
from bifurcation_unit_modules.bifurcation_mesh_class import bifurcation_mesh
from bifurcation_tree_modules.single_pathway_bifurcating_model_class import single_pathway_bifurcating_model
# import utility_functions_bifurcation as ufb

test_branching_angle = 35*np.pi/180
test_carina_angle = 3.6*np.pi/180

#Create 4 generation tree from 15 to 18

gen_15 = [420, 2010, 335, 2010, 1675, test_branching_angle, test_carina_angle, 0]
gen_16 = [335, 0, 270, 1620, 1350, test_branching_angle, test_carina_angle, 1]
gen_17 = [270, 0, 215, 1290, 1075, test_branching_angle, test_carina_angle, 1]
gen_18 = [215, 0, 175, 1050,  875, test_branching_angle, test_carina_angle, 1]
tree_parameters = [gen_15, gen_16, gen_17, gen_18]
axial_rotations = [0, np.pi/2, 0, np.pi/2]

# bifurcating_tree = single_pathway_bifurcating_model(tree_parameters, axial_rotations)
# combined_mesh = bifurcating_tree.tree_mesh
mesh_class = bifurcation_mesh(gen_15)
combined_mesh = mesh_class.initial_capped_mesh
pcd = mesh_class.pointcloud

combined_mesh.remove_duplicated_vertices()
combined_mesh.remove_degenerate_triangles()

visualization = [combined_mesh, pcd]
combined_mesh.compute_triangle_normals()

o3d.visualization.draw_geometries(visualization, window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)

#o3d.io.write_triangle_mesh("Four_Generations_Mesh.stl", combined_mesh)
