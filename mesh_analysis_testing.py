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

gen_16 = [420, 0, 335, 2010, 1675, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 0]

mesh_class = bifurcation_mesh(gen_16)
mesh = mesh_class.initial_capped_mesh
print("\n#######################   Done generating mesh   #########################\n")

analysis = mesh_class.perform_mesh_analysis()
volume = analysis[0]*((1E-6)**3)
surface_area = analysis[1]*((1E-6)**2)
print(f"Bifurcation Unit Volume: {volume}\n")
print(f"Bifurcation Unit Area:   {surface_area}\n")

mesh.compute_triangle_normals()
print(mesh.is_watertight())
#mesh.get_volume()
o3d.visualization.draw_geometries([mesh], window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)

o3d.io.write_triangle_mesh("Test_Mesh_For_Analysis.stl", mesh)
