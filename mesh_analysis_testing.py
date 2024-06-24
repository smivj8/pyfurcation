import csv
import numpy as np
import open3d as o3d
from bifurcation_mesh_class import bifurcation_mesh
from utility_functions_bifurcation import *

iota_b_test = 35*np.pi/180
delta_alpha_test = 3.6*np.pi/180

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
    parent_radius_j = lung_dimensions[j-1, 1] * 1000
    daughter_radius_j = lung_dimensions[j,1] * 1000
    daughter_length_j = lung_dimensions[j,2] * 1000
    outer_radius_j = daughter_radius_j * 5
    gen_j_parameters = [parent_radius_j, 0, daughter_radius_j, daughter_length_j, outer_radius_j ,branching_angle_test, carina_angle_test, 0]
    lung_parameters[j-1] = gen_j_parameters
    axial_rotations[j-1] = np.pi * np.mod(j-1, 2 ) - np.pi/2 + (np.random.rand()*np.pi/2 - (np.pi/4))

mesh_volumes = open("bifurcation_unit_volume_and_area.txt", "w")
for j, unit in enumerate(lung_parameters):
    mesh_class = bifurcation_mesh(unit)
    mesh = mesh_class.initial_capped_mesh
    print(f"\n#######################   Done generating mesh Gen {j} #########################\n")
    # o3d.visualization.draw_geometries([mesh], window_name = "Visualization",
    #                                     mesh_show_wireframe = True, mesh_show_back_face = True)
    # print(mesh.is_watertight())
    analysis = mesh_class.perform_mesh_analysis()
    volume = analysis[0]*((1E-6)**3)
    surface_area = analysis[1]*((1E-6)**2)
    mesh_volumes.write(f"Volume Gen {j}: {volume}\n")
    mesh_volumes.write(f"SurfArea Gen {j}: {surface_area}\n\n")
    print(volume)
    print(surface_area)

mesh_volumes.close()
#mesh.get_volume()
