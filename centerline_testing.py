"""
A script for testing the centerline functionality
"""
import csv
from time import time
import numpy as np
import open3d as o3d
from bifurcation_unit_centerline_class import bifurcation_unit_centerline

t0 = time()
#pull values from .csv file and convert to numpy array
with open("lung_dimensions_mm.csv") as fp:
    reader = csv.reader(fp, delimiter=",", quotechar='"')
    next(reader, None)  # skip the headers
    data_read = [row for row in reader]
lung_dimensions = np.array(data_read).astype(np.float64)
lung_dimensions[:,1:] = lung_dimensions[:,1:]
n_gens = np.shape(lung_dimensions)[0]
t1 = time()
print(f"Time to Read .csv file: {t1 - t0}")
print(f"Total Time: {t1-t0}\n")
#Get numpy parameter data into lists
branching_angle_test = 35*np.pi/180
carina_angle_test = 3.6*np.pi/180
lung_parameters = [0] * (n_gens - 1)
axial_rotations = [0] * (n_gens - 1)
for j in range(1, n_gens):
    parent_radius_j = lung_dimensions[j-1, 1]
    daughter_radius_j = lung_dimensions[j,1]
    daughter_length_j = lung_dimensions[j,2]
    outer_radius_j = daughter_radius_j * 5
    gen_j_parameters = [parent_radius_j, 0, daughter_radius_j, daughter_length_j, outer_radius_j ,branching_angle_test, carina_angle_test, 1]
    lung_parameters[j-1] = gen_j_parameters
    axial_rotations[j-1] = np.pi * np.mod(j-1, 2 ) - np.pi/2 + (np.random.rand()*np.pi/2 - (np.pi/4))
axial_rotations[0] = 0
t2 = time()
print(f"Time to construct lung dimensions list: {t2 - t1}")
print(f"Total Time: {t2 - t0}\n")
#test centerline function
centerline_module = bifurcation_unit_centerline(lung_parameters[0])
centerline_pcd = centerline_module.centerline_pointcloud
visualization = [centerline_pcd]
t3 = time()
print(f"Time to construct centerline pointcloud: {t3 - t2}")
print(f"Total Time: {t3 - t0}\n")

print(lung_dimensions[:, 0])

#o3d.visualization.draw_geometries(visualization, window_name = "Visualization")