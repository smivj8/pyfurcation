"""
Function module for creating the pointcloud for a bifurcation unit using the functions
defined in stylianou_geometry_equations.py.

Construction of the bifurcation unit is dependent on a list of parameters. The parameters
argument completely defines the dimensions of a bifurcation unit. The parameters argument is
a list structured as follows: 

    parameters = [parent_radius, parent_length, daughter_radius, daughter_length, 
                  outer_radius, branching_angle, carina_angle, n_cont_outlets]

n_cont_outlets is the number of daughter outlets (0, 1, or 2) that will continue
on and attach to a successive bifurcation unit, and is not used in the creation of
the bifurcation unit pointcloud. As such, it is not accessed in this script.

The density of the pointcloud is determined by four additional parameters, defined
as follows:

    n_streamline (int, optional): number of points along length of bifurcation
        unit, default 120.

    n_circ (int, optional): number of points around the circumference of the 
        bifurcation unit, default 100.

    n_car_range (int, optional ): number of points along the arc of the carina
        connecting the daughter outlet pathways, default 80.

    n_car_rad (int, optional): number of points along the radius of the carina 
        connecting the daughter outlet pathways, default 8.

As stated, these parameters have default values for a high density pointcloud which have been
successfuly tested to produce appropriate meshes. These can be changed when calling the function.

"""
import numpy as np
import stylianou_geometry_equations as sge

def get_pipe_section_pointcloud_sizes(parameters, n_streamline, n_circ):
    """
    A method for determining the exact size of the pointcloud of the different sections
    (parent, branching) of the bifurcation unit. 

    Args:
        parameters (list): parameters of bifurcation unit, see docstring

        n_streamline (int, optional): number of points along length of bifurcation
            unit, default 120.

        n_circ (int, optional): number of points around the circumference of the 
            bifurcation unit, default 100.
    """
    #Get total range of phi_s
    phi_s_min = sge.get_phi_s_min(parameters)
    phi_s_max = sge.get_phi_s_max(parameters)
    phi_s_linspace = np.linspace(phi_s_min, phi_s_max, n_streamline)
    #determine (i.e. set) the number of points in the pos/neg branching pointclouds
    n_parent_streamline = 0
    while phi_s_linspace[n_parent_streamline] < 0:
        n_parent_streamline += 1
    parent_pointcloud_size = n_parent_streamline * n_circ
    #determine (i.e. set) the number of points in the pos/neg branching pointclouds
    branching_pointcloud_size = 0
    for j in range(n_parent_streamline, n_streamline):
        phi_c_pos_min = sge.get_phi_c_pos_min(phi_s_linspace[j], parameters)
        phi_c_pos_max = sge.get_phi_c_pos_max(phi_c_pos_min)
        phi_c_range = np.abs(phi_c_pos_max - phi_c_pos_min)
        n_circ_phi = int(np.round((phi_c_range/(2*np.pi))*n_circ))
        branching_pointcloud_size += n_circ_phi
    return parent_pointcloud_size, branching_pointcloud_size

def construct_pipe_section_pointcloud(parameters, n_streamline = 120, n_circ = 100):
    """
    A method to iterate through the phi_s and phi_c ranges to construct a high
    density pointcloud of the parent and branching pipe sections.

    Args:
        parameters (list): parameters of bifurcation unit, see docstring

        n_streamline (int, optional): number of points along length of bifurcation
            unit, default 120.

        n_circ (int, optional): number of points around the circumference of the 
            bifurcation unit, default 100.
    """
    #get the pointcloud section sizes 
    n_parent, n_branching= get_pipe_section_pointcloud_sizes(parameters, n_streamline, n_circ)
    #Initialize pointclouds
    parent_pointcloud = np.zeros((n_parent, 3))
    branching_pos_pointcloud = np.zeros((n_branching, 3))
    branching_neg_pointcloud = np.zeros((n_branching, 3))
    #get phi_s range
    phi_s_min = sge.get_phi_s_min(parameters)
    phi_s_max = sge.get_phi_s_max(parameters)
    phi_s_linspace = np.linspace(phi_s_min, phi_s_max, n_streamline)
    #Populate pipe (parent and branching) sections pointclouds
    pipe_section_index = 0
    for j in range(n_streamline):
        if phi_s_linspace[j] < 0 and n_parent > 0:   #populate the parent pipe section
            phi_c_pos_min = sge.get_phi_c_pos_min(phi_s_linspace[j], parameters)
            phi_c_pos_max = sge.get_phi_c_pos_max(phi_c_pos_min)
            phi_c_pos_linspace = np.linspace(phi_c_pos_min, phi_c_pos_max, n_circ)
            for k in range(n_circ):
                parent_pointcloud[pipe_section_index] = \
                    sge.get_r_parent(phi_s_linspace[j], phi_c_pos_linspace[k], parameters)
                pipe_section_index += 1
        else: #populate the positive and negative branching pipe sections!
            #positive range
            phi_c_pos_min = sge.get_phi_c_pos_min(phi_s_linspace[j], parameters)
            phi_c_pos_max = sge.get_phi_c_pos_max(phi_c_pos_min)
            phi_c_range = np.abs(phi_c_pos_max - phi_c_pos_min)
            n_circ_phi_s = int(np.round((phi_c_range/(2*np.pi))*n_circ))
            #negative range
            phi_c_neg_min = sge.get_phi_c_neg_min(phi_c_pos_min)
            phi_c_neg_max = sge.get_phi_c_neg_max(phi_c_pos_min)
            #create linspaces
            phi_c_pos_linspace = np.linspace(phi_c_pos_min, phi_c_pos_max, n_circ_phi_s)
            phi_c_neg_linspace = np.linspace(phi_c_neg_min, phi_c_neg_max, n_circ_phi_s)
            for k in range(n_circ_phi_s):
                branching_pos_pointcloud[pipe_section_index - n_parent] = \
                    sge.get_r_positive(phi_s_linspace[j], phi_c_pos_linspace[k], parameters)
                branching_neg_pointcloud[pipe_section_index - n_parent] = \
                    sge.get_r_negative(phi_s_linspace[j], phi_c_neg_linspace[k], parameters)
                pipe_section_index += 1
    #Combine pointclouds into single numpy array
    pipe_section_pointcloud = np.concatenate((parent_pointcloud, branching_pos_pointcloud, branching_neg_pointcloud))
    return pipe_section_pointcloud

def construct_carina_pointcloud(parameters, n_carina_range = 80, n_carina_rad = 8):
    """
    A method to iterate through the phi_s and psi ranges to construct a high
    density pointcloud of the carina. 

    Args:
        parameters (list): parameters of bifurcation unit, see docstring

        n_carina_range (int, optional ): number of points along the arc of the carina
            connecting the daughter outlet pathways, default 80.

        n_carina_rad (int, optional): number of points along the radius of the carina 
            connecting the daughter outlet pathways, default 8.
    """
    #Initialize carina pointcloud
    carina_pos_pointcloud = np.zeros((n_carina_range*n_carina_rad, 3))
    carina_neg_pointcloud = np.zeros((n_carina_range*n_carina_rad, 3))
    #populate carina pointcloud
    angle_phi_gamma = sge.get_angle_phi_gamma(parameters)
    phi_s_carina_linspace = np.linspace(0, angle_phi_gamma, n_carina_range)
    carina_index = 0
    for m in range(n_carina_range):
        _, carina_angle, _, _ = sge.get_local_carina_radius_and_angles(phi_s_carina_linspace[m], parameters)
        psi_linspace = np.linspace(-1*carina_angle, carina_angle, n_carina_rad)
        for n in range(n_carina_rad):
            carina_pos_pointcloud[carina_index], carina_neg_pointcloud[carina_index] = \
                sge.get_carina_position(phi_s_carina_linspace[m], psi_linspace[n], parameters)
            carina_index += 1
    #combine pointclouds into single numpy array
    carina_pointcloud = np.concatenate((carina_pos_pointcloud, carina_neg_pointcloud))
    return carina_pointcloud
    
################################################################################################
#                                        TESTING                                               #
################################################################################################
test_branching_angle = 35*np.pi/180
test_carina_angle = 3.6*np.pi/180
test_parameters = gen_15 = [420, 0.05 * 2010, 335, 2010, 1675, test_branching_angle, test_carina_angle, 1]
pipe_section_pcd = construct_pipe_section_pointcloud(test_parameters)
#carina_pcd = construct_carina_pointcloud(test_parameters)

import open3d as o3d
pointcloud = o3d.geometry.PointCloud()
pointcloud.points = o3d.utility.Vector3dVector(pipe_section_pcd)#np.concatenate(pipe_section_pcd, carina_pcd))
o3d.visualization.draw_geometries([pointcloud], window_name = "Visualization")