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

def get_pipe_section_pointcloud_size(parameters, n_streamline, n_circ):
    """
    A method for determining the exact size of the pointcloud of the pipe section
    of the bifurcation unit. 

    Args:
        parameters (list): parameters of bifurcation unit, see docstring

        n_streamline (int, optional): number of points along length of bifurcation
            unit, default 120.

        n_circ (int, optional): number of points around the circumference of the 
            bifurcation unit, default 100.
    """
    #Unpack parameters list for any required parameters in this function
    parent_radius = parameters[0]
    daughter_radius = parameters[2]
    branching_angle = parameters[5]
    #Get total range of phi_s
    phi_s_min = sge.get_phi_s_min(parameters)
    phi_s_max = sge.get_phi_s_max(parameters)
    phi_s_range = np.abs((phi_s_max - phi_s_min))
    #determine number of points per section by divying up the total number of 
    #points into each section
    n_parent_streamline = np.max([2, (np.abs(phi_s_min)/phi_s_range * n_streamline)])
    n_merging_streamline = branching_angle/phi_s_range * n_streamline

    return

def get_total_pointcloud_size(parameters, n_circ, n_streamline,
                        n_car_range, n_car_rad):
    """
    A method for determining the exact size of the pointcloud for a bifurcation unit using 
    default number a points along each range. 

    Args:
        parameters (list): parameters of bifurcation unit, see docstring

        n_streamline (int, optional): number of points along length of bifurcation
            unit, default 120.

        n_circ (int, optional): number of points around the circumference of the 
            bifurcation unit, default 100.

        n_car_range (int, optional ): number of points along the arc of the carina
            connecting the daughter outlet pathways, default 80.

        n_car_rad (int, optional): number of points along the radius of the carina 
            connecting the daughter outlet pathways, default 8.
    """
    return

def construct_total_pointcloud(parameters, n_circ = 100, n_streamline = 120,
                        n_car_range = 80, n_car_rad = 8):
    return
