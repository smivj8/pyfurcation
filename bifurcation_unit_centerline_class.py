import numpy as np
import open3d as o3d
from bifurcation_vectors_class import bifurcation_vectors

class bifurcation_unit_centerline:

    def __init__(self, parameters, centerline_density = 100):
        """
        A class for the generation of an APPROXIMATE centerline diagram of a bifurcation
        unit. To be utilized in quickly generating the centerline diagram of entire tree
        without generating an entire tree mesh. 

        [1] T. Heistracher and W. Hofmann, “Physiologically realistic models of bronchial airway 
        bifurcations,” Journal of Aerosol Science, vol. 26, no. 3, pp. 497-509, Apr. 1995, 
        doi: 10.1016/0021-8502(94)00113-D.

        [2] F. S. Stylianou, J. Sznitman, and S. C. Kassinos, “Direct numerical simulation of
        particle laden flow in a human airway bifurcation model,” International Journal of Heat and 
        Fluid Flow, vol. 61, pp. 677-710, Oct. 2016, doi: 10.1016/j.ijheatfluidflow.2016.07.013.

        Args:
            parameters (list): bifurcation unit parameters, see github README for more information
        
        """
        self.parameters = parameters
        self.parent_radius = parameters[0]
        self.daughter_radius = parameters[2]
        self.daughter_length = parameters[3]
        self.outer_radius = parameters[4]
        self.branching_angle = parameters[5]
        self.carina_angle = parameters[6]
        self.n_cont_outlets = parameters[7]
        self.centerline_density = centerline_density
        self.generate_unit_centerline()

    def generate_unit_centerline(self):
        """
        A method to generate the centerline diagram. Returns the xyz coordinates of the centerlines
        as both a numpy array and an open3d pointcloud object.
        """
        #Get position and normal vectors of the bifurcation unit
        vectors = bifurcation_vectors(self.parameters)
        outlet_pos = vectors.r_pos
        outlet_neg = vectors.r_neg
        normal_pos = vectors.n_pos
        # Initialize the origin vector
        origin = np.zeros(3)
        # dtermine the location of the bifurcation position, the vector that connects the 
        # x_axis to the positive outlet, along the vector normal to the outlet
        k =  outlet_pos[2]/normal_pos[2]
        bifurcation_x = outlet_pos[0] - k*(normal_pos[0])
        bifurcation_position = np.array([bifurcation_x, 0, 0])
        #Construct the origin to bifurcation centerline
        n_centerline_parent = int(np.round(bifurcation_x/(k + bifurcation_x) * self.centerline_density))
        t_parameter = np.linspace(0, 1, n_centerline_parent, endpoint = False)
        parent_centerline = np.zeros((n_centerline_parent, 3))
        for ind, t in enumerate(t_parameter):
            parent_centerline[ind] = origin + t * (bifurcation_position - origin)
        # Construct the daughter centerlines
        n_centerline_daughter = self.centerline_density - n_centerline_parent
        s_parameter = np.linspace(0, 1, n_centerline_daughter)
        positive_daughter_centerline = np.zeros((n_centerline_daughter, 3))
        negative_daughter_centerline = np.zeros((n_centerline_daughter, 3))
        for ind, s in enumerate(s_parameter):
            positive_daughter_centerline[ind] = bifurcation_position + s * (outlet_pos - bifurcation_position)
            negative_daughter_centerline[ind] = bifurcation_position + s * (outlet_neg - bifurcation_position)
        self.centerline_array = np.vstack([parent_centerline, positive_daughter_centerline, negative_daughter_centerline])
        self.centerline_pointcloud = o3d.geometry.PointCloud()
        self.centerline_pointcloud.points = o3d.utility.Vector3dVector(self.centerline_array)
    
    # def rotate_
