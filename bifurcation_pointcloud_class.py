import numpy as np
import stylianou_geometry_equations as sge
class bifurcation_xyz:

    def __init__(self, parameters):
        """
        A class for the generation of the pointcloud for a SYMMETRIC bifurcation unit
        Generates high density xyz point cloud using bifurcation unit parameters from [1] and [2].

        [1] T. Heistracher and W. Hofmann, “Physiologically realistic models of bronchial airway 
        bifurcations,” Journal of Aerosol Science, vol. 26, no. 3, pp. 497–509, Apr. 1995, 
        doi: 10.1016/0021-8502(94)00113-D.

        [2] F. S. Stylianou, J. Sznitman, and S. C. Kassinos, “Direct numerical simulation of 
        particle laden flow in a human airway bifurcation model,” International Journal of Heat 
        and Fluid Flow, vol. 61, pp. 677–710, Oct. 2016, doi: 10.1016/j.ijheatfluidflow.2016.07.013.

        Args:
            parameters (list): bifurcation unit parameters, see github README for more information
        """
        self.parameters = parameters
        self.parent_radius = parameters[0]
        self.parent_length = parameters[1]
        self.daughter_radius = parameters[2]
        self.daughter_length = parameters[3]
        self.outer_radius = parameters[4]
        self.branching_angle = parameters[5]
        self.carina_angle = parameters[6]
        self.n_cont_outlets = parameters[7]

    



