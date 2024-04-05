from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_mesh_class import bifurcation_mesh
from bifurcation_free_vertices_class import bifurcation_free_vertices
import numpy as np
import open3d as o3d

class bifurcation_unit:

    def __init__(self, parameters):
        self.R_p = parameters[0]; self.L_p = parameters[1]; self.R_d = parameters[2]
        self.L_d = parameters[3]; self.R_o = parameters[4]; self.R_i = parameters[5]
        self.iota_b = parameters[6]; self.delta_alpha = parameters[7]
        self.iota_gamma = parameters[8]; self.n_cont_outlets = parameters[9]
        self.parameters = parameters
        self.r_pos, self.r_neg = bifurcation_vectors(self.parameters, False).calculate_outlet_positions()
        self.n_pos, self.n_neg = bifurcation_vectors(self.parameters, True).calculate_outlet_normals()

