from stylianou_geometry_equations import *
import numpy as np
import open3d as o3d


class bifurcation_vectors:

    def __init__(self, parameters):
        self.R_p = parameters[0]; self.L_p = parameters[1]; self.R_d = parameters[2]
        self.L_d = parameters[3]; self.R_o = parameters[4]; self.R_i = parameters[5]
        self.iota_b = parameters[6]; self.delta_alpha = parameters[7]
        self.iota_gamma = parameters[8]; self.n_cont_outlets = parameters[9]
        self.parameters = parameters
        # is_normal: True -> generate normal vector, False -> generate positional vector
        self.calculate_outlet_normals()
        self.calculate_outlet_positions()

    def calculate_outlet_positions(self):
        eta = np.tan(self.iota_gamma - self.iota_b)
        phi_s_max = self.iota_b + np.arctan(eta + (self.L_d/self.R_o))
        def phi_c_pos_min_lambda(phi_s, parameters, scalar):
            return -1*np.pi/2
        def phi_c_pos_max_lambda(phi_s, parameters, scalar):
            return 3*np.pi/2
        r0, _ =  r_positive(phi_s_max, 0, phi_c_pos_min_lambda, phi_c_pos_max_lambda, self.parameters, 1.1)
        rpi, _ = r_positive(phi_s_max, np.pi, phi_c_pos_min_lambda, phi_c_pos_max_lambda, self.parameters, 1.1)
        r_pos = ((r0 + rpi)/2)
        
        r_neg = np.zeros(3)
        
        r_neg[0] = r_pos[0]
        r_neg[1] = r_pos[1]
        r_neg[2] = -1*r_pos[2]
        
        self.r_pos = r_pos
        self.r_neg = r_neg
        return self.r_pos, self.r_neg

    def calculate_outlet_normals(self):
        eta = np.tan(self.iota_gamma - self.iota_b)
        phi_s_max = self.iota_b + np.arctan(eta + (self.L_d/self.R_o))
        def phi_c_pos_min_lambda(phi_s, parameters, scalar):
            return -1*np.pi/2
        def phi_c_pos_max_lambda(phi_s, parameters, scalar):
            return 3*np.pi/2
        r_C, _ = r_positive(self.iota_b, 0, phi_c_pos_min_lambda, phi_c_pos_max_lambda, self.parameters, 1.1)
        r_F, _ = r_positive(phi_s_max, 0, phi_c_pos_min_lambda, phi_c_pos_max_lambda, self.parameters, 1.1)
        n_pos = ((r_F - r_C)/np.linalg.norm((r_F - r_C)))
        
        n_neg = np.zeros(3)
        
        n_neg[0] = n_pos[0]
        n_neg[1] = n_pos[1]
        n_neg[2] = -1*n_pos[2]
        
        self.n_pos = n_pos
        self.n_neg = n_neg
        
        return 
    
    def rotate_vector(self, matrix):
        self.n_pos = matrix @ self.n_pos
        self.n_neg = matrix @ self.n_neg
        self.r_pos = matrix @ self.r_pos
        self.r_neg = matrix @ self.r_neg

    def translate_vector(self, translation_vector):
        self.n_pos = self.n_pos + translation_vector
        self.n_neg = self.n_neg + translation_vector
        self.r_pos = self.r_pos + translation_vector
        self.r_neg = self.r_neg + translation_vector