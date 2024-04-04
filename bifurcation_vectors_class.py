from stylianou_geometry_equations import *
import numpy as np
import open3d as o3d


class bifurcation_vectors:

    def __init__(self, parameters, is_normal_in):
        self.R_p = parameters[0]; self.L_p = parameters[1]; self.R_d = parameters[2]
        self.L_d = parameters[3]; self.R_o = parameters[4]; self.R_i = parameters[5]
        self.iota_b = parameters[6]; self.delta_alpha = parameters[7]
        self.iota_gamma = parameters[8]; self.n_cont_outlets = parameters[9]
        self.is_normal = is_normal_in; self.parameters = parameters
        # is_normal: True -> generate normal vector, False -> generate positional vector
        if self.is_normal:
            self.calculate_outlet_normals()
        else:
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
        r_D_pos = ((r0 + rpi)/2)
        
        r_D_neg = np.zeros(3)
        
        r_D_neg[0] = r_D_pos[0]
        r_D_neg[1] = r_D_pos[1]
        r_D_neg[2] = -1*r_D_pos[2]
        
        self.r_D_pos = r_D_pos
        self.r_D_neg = r_D_neg
        return self.r_D_pos, self.r_D_neg

    def calculate_outlet_normals(self):
        eta = np.tan(self.iota_gamma - self.iota_b)
        phi_s_max = self.iota_b + np.arctan(eta + (self.L_d/self.R_o))
        def phi_c_pos_min_lambda(phi_s, parameters, scalar):
            return -1*np.pi/2
        def phi_c_pos_max_lambda(phi_s, parameters, scalar):
            return 3*np.pi/2
        r_C, _ = r_positive(self.iota_b, 0, phi_c_pos_min_lambda, phi_c_pos_max_lambda, self.parameters, 1.1)
        r_F, _ = r_positive(phi_s_max, 0, phi_c_pos_min_lambda, phi_c_pos_max_lambda, self.parameters, 1.1)
        n_D_pos = ((r_F - r_C)/np.linalg.norm((r_F - r_C)))
        
        n_D_neg = np.zeros(3)
        
        n_D_neg[0] = n_D_pos[0]
        n_D_neg[1] = n_D_pos[1]
        n_D_neg[2] = -1*n_D_pos[2]
        
        self.n_D_pos = n_D_pos
        self.n_D_neg = n_D_neg
        
        return self.n_D_pos, self.n_D_neg
    
    def rotate(self, matrix):
        if self.is_normal:
            self.n_D_pos = matrix @ self.n_D_pos
            self.n_D_neg = matrix @ self.n_D_neg
            return self.n_D_pos, self.n_D_neg
        else:
            self.r_D_pos = matrix @ self.r_D_pos
            self.r_D_pos = matrix @ self.r_D_neg
            return self.r_D_pos, self.r_D_neg

    def translate(self, translation_vector):
        if self.is_normal:
            # I dont think I'm going to translate normal vectors for the purpose
            self.n_D_pos = self.n_D_pos + translation_vector
            self.n_D_neg = self.n_D_neg + translation_vector
            return self.n_D_pos, self.n_D_neg
        else:
            self.r_D_pos = self.r_D_pos + translation_vector
            self.r_D_neg = self.r_D_neg + translation_vector
            return self.r_D_pos, self.r_D_neg