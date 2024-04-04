"""
Samuel Smith
MATH6019 - Nonlinear Programming
Date: Thu Feb  1 18:28:00 2024

Assignment
Problem

Description:

Author: smit3sl@mail.uc.edu
"""
import numpy as np
#import matplotlib.pyplot as plt
#import scipy.spatial as spatial
#from stl import mesh
#from matplotlib import cm
#from mpl_toolkits.mplot3d import Axes3D
#import surf2stl
import open3d as o3d


#My Parameters List:
    
# R_p = parameters[0]; L_p = parameters[1]; R_d = parameters[2]
# L_d = parameters[3]; R_o = parameters[4]; R_i = parameters[5]
# iota_b = parameters[6]; delta_alpha = parameters[7]
# iota_gamma = parameters[8]


def parent_pipe_section(phi_s, phi_c, R_p, L_p, R_o, scalar):
    xi = L_p + R_o * np.tan(phi_s)
    x_p = xi
    y_p = R_p * np.cos(phi_c)
    z_p = R_p * np.sin(phi_c)
    xyz_p = np.array([x_p, y_p, z_p])
    
    parent_dummy = np.array([x_p, scalar*y_p, scalar*z_p])
    normals_p = parent_dummy - xyz_p
    
    return xyz_p, normals_p * (1/np.linalg.norm(normals_p))

def R_i(phi_s, R_p, R_d, iota_b):
    
    def sigmoidal(phi_s, iota_b):
        if phi_s < 0:
            value = 0
        elif phi_s >= 0 and phi_s <= iota_b:
            value = 3 * ((phi_s/iota_b)**2) - 2 * ((phi_s/iota_b)**3)
        else:
            value = 1
        return value
    
    R_i_out = R_p - (R_p - R_d)*sigmoidal(phi_s, iota_b)
    
    return R_i_out


def bending_pipe_positive_z(phi_s, phi_c, L_p, R_o, R_i, R_p, R_d, iota_b, scalar):
    x_m_pos = L_p + R_o * np.sin(phi_s) - np.sin(phi_s) * R_i(phi_s, R_p, R_d, iota_b) * np.sin(phi_c)
    y_m_pos = R_i(phi_s, R_p, R_d, iota_b) * np.cos(phi_c)
    z_m_pos = R_o * (1-np.cos(phi_s)) + np.cos(phi_s)*R_i(phi_s, R_p, R_d, iota_b)*np.sin(phi_c)
    xyz_m = np.array([x_m_pos, y_m_pos, z_m_pos])
    
    dummy_x_m_pos = L_p + R_o * np.sin(phi_s) - np.sin(phi_s) * R_i(phi_s, R_p*scalar, R_d*scalar, iota_b) * np.sin(phi_c)
    dummy_y_m_pos = R_i(phi_s, R_p*scalar, R_d*scalar, iota_b) * np.cos(phi_c)
    dummy_z_m_pos = R_o * (1-np.cos(phi_s)) + np.cos(phi_s)*R_i(phi_s, R_p*scalar, R_d*scalar, iota_b)*np.sin(phi_c)
    
    bending_dummy = np.array([dummy_x_m_pos, dummy_y_m_pos, dummy_z_m_pos])
    normals_b = bending_dummy - xyz_m

    return xyz_m, normals_b * (1/np.linalg.norm(normals_b))

def daughter_positive_z(phi_s, phi_c, iota_b, L_p, R_i, R_o, R_p, R_d, scalar):
    xyz_m_pos, _ = bending_pipe_positive_z(iota_b, phi_c, L_p, R_o, R_i, R_p, R_d, iota_b, scalar)
    
    zeta = R_o * np.tan(phi_s - iota_b)
    x_d_pos = xyz_m_pos[0] + zeta*np.cos(iota_b)
    y_d_pos = xyz_m_pos[1]
    z_d_pos = xyz_m_pos[2] + zeta*np.sin(iota_b)
    
    xyz_d = np.array([x_d_pos, y_d_pos, z_d_pos])
    
    xyz_m_dummy, _ = bending_pipe_positive_z(iota_b, phi_c, L_p, R_o, R_i, R_p*scalar, R_d*scalar, iota_b, scalar)
    dummy_x_d_pos = xyz_m_dummy[0] + zeta*np.cos(iota_b)
    dummy_y_d_pos = xyz_m_dummy[1]
    dummy_z_d_pos = xyz_m_dummy[2] + zeta*np.sin(iota_b)
    
    daughter_dummy = np.array([dummy_x_d_pos, dummy_y_d_pos, dummy_z_d_pos])
    normals_d = daughter_dummy - xyz_d
   
    return xyz_d, normals_d*(1/np.linalg.norm(normals_d))

def r_positive(phi_s, phi_c, phi_c_pos_min, phi_c_pos_max, parameters, scalar):
    R_p = parameters[0]; L_p = parameters[1]; R_d = parameters[2]
    L_d = parameters[3]; R_o = parameters[4]; R_i = parameters[5]
    iota_b = parameters[6]; delta_alpha = parameters[7]
    iota_gamma = parameters[8]
    
    phi_c_min = phi_c_pos_min(phi_s, parameters, scalar)
    phi_c_max = phi_c_pos_max(phi_s, parameters, scalar)
    
    
    if phi_s < 0 and (phi_c >= (-1*np.pi/2) and phi_c <= 3*np.pi/2):
        xyz, normals = parent_pipe_section(phi_s, phi_c, R_p, L_p, R_o, scalar)
        return xyz, normals
    elif (phi_s >= 0 and phi_s < iota_b) and (phi_c >= phi_c_min and phi_c <= phi_c_max):
        xyz, normals = bending_pipe_positive_z(phi_s, phi_c, L_p, R_o, R_i, R_p, R_d, iota_b, scalar)
        return xyz, normals
    elif phi_s >= iota_b and (phi_c >= phi_c_min and phi_c <= phi_c_max):
        xyz, normals = daughter_positive_z(phi_s, phi_c, iota_b, L_p, R_i, R_o, R_p, R_d, scalar)
        return xyz, normals
    else:
        return None, None
    
def position_A(parameters, delta_alpha, scalar):
    def phi_c_pos_min_lambda(phi_s, parameters, scalar):
        return -1*np.pi/2
    def phi_c_pos_max_lambda(phi_s, parameters, scalar):
        return 3*np.pi/2
    
    
    xyz_A, _ = r_positive(0, delta_alpha/2, phi_c_pos_min_lambda, phi_c_pos_max_lambda, parameters, scalar)
    return xyz_A

def position_gamma(parameters, iota_gamma, scalar):
    def phi_c_pos_min_lambda(phi_s, parameters, scalar):
        return -1*np.pi/2
    def phi_c_pos_max_lambda(phi_s, parameters, scalar):
        return 3*np.pi/2
    
    xyz_gamma, _ = r_positive(iota_gamma, -1*np.pi/2, phi_c_pos_min_lambda, phi_c_pos_max_lambda, parameters, scalar)
    return xyz_gamma

def phi_c_pos_min(phi_s, parameters, scalar):
    R_p = parameters[0]; L_p = parameters[1]; R_d = parameters[2]
    R_o = parameters[4]; R_i = parameters[5]
    iota_b = parameters[6]; delta_alpha = parameters[7]
    iota_gamma = parameters[8]
    
    xyz_A = position_A(parameters, delta_alpha, scalar)
    xyz_gamma = position_gamma(parameters, iota_gamma, scalar)
    
    a = (xyz_gamma[2] - xyz_A[2]) / (xyz_gamma[0] - xyz_A[0])
    b = xyz_A[2] - ((xyz_gamma[2] - xyz_A[2]) / (xyz_gamma[0] - xyz_A[0])) * xyz_A[0]
    
    if phi_s >= 0 and phi_s <= iota_b:
        numerator = b + (a*(L_p + R_o *np.sin(phi_s))) - (R_o * (1-np.cos(phi_s)))
        denominator = (np.cos(phi_s) + a * np.sin(phi_s)) * R_i(phi_s, R_p, R_d, iota_b)
        return np.arcsin(numerator/denominator)
    elif phi_s > iota_b and phi_s <= iota_gamma:
        numerator_1 =  b + (a*(L_p + R_o *np.sin(iota_b))) - (R_o * (1-np.cos(iota_b)))
        denominator = (np.cos(iota_b) + a * np.sin(iota_b)) * R_i(iota_b, R_p, R_d, iota_b)
        
        numerator_2 = (a * np.cos(iota_b) - np.sin(iota_b)) * R_o * np.tan(phi_s - iota_b)
        # TODO this is where the error is being thrown and "NoneType" is being returned? 
        # but tbh that does not make any sense
        return np.arcsin((numerator_1/denominator) + (numerator_2/denominator))
    elif phi_s > iota_gamma:
        return -1*np.pi/2
    else:
        return -1*np.pi/2
    
def phi_c_pos_max(phi_s, parameters, scalar):
    angle = np.pi - phi_c_pos_min(phi_s, parameters, scalar)
    return angle

def phi_c_neg_min(phi_s, parameters, scalar):
    angle = phi_c_pos_min(phi_s, parameters, scalar) + np.pi
    return angle 

def phi_c_neg_max(phi_s, parameters, scalar):
    angle = phi_c_pos_max(phi_s, parameters, scalar) + np.pi
    return angle


#Carina Geometry

#carina_parameters = [R_p, R_d, R_o, R_i, iota_b]

def carina_unit_vectors(phi_s, phi_c, carina_parameters):
    iota_b = carina_parameters[4]
    
    phi_tilde = np.max([0, np.min([phi_s, iota_b])])
    
    e_circ = np.array([-1*np.sin(phi_tilde) * np.cos(phi_c),
                       -1*np.sin(phi_c),
                       np.cos(phi_tilde) * np.cos(phi_c)])
    
    e_stream = np.array([np.cos(phi_tilde),
                         0,
                         np.sin(phi_tilde)])
    
    normal = np.array([-1*np.sin(phi_tilde) * np.sin(phi_c),
                       np.cos(phi_c),
                       np.cos(phi_tilde) * np.sin(phi_c)])
    
    return e_circ, e_stream, normal

def unit_vector_transform(vector, phi_s, phi_c, carina_parameters):
    upsilon_angle = upsilon(phi_s, carina_parameters)
    e_circ, e_stream, normal = carina_unit_vectors(phi_s, phi_c, carina_parameters)
       
    A = np.array([[0,            -1*e_circ[2], e_circ[1]],
                  [e_circ[2],    0,             -1*e_circ[0]],
                  [-1*e_circ[1], e_circ[0],     0]])
    
    
    Q = np.identity(3) + np.sin(upsilon_angle)*A + (1-np.cos(upsilon_angle)) * A@A
    
    return Q @ vector

def rho_omega_theta_parameters(phi_s, parameters, carina_parameters, scalar):
    phi_min = phi_c_pos_min(phi_s, parameters, scalar)
    def phi_c_pos_min_lambda(phi_s, parameters, scalar):
        return -np.infty
    def phi_c_pos_max_lambda(phi_s, parameters, scalar):
        return np.infty
    
    r_lambda, _ = r_positive(phi_s, phi_min, phi_c_pos_min_lambda, phi_c_pos_max_lambda, parameters, scalar)
    _, _, n_prime = carina_unit_vectors(phi_s, phi_min, carina_parameters)
    n_lambda = unit_vector_transform(n_prime, phi_s, phi_min, carina_parameters)
    
    
    if type(r_lambda) == type(None):
        print(phi_s)
        print(phi_min)
        print(parameters)
        print(scalar)
        
    assert (type(r_lambda)) != type(None)
    
    r_kappa = r_lambda - (r_lambda[2] / n_lambda[2])*n_lambda
    r_M     = r_lambda - (r_lambda[1] / n_lambda[1])*n_lambda
    r_N = np.array([r_M[0], 0, 0])
    
    rho = np.linalg.norm((r_kappa - r_lambda))
    omega = np.arctan(r_M[2]/(np.linalg.norm((r_kappa - r_N))) )
    
    sin_theta = (r_kappa[0] - r_N[0])/(np.linalg.norm((r_kappa - r_N)))
    cos_theta = r_kappa[1]/(np.linalg.norm((r_kappa - r_N)))
    
    return rho, omega, sin_theta, cos_theta

def upsilon(phi_s, carina_parameters):
    R_p = carina_parameters[0]; R_d = carina_parameters[1]
    R_o = carina_parameters[2]; iota_b = carina_parameters[4]
    
    a = (iota_b - phi_s)/iota_b
    b = (R_p - R_d)/(R_o * iota_b)
    
    upsilon = np.arctan(-6 * (phi_s/iota_b)*a*b)
    
    return upsilon

def carina_positive(phi_s, psi, parameters, carina_parameters, scalar):
    phi_min = phi_c_pos_min(phi_s, parameters, scalar)
    def phi_c_pos_min_lambda(phi_s, parameters, scalar):
        return -1*np.pi/2
    def phi_c_pos_max_lambda(phi_s, parameters, scalar):
        return 5*np.pi/2
    
    r_lambda,_ = r_positive(phi_s, phi_min, phi_c_pos_min_lambda, phi_c_pos_max_lambda, parameters, scalar)
    _, _, n_prime = carina_unit_vectors(phi_s, phi_min, carina_parameters)
    n_lambda = unit_vector_transform(n_prime, phi_s, phi_min, carina_parameters)
    r_kappa = r_lambda - (r_lambda[2] / n_lambda[2])*n_lambda
    rho, omega, sin_theta, cos_theta = rho_omega_theta_parameters(phi_s, parameters, carina_parameters, scalar)
    
    x_carina_pos = r_kappa[0] - (rho * np.cos(psi)*sin_theta)
    y_carina_pos = r_kappa[1] - (rho * np.cos(psi)*cos_theta)
    z_carina_pos = rho*np.sin(psi)
    
    xyz_car_pos = np.array([x_carina_pos, y_carina_pos, z_carina_pos])
    
    
    normal_car_pos = -1*(r_kappa - xyz_car_pos) / np.linalg.norm((r_kappa - xyz_car_pos))*np.sign(n_lambda[2])
    
    
    return xyz_car_pos, normal_car_pos

def carina_negative(phi_s, psi, parameters, carina_parameters, scalar):
    carina_pos, pos_normal = carina_positive(phi_s, psi, parameters, carina_parameters, scalar)
    
    x_carina_neg = carina_pos[0]
    y_carina_neg = -1*carina_pos[1]
    z_carina_neg = carina_pos[2]
    
    xyz_car_neg = np.array([x_carina_neg, y_carina_neg, z_carina_neg])
    
    neg_normal = np.zeros(np.shape(pos_normal))
    neg_normal[0] = pos_normal[0]
    neg_normal[1] = -1*pos_normal[1]
    neg_normal[2] = pos_normal[2]
    
    return xyz_car_neg, neg_normal
