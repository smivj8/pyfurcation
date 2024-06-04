"""
Samuel Smith
MATH6019 - Nonlinear Programming
Date: Thu Feb  1 18:28:00 2024

Assignment
Problem

Description:

Author: smit3sl@mail.uc.edu
"""

###################################################################
#           DO NOT USE THIS ANY MORE. OUTDATED                    #
###################################################################


import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial as spatial
#from stl import mesh
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
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

# def bending_pipe_negative_z(phi_s, phi_c, L_p, R_o, R_i, R_p, R_d, iota_b):
#     x_m_neg = L_p + R_o * np.sin(phi_s) + np.sin(phi_s) * R_i(phi_s, R_p, R_d, iota_b) * np.sin(phi_c)
#     y_m_neg = R_i(phi_s, R_p, R_d, iota_b) * np.cos(phi_c)
#     z_m_neg = -1 * R_o * (1-np.cos(phi_s)) + np.cos(phi_s)*R_i(phi_s, R_p, R_d, iota_b)*np.sin(phi_c)
#     return np.array([x_m_neg, y_m_neg, z_m_neg])

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

# def daughter_negative_z(phi_s, phi_c, iota_b, L_p, R_i, R_o, R_p, R_d):
#     xyz_m_neg = bending_pipe_negative_z(iota_b, phi_c, L_p, R_o, R_i, R_p, R_d, iota_b)
    
#     zeta = R_o * np.tan(phi_s - iota_b)
#     x_d_neg = xyz_m_neg[0] + zeta*np.cos(iota_b)
#     y_d_neg = xyz_m_neg[1]
#     z_d_neg = xyz_m_neg[2] - zeta*np.sin(iota_b)
    
#     return np.array([x_d_neg, y_d_neg, z_d_neg])

# Parameters = [R_p, L_p, R_d, L_d, R_o, R_i, iota_b, delta_alpha, iota_gamma]

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
    
# def r_negative(phi_s, phi_c, phi_c_neg_min, phi_c_neg_max, parameters):
#     R_p = parameters[0]; L_p = parameters[1]; R_d = parameters[2]
#     L_d = parameters[3]; R_o = parameters[4]; R_i = parameters[5]
#     iota_b = parameters[6]; delta_alpha = parameters[7]
#     iota_gamma = parameters[8]
    
#     phi_c_min = phi_c_neg_min(phi_s, parameters)
#     phi_c_max = phi_c_neg_max(phi_s, parameters)
    
    
#     if phi_s < 0 and (phi_c >= (-1*np.pi/2) and phi_c <= 3*np.pi/2):
#         xyz = parent_pipe_section(phi_s, phi_c, R_p, L_p, R_o)
#         return xyz
#     elif (phi_s >= 0 and phi_s < iota_b) and (phi_c >= phi_c_min and phi_c <= phi_c_max):
#         xyz = bending_pipe_negative_z(phi_s, phi_c, L_p, R_o, R_i, R_p, R_d, iota_b)
#         return xyz
#     elif phi_s >= iota_b and (phi_c >= phi_c_min and phi_c <= phi_c_max):
#         xyz = daughter_negative_z(phi_s, phi_c, iota_b, L_p, R_i, R_o, R_p, R_d)
#         return xyz
#     else:
#         return
    
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

# TODO Create cap for z-neg daughter pipe 
    
def generate_pointclouds(parameters, n_circ = 100, n_streamline = 120, 
                           n_car_range = 80, n_car_rad = 8, scalar = 1.1):
    R_p = parameters[0]; L_p = parameters[1]; R_d = parameters[2]
    L_d = parameters[3]; R_o = parameters[4]; R_i = parameters[5]
    iota_b = parameters[6]; iota_gamma = parameters[8]
    
    L_p = 0.05*L_d
    L_d = 1.05*L_d
    
    carina_parameters = [R_p, R_d, R_o, R_i, iota_b]

    phi_s_pos_min = np.arctan(-1*L_p/R_o)  
    #phi_s_neg_min = 0

    eta = np.tan(iota_gamma - iota_b)
    phi_s_pos_max = iota_b + np.arctan(eta + (L_d/R_o))
    #phi_s_neg_max = phi_s_pos_max
    
    phi_s_pos = np.linspace(phi_s_pos_min, phi_s_pos_max, n_streamline)
    # phi_s_neg = np.linspace(phi_s_neg_min, phi_s_neg_max, n_streamline)
    phi_c = np.linspace(-1*np.pi/2, 3*np.pi/2, n_circ)
    
    phi_s_carina_max = iota_gamma
    phi_s_carina = np.linspace(0, phi_s_carina_max, n_car_range)
    
    
    xyz_pos = np.zeros([n_circ*n_streamline, 3])
    xyz_pos_normals = np.zeros([n_circ*n_streamline, 3])
    
    car_pos = np.zeros([n_car_range*n_car_rad, 3])
    car_pos_normals =  np.zeros([n_car_range*n_car_rad, 3])
    
    car_neg = np.zeros([n_car_range*n_car_rad, 3])
    car_neg_normals = np.zeros([n_car_range*n_car_rad, 3])


    # for j in range(len(phi_c)):
    #     for i in range(len(phi_s_pos)):
    #         xyz_pos[j*n_streamline+i,:], xyz_pos_normals[j*n_streamline+i,:] = r_positive(phi_s_pos[i], phi_c[j], 
    #                                          phi_c_pos_min, phi_c_pos_max, 
    #                                          parameters, scalar)
            
    # Radial Test
    iteration = 0
    for i in range(len(phi_s_pos)):
        for j in range(len(phi_c)):
            xyz_pos[iteration,:], xyz_pos_normals[iteration,:] = r_positive(phi_s_pos[i], phi_c[j], 
                                             phi_c_pos_min, phi_c_pos_max, 
                                             parameters, scalar)
            iteration += 1
    
    for m in range(len(phi_s_carina)):
        _, omega, _, _ = rho_omega_theta_parameters(phi_s_carina[m], parameters, carina_parameters, scalar)
        psi_range = np.linspace(-1*omega, omega, n_car_rad)
        for n in range(len(psi_range)):
            car_pos[n*n_car_range + m, :], car_pos_normals[n*n_car_range + m, :] = carina_positive(
                phi_s_carina[m], psi_range[n], parameters, carina_parameters, scalar)
            car_neg[n*n_car_range + m, :], car_neg_normals[n*n_car_range + m, :] = carina_negative(
                phi_s_carina[m], psi_range[n], parameters, carina_parameters, scalar)
    
    
    xyz_pc_positive = np.concatenate((xyz_pos, xyz_pos_normals), axis = 1)
    
    carina_pc_positive = np.concatenate((car_pos, car_pos_normals), axis = 1)
    carina_pc_negative = np.concatenate((car_neg, car_neg_normals), axis = 1)
    
    xyz_pos_trimmed = xyz_pc_positive[~np.isnan(xyz_pc_positive).any(axis=1),:] 
    xyz_neg_trimmed = np.zeros(np.shape(xyz_pos_trimmed))

    xyz_neg_trimmed[:,0] = xyz_pos_trimmed[:,0]
    xyz_neg_trimmed[:,1] = xyz_pos_trimmed[:,1]
    xyz_neg_trimmed[:,2] = -1*xyz_pos_trimmed[:,2]
    xyz_neg_trimmed[:,3] = xyz_pos_trimmed[:,3]
    xyz_neg_trimmed[:,4] = xyz_pos_trimmed[:,4]
    xyz_neg_trimmed[:,5] = -1*xyz_pos_trimmed[:,5]
    
    car_pos_trimmed = carina_pc_positive[~np.isnan(carina_pc_positive).any(axis=1),:]
    car_neg_trimmed = carina_pc_negative[~np.isnan(carina_pc_negative).any(axis=1),:]
    
    # bifurcation_pointcloud = np.concatenate((xyz_pos_trimmed, xyz_neg_trimmed))
    # carina_pointcloud = np.concatenate((car_pos_trimmed, car_neg_trimmed))
    
    
    pointcloud = np.concatenate((xyz_pos_trimmed, xyz_neg_trimmed, car_pos_trimmed, car_neg_trimmed))
    
    #Shift pointcloud to x: 0
    x_min = np.min(pointcloud[:,0])
    
    pointcloud[:,0] = pointcloud[:,0] - x_min - L_p
    

    
    return pointcloud



def pointcloud_to_open3d(pointcloud):
    pcd = o3d.geometry.PointCloud()
    if len(pointcloud[0]) == 6:
        pcd.points = o3d.utility.Vector3dVector(pointcloud[:,:3])
        pcd.normals = o3d.utility.Vector3dVector(pointcloud[:,3:])
    else:
        pcd.points = o3d.utility.Vector3dVector(pointcloud)
    return pcd

def open3d_pointcloud_visualization(pointcloud):
    pcd = pointcloud_to_open3d(pointcloud)

    o3d.visualization.draw_geometries([pcd], window_name = "PointCloud")
    return

def ball_pivoting_algorithm(pointcloud, radius_scale = 3, optimization = False, visualization = False):
    #REQUIRES NORMALS WHICH I DO NOT HAVE??? HOW DO I GENERATE THEM
    
    pcd = pointcloud_to_open3d(pointcloud)
    
    distances = pcd.compute_nearest_neighbor_distance()
    avg_distance = np.mean(distances)
    radius = radius_scale * avg_distance
    
    bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd, o3d.utility.DoubleVector([radius, radius*2]))
    
    # if optimization:
        
    
    if visualization:
        o3d.visualization.draw_geometries([bpa_mesh], window_name = "BPA_Meshing")
    
    return bpa_mesh, pcd

def poisson_algorithm(pointcloud, depth = 7, width = 0, scale = 1.1, 
                      linear_fit = False, visualization = False):
    #Also requires normals
    
    pcd = pointcloud_to_open3d(pointcloud)
    
    #computing the mesh
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth = depth, width = width,
                                                                             scale = scale, linear_fit = linear_fit)[0]
    
    
    #visualization???
    #poisson_mesh.compute_triangle_normals()
    
    #decimate to 100,000 triangles
    #poisson_mesh.simplify_quadric_decimation(100000)
    
    #visualize
    if visualization:
        o3d.visualization.draw_geometries([poisson_mesh], window_name = "poisson_Meshing")
    
    return poisson_mesh, pcd



def rotation_matrix_about_y(theta):
    #If anyone ever reads this source code, I accidentally had a fucking 1 
    #instead of a zero in the last row of this 
    R_y = np.array([[np.cos(theta),      0, np.sin(theta)],
                    [0,                  1, 0            ],
                    [-1* np.sin(theta), 0, np.cos(theta)]])
    
    return R_y




#CLEAN CROPPING?? Praying to the stack overflow gods that I can get this to work

def sliceplane(mesh, axis, value, direction):
    # axis can be 0,1,2 (which corresponds to x,y,z)
    # value where the plane is on that axis
    # direction can be True or False (True means remove everything that is
    # greater, False means less
    # than)

    vertices = np.asarray(mesh.vertices)
    triangles = np.asarray(mesh.triangles)
    new_vertices = list(vertices)
    new_triangles = []

    # (a, b) -> c
    # c refers to index of new vertex that sits at the intersection between a,b
    # and the boundingbox edge
    # a is always inside and b is always outside
    intersection_edges = dict()

    # find axes to compute
    axes_compute = [0,1,2]
    # remove axis that the plane is on
    axes_compute.remove(axis)

    def compute_intersection(vertex_in_index, vertex_out_index):
        vertex_in = vertices[vertex_in_index]
        vertex_out = vertices[vertex_out_index]
        if (vertex_in_index, vertex_out_index) in intersection_edges:
            intersection_index = intersection_edges[(vertex_in_index, vertex_out_index)]
            intersection = new_vertices[intersection_index]
        else:
            intersection = [None, None, None]
            intersection[axis] = value
            const_1 = (value - vertex_in[axis])/(vertex_out[axis] - vertex_in[axis])
            c = axes_compute[0]
            intersection[c] = (const_1 * (vertex_out[c] - vertex_in[c])) + vertex_in[c]
            c = axes_compute[1]
            intersection[c] = (const_1 * (vertex_out[c] - vertex_in[c])) + vertex_in[c]
            assert not (None in intersection)
            # save new vertice and remember that this intersection already added an edge
            new_vertices.append(intersection)
            intersection_index = len(new_vertices) - 1
            intersection_edges[(vertex_in_index, vertex_out_index)] = intersection_index

        return intersection_index

    for t in triangles:
        v1, v2, v3 = t
        if direction:
            v1_out = vertices[v1][axis] > value
            v2_out = vertices[v2][axis] > value
            v3_out = vertices[v3][axis] > value
        else: 
            v1_out = vertices[v1][axis] < value
            v2_out = vertices[v2][axis] < value
            v3_out = vertices[v3][axis] < value

        bool_sum = sum([v1_out, v2_out, v3_out])
        # print(f"{v1_out=}, {v2_out=}, {v3_out=}, {bool_sum=}")

        if bool_sum == 0:
            # triangle completely inside --> add and continue
            new_triangles.append(t)
        elif bool_sum == 3:
            # triangle completely outside --> skip
            continue
        elif bool_sum == 2:
            # two vertices outside 
            # add triangle using both intersections
            vertex_in_index = v1 if (not v1_out) else (v2 if (not v2_out) else v3)
            vertex_out_1_index = v1 if v1_out else (v2 if v2_out else v3)
            vertex_out_2_index = v3 if v3_out else (v2 if v2_out else v1)
            # print(f"{vertex_in_index=}, {vertex_out_1_index=}, {vertex_out_2_index=}")
            # small sanity check if indices sum matches
            assert sum([vertex_in_index, vertex_out_1_index, vertex_out_2_index]) == sum([v1,v2,v3])

            # add new triangle 
            new_triangles.append([vertex_in_index, compute_intersection(vertex_in_index, vertex_out_1_index), 
                compute_intersection(vertex_in_index, vertex_out_2_index)])

        elif bool_sum == 1:
            # one vertice outside
            # add three triangles
            vertex_out_index = v1 if v1_out else (v2 if v2_out else v3)
            vertex_in_1_index = v1 if (not v1_out) else (v2 if (not v2_out) else v3)
            vertex_in_2_index = v3 if (not v3_out) else (v2 if (not v2_out) else v1)
            # print(f"{vertex_out_index=}, {vertex_in_1_index=}, {vertex_in_2_index=}")
            # small sanity check if outdices sum matches
            assert sum([vertex_out_index, vertex_in_1_index, vertex_in_2_index]) == sum([v1,v2,v3])

            new_triangles.append([vertex_in_1_index, compute_intersection(vertex_in_1_index, vertex_out_index), vertex_in_2_index])
            new_triangles.append([compute_intersection(vertex_in_1_index, vertex_out_index), 
                compute_intersection(vertex_in_2_index, vertex_out_index), vertex_in_2_index])

        else:
            assert False

    # remap indices and remove unused 

    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(np.array(new_vertices))
    mesh.triangles = o3d.utility.Vector3iVector(np.array(new_triangles))
    return mesh

def clean_crop_x(mesh, min_x = None, max_x = None, min_slice = False, max_slice = False):
    #min_x = min(min_corner[0], max_corner[0])
    #min_y = min(min_corner[1], max_corner[1])
    #max_x = max(min_corner[0], max_corner[0])
    #max_y = max(min_corner[1], max_corner[1])
    mesh_sliced = mesh
    # mesh = sliceplane(mesh, 0, min_x, False)
    if max_slice:
        mesh_sliced = sliceplane(mesh_sliced, 0, max_x, True)
    if min_slice:    
        mesh_sliced = sliceplane(mesh_sliced, 0, min_x, False)
    #mesh_sliced = sliceplane(mesh_sliced, 1, max_y, True)
    #mesh_sliced = sliceplane(mesh_sliced, 1, min_y, False)
    # mesh_sliced = mesh_sliced.paint_uniform_color([0,0,1])
    return mesh_sliced



            

# def fuse_two_meshes_archive(mesh_A, mesh_B, mesh_B_parameters, fuse_location, fuse_normal, fuse_range):
#     #  fuse distance might be wrong?
#     # is the correct fuse distance the x value of the rotated fuse vector? 
#     # wait duh yes that makes sense. leaving this comment here in case I was indeed 
#     # wrong and need to redo it
    
#     #Rotate meshes from original to crop
    
#     R1 = unit_vector_match_rotation(fuse_normal, np.array([1,0,0]))
#     fuse_distance = (R1 @ fuse_location)[0]
#     mesh_A.rotate(R1, center = np.zeros(3))
#     mesh_B.rotate(R1, center = np.zeros(3))
    
#     #Initialize Variables from meshes
#     mesh_A_vertices = np.asarray(mesh_A.vertices)
#     mesh_A_triangles = np.asarray(mesh_A.triangles)
    
#     mesh_B_vertices = np.asarray(mesh_B.vertices)
#     mesh_B_triangles = np.asarray(mesh_B.triangles)
    

    
#     new_A_triangles = []
#     new_B_triangles = []
    
#     fuse_A_vertices = []
#     fuse_B_vertices = []
    
#     # TODO add bounds to crop plane (of mesh A only, mesh B is a single bifurcation unit
#     # so it should be fine. I don't want to crop earlier generations on accident!!!
#     # Need some sort of metric to measure the yz coordinates of the triangle in relation
#     # to the fuse location. Radius of mesh B would be perfect
#     mesh_radius = mesh_B_parameters[0]
#     fuse_y = (R1 @ fuse_location)[1]
#     fuse_z = (R1 @ fuse_location)[2]

    
#     # find axes to compute
#     axis = 0    
    
#     # Evaluate mesh A triangles
#     value_A = fuse_distance - (fuse_range/2)
    
#     mesh_A = o3d.geometry.TriangleMesh()
#     mesh_A.vertices = o3d.utility.Vector3dVector(mesh_A_vertices)
#     mesh_A.triangles = o3d.utility.Vector3iVector(np.array(new_A_triangles))
    
#     mesh_B = o3d.geometry.TriangleMesh()
#     mesh_B.vertices = o3d.utility.Vector3dVector(mesh_B_vertices)
#     mesh_B.triangles = o3d.utility.Vector3iVector(np.array(new_B_triangles))
    
#     # TODO actually fuse the two meshes together !!!!!!
#     #convert list of vertices to numpy array
#     fuse_A_vertices = np.array(fuse_A_vertices)
#     fuse_B_vertices = np.array(fuse_B_vertices)
    
    
#     #Shift y,z to be about origin by subtracting from the list of vertices 
#     #the shifted fuse location. This should be right, and should be independent 
#     #of mesh radius? Think cartesian to polar using arcos (range (0, 2pi))
#     shifted_A_vertices = fuse_A_vertices - (R1 @ fuse_location)
#     shifted_B_vertices = fuse_B_vertices - (R1 @ fuse_location)

#     n_A = np.shape(fuse_A_vertices)[0]
#     n_B = np.shape(fuse_B_vertices)[0]
    
#     theta_A = np.zeros(n_A)
#     theta_B = np.zeros(n_B)
    
#     #There should be EXACTLY this many triangles I think!!
#     fuse_triangles = np.zeros((n_A + n_B, 3))
    

#     #Parameterize the yz coordinate by theta about (x, 0, 0). This gets a 1D 
#     #"ordering" of the points
#     for j in range(n_A):
#         z = shifted_A_vertices[j,2]
#         y = shifted_A_vertices[j,1]
#         H = np.sqrt(shifted_A_vertices[j,1]**2 + shifted_A_vertices[j,2]**2)
#         assert H != 0.0
#         theta_A[j] = np.arccos(z/H) + np.pi * bool(y < 0)
#     for i in range(n_B):
#         z = shifted_B_vertices[i,2]
#         y = shifted_B_vertices[i,1]
#         H = np.sqrt(shifted_B_vertices[i,1]**2 + shifted_B_vertices[i,2]**2)
#         assert H != 0.0
#         theta_B[i] = np.arccos(z/H) + np.pi * bool(y < 0)
    
#     #Sort the fuseA and fuseB vertices by thetaA and thetaB respectively
#     #also sort thetaA and theta B of course
#     sort_index_A = theta_A.argsort()
#     sort_index_B = theta_B.argsort()
    
#     sorted_A_vertices = fuse_A_vertices[sort_index_A]
#     theta_A.sort()
#     sorted_B_vertices = fuse_B_vertices[sort_index_B]
#     theta_B.sort()
#     fuse_vertices = np.concatenate((sorted_A_vertices, sorted_B_vertices))
#     index_A = 0
#     index_B = 0
#     downward_indices = np.zeros(n_A)
    
#     # print(theta_A)
#     # print(theta_B)
    
#     #Create DOWNWARD TRIANGLES (there will be nA total, create the first nA-1)
#     for d in range(n_A - 1):
#         v1_index = d
#         v2_index = d+1
#         #print(d); print(index_B); print("\n")
#         index_B = 0
#         while theta_B[index_B] <= theta_A[d]:
#             index_B += 1
#         v3_index = n_A + index_B
#         downward_indices[d] = index_B
#         fuse_triangles[d] = np.array([v1_index, v2_index, v3_index])
    
#     #create the final downward triangle
#     v1_index = 0; v2_index = (n_A - 1); v3_index = (n_A + n_B - 1)
#     downward_indices[n_A - 1] = n_B - 1
#     fuse_triangles[n_A - 1] = np.array([v1_index, v2_index, v3_index])
    
    
    
#     #Okay, I have them sorted. Now, I need to finish up and create the triangles
#     #I am tired and don't feel like thinking through this logic right now tbh
#     #so problem for tomorrow!
    
#     mesh_fuse = o3d.geometry.TriangleMesh()
#     mesh_fuse.vertices = o3d.utility.Vector3dVector(fuse_vertices)
#     mesh_fuse.triangles = o3d.utility.Vector3iVector(fuse_triangles)
    
#     mesh = mesh_A + mesh_B + mesh_fuse
    
#     R2 = unit_vector_match_rotation(np.array([1,0,0]), fuse_normal)
#     mesh.rotate(R2, center = np.zeros(3))
   

#     return mesh


def get_free_edges(mesh):
    # Highlight Boundaries!
    
    #mesh.compute_vertex_normals()
    geoms = [mesh]
    free_edges = mesh.get_non_manifold_edges(allow_boundary_edges=False)
    o3d.visualization.draw_geometries(geoms, window_name = "Visualization",
                                      mesh_show_wireframe = True, mesh_show_back_face = True)
    return



def outlet_positions(parameters):
    R_p = parameters[0]; L_p = parameters[1]; R_d = parameters[2]
    L_d = parameters[3]; R_o = parameters[4]; R_i = parameters[5]
    iota_b = parameters[6]; delta_alpha = parameters[7]
    iota_gamma = parameters[8]
    
    #Crop Daughter
        # z pos crop
        #find x_max
    eta = np.tan(iota_gamma - iota_b)
    phi_s_max = iota_b + np.arctan(eta + (L_d/R_o))
    def phi_c_pos_min_lambda(phi_s, parameters, scalar):
        return -1*np.pi/2
    def phi_c_pos_max_lambda(phi_s, parameters, scalar):
        return 3*np.pi/2
    r0, _ =  r_positive(phi_s_max, 0, phi_c_pos_min_lambda, phi_c_pos_max_lambda, parameters, 1.1)
    rpi, _ = r_positive(phi_s_max, np.pi, phi_c_pos_min_lambda, phi_c_pos_max_lambda, parameters, 1.1)
    r_D_pos = ((r0 + rpi)/2)
    
    r_D_neg = np.zeros(3)
    
    r_D_neg[0] = r_D_pos[0]
    r_D_neg[1] = r_D_pos[1]
    r_D_neg[2] = -1*r_D_pos[2]
    
    return r_D_pos, r_D_neg

def outlet_normals(parameters):
    """
    Function that returns the vectors normal to the z+ and z- outlet planes.
    Used to define the rotation matrices for transforming for cropping and attaching 
    subsequent bifurcation units
    Parameters
    ----------
    parameters : List
        Bifurcation Unit Parameters. In form:
        parameters = [InletRadius, 0, OutletRadii, DaughterBranchingLength,
                      OuterRadius, R_i, iota_b, delta_alpha, iota_gamma]
        See paper for explanation of variables

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    L_d = parameters[3]; R_o = parameters[4]; 
    iota_b = parameters[6]
    iota_gamma = parameters[8]
    
    eta = np.tan(iota_gamma - iota_b)
    phi_s_max = iota_b + np.arctan(eta + (L_d/R_o))
    def phi_c_pos_min_lambda(phi_s, parameters, scalar):
        return -1*np.pi/2
    def phi_c_pos_max_lambda(phi_s, parameters, scalar):
        return 3*np.pi/2
    r_C, _ = r_positive(iota_b, 0, phi_c_pos_min_lambda, phi_c_pos_max_lambda, parameters, 1.1)
    r_F, _ = r_positive(phi_s_max, 0, phi_c_pos_min_lambda, phi_c_pos_max_lambda, parameters, 1.1)
    n_D_pos = ((r_F - r_C)/np.linalg.norm((r_F - r_C)))
    
    n_D_neg = np.zeros(3)
    
    n_D_neg[0] = n_D_pos[0]
    n_D_neg[1] = n_D_pos[1]
    n_D_neg[2] = -1*n_D_pos[2]
    
    return n_D_pos, n_D_neg

def unit_vector_match_rotation(a, b):
    #ROTATES UNIT VECTOR A ONTO UNIT VECTOR B
    u = np.cross(a, b)
    s = np.linalg.norm(u)
    c = np.dot(a, b)
    
    if np.linalg.norm(b - a) == 0:
        R = np.identity(3)
    else:
        U_mat = np.array([[0,      -1*u[2], u[1]   ],
                          [u[2],    0,      -1*u[0]],
                          [-1*u[1], u[0],   0      ]])
    
    #V_mat = np.cross(v, v)
    
        R = np.identity(3) + U_mat + (U_mat @ U_mat) * ((1-c)/s**2)
    
    return R

def rotation_about_vector_u(u, theta):
    U_mat = np.array([[0,      -1*u[2], u[1]   ],
                      [u[2],    0,      -1*u[0]],
                      [-1*u[1], u[0],   0      ]])
    
    R = np.cos(theta)*np.identity(3) + np.sin(theta)*U_mat + (1-np.cos(theta)) * np.outer(u, u)
    
    return R

def crop_bifurcation_poisson_mesh(mesh, parameters):
    iota_b = parameters[6]
    
    #crop inlet end
    mesh = clean_crop_x(mesh, min_x = 0, min_slice = True)
    
    
    #Crop Daughter ends
        # z pos crop
        #find x_max
    
    r_D_pos, _ = outlet_positions(parameters)
    n_D_pos, n_D_neg = outlet_normals(parameters)
    cos_theta = np.dot(n_D_pos, r_D_pos)/np.linalg.norm(r_D_pos)
    hypot_max = np.linalg.norm(r_D_pos)*cos_theta
        #rotate poisson cloud and do the first daughter mesh crop

    mesh.rotate((rotation_matrix_about_y(iota_b)), center = np.zeros(3))
    #mesh.rotate((unit_vector_match_rotation(np.array([1,0,0]), n_D_pos)), 
    #            center = np.zeros(3))
    mesh = clean_crop_x(mesh, max_x = hypot_max, max_slice = True)
    
    mesh.rotate((rotation_matrix_about_y(-1*iota_b)), center = np.zeros(3))
    mesh.rotate((rotation_matrix_about_y(-1*iota_b)), center = np.zeros(3))
    #mesh.rotate((unit_vector_match_rotation(n_D_pos, n_D_neg)), 
    #            center = np.zeros(3))

    
    mesh = clean_crop_x(mesh, max_x = hypot_max, max_slice = True)
    
    mesh.rotate((rotation_matrix_about_y(iota_b)), center = np.zeros(3))
    #mesh.rotate((unit_vector_match_rotation(n_D_neg, np.array([1,0,0]))), 
    #            center = np.zeros(3))


    return mesh



def attach_bifurcation_unit_to_parent_mesh(parent_mesh, bifurcation_mesh,
                                           axial_rotation_angle, r_parent, 
                                           n_parent, bifurcation_parameters):
    
    r_pos, r_neg = outlet_positions(bifurcation_parameters)
    n_pos, n_neg = outlet_normals(bifurcation_parameters)
    
    #Define rotation matrices:
        #R_1 rotates the bifurcation mesh such that it matches the outlet of the 
        #parent mesh
        #R_2 axially rotates the bifurcation about the inlet normal axis
        
    R_1 = unit_vector_match_rotation(np.array([1,0,0]), n_parent)
    R_2 = rotation_about_vector_u(n_parent, axial_rotation_angle)

    #Rotate the bifurcation mesh to match up with the parent mesh, then axial 
    #rotation about its inlet normal
    bifurcation_mesh.rotate(R_1, center = np.array([0,0,0]))
    bifurcation_mesh.rotate(R_2, center = np.array([0,0,0]))
    #Translate Bifurcation Mesh
    bifurcation_mesh.translate(r_parent)
    
    #update new outlet positions and normals
    r_outlet_pos = (R_2 @ (R_1 @ r_pos)) + r_parent
    n_outlet_pos = R_2 @ (R_1 @ n_pos)
    
    r_outlet_neg = (R_2 @ (R_1 @ r_neg)) + r_parent
    n_outlet_neg = R_2 @ (R_1 @ n_neg)


    
    mesh = parent_mesh + bifurcation_mesh
    
    return mesh, r_outlet_pos, r_outlet_neg, n_outlet_pos, n_outlet_neg


def truncate_outlet(mesh, mesh_parameters, positive = True, negative = False):
    r_pos, _ = outlet_positions(mesh_parameters)
    n_pos, _ = outlet_normals(mesh_parameters)
    iota_b = mesh_parameters[6]
    
    
    #Crop Daughter ends
        # z pos crop
        #find x_max
    

    cos_theta = np.dot(n_pos, r_pos)/np.linalg.norm(r_pos)
    hypot_max = np.linalg.norm(r_pos)*cos_theta
    truncation = hypot_max*0.975
        #rotate poisson cloud and do the first daughter mesh crop

    if positive:
        mesh.rotate((rotation_matrix_about_y(iota_b)), center = np.zeros(3))
        mesh = clean_crop_x(mesh, max_x = truncation, max_slice = True)
        mesh.rotate((rotation_matrix_about_y(-1*iota_b)), center = np.zeros(3))
        
    if negative:
        mesh.rotate((rotation_matrix_about_y(-1*iota_b)), center = np.zeros(3))
        mesh = clean_crop_x(mesh, max_x = truncation, max_slice = True)
        mesh.rotate((rotation_matrix_about_y(iota_b)), center = np.zeros(3))

    return mesh

#should I have a class before this, one that 

class bifurcation_unit:
    
    # def __init__(self, bifurcation_unit_parameters):
    #     self.bifurcation_unit_parameters = bifurcation_unit_parameters
    
    def __init__(self, bifurcation_unit_parameters):
        #Create and Crop Bifurcation Mesh
        n_cont_outlets = bifurcation_unit_parameters[-1]
        
        xyz = generate_pointclouds(bifurcation_unit_parameters)
        mesh, _ = poisson_algorithm(xyz)
        mesh = crop_bifurcation_poisson_mesh(mesh, bifurcation_unit_parameters)
        
        #ASK: do all these functions now go within this class?? Or only future
        #ones (like .rotate, for example)
        
        if n_cont_outlets == 0:
            print("TODO: Cap Outlets")
            #TODO ^^ (need the get free edges first tho)
            
        elif n_cont_outlets == 1:
            mesh = truncate_outlet(mesh, bifurcation_unit_parameters, 
                                   positive = True, negative = False)
            print("TODO: cap negative outlet")
        elif n_cont_outlets == 2:
            mesh = truncate_outlet(mesh, bifurcation_unit_parameters, 
                                   positive = True, negative = True)
        else:
            print("Invalid Input: n_cont_outlets must be == 0, 1, or 2")
            assert False
        
        
        #Get position and normal vectors
        position_vectors = outlet_positions(bifurcation_unit_parameters)
        normal_vectors = outlet_normals(bifurcation_unit_parameters)
        
        #Get Vertices for inlet and outlet free edges
        #TODO ^^
        inlet_free_edge_vertices = None
        outlet_free_edge_vertices = None
        
        self.bifurcation_unit_parameters = [mesh, position_vectors, normal_vectors,
                                            inlet_free_edge_vertices,
                                            outlet_free_edge_vertices]
        """
        How do I set this up so that I can operate on these without pulling apart the list?
        i.e. if I have a function that acts on the mesh and position vectors separately, 
        do I:
        
        def function(bifurcation_unit):
            mesh = bifurcation_mesh[0]
            #do something to mesh
            position_vectors = bifurcation_unit[1]
            #do something
            bifurcation unit[0] = mesh
            bifurcation_unit[1] = position_vectors
            return bifurcation_unit
            
        OR Can I do
        def function(bifurcation_unit):
            bifurcation_unit.mesh 
            bifurcation_unit.position_vectors
            return
        
        In which case do I have to define each mesh, etc. as their own class
        first???
            
        
        """
        
        return
        
    """
    Pseudo - Code
    Okay, here is where I want to have the functions that act on the bifurcation 
    unit as it is defined above, with the mesh, vectors, and free edges. AN easy
    example, is a rotation/translation
    
    def rotate(self, rotation_matrix):
        unit.mesh = mesh.rotate()
        unit.position = position.rotate()
        ... etc for each component of the mesh
        return unit
    
    So that I can do
    bifurcation_unit = bifurcation_unit.rotate(matrix)
    
    
    """
    

class airway_tree:
    def __init(self):
        return
        
    
    


#%% Test pointcloud function
D_0 = 100; R_0 = D_0/2; R_1 = D_0/2.48888; L_0 = 0; L_1 = 4.08888*R_1
R_o_0 = 5*R_1; iota_b_test = 35*np.pi/180; iota_gamma_test = 39.83 * np.pi/180
delta_alpha_test = 3.6*np.pi/180

D_1 = R_1*2; R_2 = D_1/2.48888; L_0 = 0; L_2 = 4.08888*R_2
R_o_1 = 5*R_2

# D_2 = R_2*2; R_3 = D_2/2.48888; L_0 = 0; L_3 = 4.08888*R_3
# R_o_2 = 5*R_3;

parameters_0 = [R_0, 0, R_1, L_1, R_o_0, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]

parameters_1 = [R_1, 0, R_2, L_2, R_o_1, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]

# parameters_2 = [R_2, 0, R_3, L_3, R_o_2, R_i, iota_b_test, delta_alpha_test, iota_gamma_test]


# n_gens_upper = 23
# upper_tree_parameters = [0]*n_gens_upper
# upper_tree_parameters[0] = parameters_0
# for j in range(1, n_gens_upper):
#     parameters = [0]*9
#     parameters[0] = upper_tree_parameters[j-1][2]
#     parameters[2] = parameters[0]*2/2.48888
#     parameters[3] = parameters[2]*4.08888
#     parameters[4] = parameters[2]*5
#     parameters[5] = R_i
#     parameters[6] = iota_b_test
#     parameters[7] = delta_alpha_test
#     parameters[8] = iota_gamma_test
#     upper_tree_parameters[j] = parameters




r0_pos, r0_neg = outlet_positions(parameters_0)
n0_pos, n0_neg = outlet_normals(parameters_0)

xyz_0 = generate_pointclouds(parameters_0)
mesh_0, pcd_0 = poisson_algorithm(xyz_0)
mesh_0 = crop_bifurcation_poisson_mesh(mesh_0, parameters_0)

xyz_1 = generate_pointclouds(parameters_1)
mesh_1, pcd_1 = poisson_algorithm(xyz_1)
mesh_1 = crop_bifurcation_poisson_mesh(mesh_1, parameters_1)

mesh_0 = truncate_outlet(mesh_0, parameters_0, positive = True)

combined_mesh, r1_pos, r1_neg, n1_pos, _ = attach_bifurcation_unit_to_parent_mesh(mesh_0, mesh_1, np.pi/2,
                                                                              r0_pos, n0_pos, 
                                                                              parameters_1)
origin = np.zeros((4,3))
origin[1] = r0_pos; origin[2] = r1_pos; origin[3] = r1_neg
origin_pcd = pointcloud_to_open3d(origin)

o3d.visualization.draw_geometries([combined_mesh, origin_pcd], window_name = "Visualization",
                                    mesh_show_wireframe = True, mesh_show_back_face = True)

# #%% Mesh Combination
# mesh = combined_mesh
# free_edges = mesh.get_non_manifold_edges(allow_boundary_edges=False)
# # o3d.visualization.draw_geometries([mesh], window_name = "Visualization",
# #                                   mesh_show_wireframe = True, mesh_show_back_face = True)

# edges = np.asarray(free_edges)
# mesh_vertices = np.asarray(mesh.vertices)

# edge_indices = np.unique(np.concatenate((edges[:,0], edges[:,1])))

# edge_vertices = mesh_vertices[edge_indices]

# edge_vertices_pcd = pointcloud_to_open3d(edge_vertices)

# target_edge_indices = []

# # o3d.visualization.draw_geometries([mesh, edge_vertices_pcd], window_name = "Visualization",
# #                                   mesh_show_wireframe = True, mesh_show_back_face = True)




# R1 = unit_vector_match_rotation(n0_pos, np.array([1,0,0]))
# fuse_distance = (R1 @ r0_pos)[0]
# mesh.rotate(R1, center = np.zeros(3))


# mesh_radius = parameters_1[0]
# fuse_y = (R1 @ r0_pos)[1]
# fuse_z = (R1 @ r0_pos)[2]


# # Evaluate mesh A triangles
# bifurcation_length = parameters_1[3]
# fuse_range = bifurcation_length * 0.05
# value_A = fuse_distance - (fuse_range/2)
# value_B = fuse_distance + (fuse_range/2)

# for j in range(len(edge_indices)):
#     if mesh_vertices[edge_indices[j]][0] >= value_A and mesh_vertices[edge_indices[j]][0] <= value_B:
#         target_edge_indices.append(edge_indices[j])

# target_edge_vertices = mesh_vertices[target_edge_indices]
# target_edge_vertices_pcd = pointcloud_to_open3d(np.array(target_edge_vertices))

# distances = target_edge_vertices_pcd.compute_nearest_neighbor_distance()
# avg_distance = np.mean(distances)
# eps = avg_distance*0.75

# o3d.visualization.draw_geometries([mesh, target_edge_vertices_pcd], window_name = "Visualization",
#                                   mesh_show_wireframe = True, mesh_show_back_face = True)

