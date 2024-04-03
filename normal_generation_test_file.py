"""
Created on Tue Feb 27 19:56:55 2024

@author: SMIVJ8

Samuel Smith
MATH6019 - Nonlinear Programming
Date: Thu Feb  1 18:28:00 2024

Assignment
Problem

Description:

Author: smit3sl@mail.uc.edu
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial as spatial
#from stl import mesh
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
#import surf2stl
import open3d as o3d

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
    L_d = parameters[3]; R_o = parameters[4]; R_i = parameters[5]
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

#carina_parameters = [R_p, R_d, R_o, R_i, iota_b, scalar]

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
        return -1*np.pi/2
    def phi_c_pos_max_lambda(phi_s, parameters, scalar):
        return 5*np.pi/2
    
    r_lambda, _ = r_positive(phi_s, phi_min, phi_c_pos_min_lambda, phi_c_pos_max_lambda, parameters, scalar)
    _, _, n_prime = carina_unit_vectors(phi_s, phi_min, carina_parameters)
    n_lambda = unit_vector_transform(n_prime, phi_s, phi_min, carina_parameters)
    
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
    R_o = carina_parameters[2]; R_i = carina_parameters[3]
    iota_b = carina_parameters[4]
    
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
    
    return xyz_car_pos, n_lambda*(1/np.linalg.norm(n_lambda))

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


    
    
def generate_pointclouds(parameters, n_circ = 100, n_streamline = 150, 
                           n_car_range = 48, n_car_rad = 8, scalar = 1.1):
    R_p = parameters[0]; L_p = parameters[1]; R_d = parameters[2]
    L_d = parameters[3]; R_o = parameters[4]; R_i = parameters[5]
    iota_b = parameters[6]; delta_alpha = parameters[7]
    iota_gamma = parameters[8]
    
    carina_parameters = [R_p, R_d, R_o, R_i, iota_b]

    phi_s_pos_min = 0 #np.arctan(-1*L_p/R_o)  
    phi_s_neg_min = 0

    eta = np.tan(iota_gamma - iota_b)
    phi_s_pos_max = iota_b + np.arctan(eta + (L_d/R_o))
    phi_s_neg_max = phi_s_pos_max
    
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

    # for i in range(len(phi_s_pos)):
    #     for j in range(len(phi_c)):
    #         xyz_pos[j*n_streamline+i,:] = r_positive(phi_s_pos[i], phi_c[j], 
    #                                        phi_c_pos_min, phi_c_pos_max, 
    #                                        parameters)
            
    for j in range(len(phi_c)):
        for i in range(len(phi_s_pos)):
            xyz_pos[j*n_streamline+i,:], xyz_pos_normals[j*n_streamline+i,:] = r_positive(phi_s_pos[i], phi_c[j], 
                                             phi_c_pos_min, phi_c_pos_max, 
                                             parameters, scalar)
    
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
    
    return pointcloud

# def delaunay_surface_triangulation(parameters,  n_circ = 80, n_streamline = 80, 
#                            n_car_range = 48, n_car_rad = 16):
    
#     xyz_pos, xyz_neg, car_pos, car_neg = generate_pointclouds(parameters, 
#                                                                n_circ = n_circ, 
#                                                                n_streamline = n_streamline, 
#                                                                n_car_range = n_car_range, 
#                                                                n_car_rad = n_car_rad)
    
#     len_xyz = np.shape(xyz_pos)[0]
    
#     y_pos = np.concatenate((xyz_pos[0:int(len_xyz/2), :], xyz_neg[0:int(len_xyz/2),:], car_pos))
#     y_neg = np.concatenate((xyz_pos[int(len_xyz/2):, :], xyz_neg[int(len_xyz/2):,:], car_neg))

#     pos_sorted_index = np.lexsort((y_pos[:,0], y_pos[:,2], y_pos[:,1]))
#     neg_sorted_index = np.lexsort((y_neg[:,0], y_neg[:,2], -1*y_neg[:,1]))

#     y_pos = y_pos[pos_sorted_index]
#     y_neg = y_neg[neg_sorted_index]
    
#     #Don't forget to duplicate some of the vertices HERE, but skipping for now
#     #to try to simply get the surface triangulation
    
    
#     xz_pos_projection = np.zeros((np.shape(y_pos)[0], 2))
#     xz_neg_projection = np.zeros((np.shape(y_neg)[0], 2))
#     xz_pos_projection[:,0] = y_pos[:,0]
#     xz_pos_projection[:,1] = y_pos[:,2]
#     xz_neg_projection[:,0] = y_neg[:,0]
#     xz_neg_projection[:,1] = y_neg[:,2]
    
#     # fig = plt.figure(figsize =(10,10))
#     # ax = fig.add_subplot(111, projection='3d')
    
#     tri_pos = spatial.Delaunay(xz_pos_projection)
#     tri_neg = spatial.Delaunay(xz_neg_projection)
    
#     # ax.plot_trisurf(
#     #     y_pos[:,0], y_pos[:,1], y_pos[:,2],
#     #     triangles=tri_pos.simplices, cmap="magma"
#     # )
#     # plt.show()
    
#     # ax.plot_trisurf(
#     #     y_neg[:,0], y_neg[:,1], y_neg[:,2],
#     #     triangles=tri_neg.simplices, cmap="magma"
#     # )
#     # plt.show()
    
#     return tri_pos, tri_neg, y_pos, y_neg

def pointcloud_to_open3d(pointcloud):
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(pointcloud[:,:3])
    pcd.normals = o3d.utility.Vector3dVector(pointcloud[:,3:])
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
    
    return bpa_mesh

def poisson_algorithm(pointcloud, depth = 8, width = 0, scale = 1.1, 
                      linear_fit = False, visualization = False):
    
    pcd = pointcloud_to_open3d(pointcloud)
    
    #computing the mesh
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth = depth, width = width,
                                                                             scale = scale, linear_fit = linear_fit)
    #visualize
    if visualization:
        o3d.visualization.draw_geometries([poisson_mesh], window_name = "poisson_Meshing")
    
    return

#%% Test pointcloud function
D_p = 10; R_p = D_p/2; R_d = D_p/2.48888; L_p = 3.92857*R_p; L_d = 4.08888*R_d
R_o = 5*R_d; iota_b = 35*np.pi/180; iota_gamma = 39.83 * np.pi/180
delta_alpha = 3.6*np.pi/180

parameters_test = [R_p, L_p, R_d, L_d, R_o, R_i, iota_b, delta_alpha, iota_gamma]


#tri_pos, tri_neg, y_pos, y_neg = delaunay_surface_triangulation(parameters_test)

# pointcloud = generate_pointclouds(parameters_test)
# xyz = np.concatenate((pointcloud[0], pointcloud[1], pointcloud[2], pointcloud[3]))


#test pcloud w normals
ptc = generate_pointclouds(parameters_test)


# np.savetxt('control_points.csv', np.concatenate((pointcloud)), delimiter = ',', comments = '')


# # surf2stl.tri_write('export_test_y_pos.stl', y_pos[:,0], y_pos[:,1], y_pos[:,2], tri_pos)
# # surf2stl.tri_write('export_test_y_neg.stl', y_neg[:,0], y_neg[:,1], y_neg[:,2], tri_neg)

# color = np.arange(0, len(ptc))

# fig = plt.figure(figsize =(10,10))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(ptc[:,0],ptc[:,1],ptc[:,2], c = 'b')
# #ax.scatter(dummy[:,0],dummy[:,1],dummy[:,2], c = 'r')
# # for k in range(6):
# #     ax.text(xn[k], yn[k], zn[k], str(k))
# #ax.plot_trisurf(xn,yn,zn, triangles = tri.simplices)
# plt.title('Test')

# fig = plt.figure(figsize =(10,10))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(y_pos_unsorted[:,0],y_pos_unsorted[:,1],y_pos_unsorted[:,2])
# # for k in range(6):
# #     ax.text(xn[k], yn[k], zn[k], str(k))
# #ax.plot_trisurf(xn,yn,zn, triangles = tri.simplices)
# plt.title('Test')
# plt.show()

# fig = plt.figure(figsize =(10,10))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(y_neg_unsorted[:,0],y_neg_unsorted[:,1],y_neg_unsorted[:,2])
# # for k in range(6):
# #     ax.text(xn[k], yn[k], zn[k], str(k))
# #ax.plot_trisurf(xn,yn,zn, triangles = tri.simplices)
# plt.title('Test')


#%% Test o3d and surface triangulation algorithms

#Pointclkoud_visualization
open3d_pointcloud_visualization(ptc)

#Ball Pivoting Mesh Visualization Trial 1

#REQUIRES NORMALS
#bpa_mesh = ball_pivoting_algorithm(ptc, radius_scale = 3, visualization=False)

# #POssion Reconstruction Visualization Trial 1
# poisson_algorithm(ptc, depth = 8, width = 0, scale = 1.1, 
#                   linear_fit = False, visualization = True)


#%%TEST EXAMPLE

# xyz = np.concatenate((xyz_pos, xyz_neg))
#xyz = xyz_pos

# fig = plt.figure(figsize =(10,10))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])
# # for k in range(6):
# #     ax.text(xn[k], yn[k], zn[k], str(k))
# #ax.plot_trisurf(xn,yn,zn, triangles = tri.simplices)
# plt.title('Test')

# plt.show()

# convexhull = spatial.ConvexHull(xyz)
# faces = convexhull.simplices

# bifurcation_mesh = mesh.Mesh(
#   np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype)
# )
# for i, f in enumerate(faces):
#   for j in range(3):
#     bifurcation_mesh.vectors[i][j] = xyz[f[j],:]
    
# bifurcation_mesh.save('numpy_stl_example_03.stl')





