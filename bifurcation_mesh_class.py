from bifurcation_vectors_class import bifurcation_vectors
import numpy as np
import open3d as o3d
from utility_functions_bifurcation import *
from stylianou_geometry_equations import *


class bifurcation_mesh:

    def __init__(self, parameters):
        self.R_p = parameters[0]; self.L_p = parameters[1]; self.R_d = parameters[2]
        self.L_d = parameters[3]; self.R_o = parameters[4]; self.R_i = parameters[5]
        self.iota_b = parameters[6]; self.delta_alpha = parameters[7]
        self.iota_gamma = parameters[8]; self.n_cont_outlets = parameters[9]
        self.parameters = parameters
        self.generate_parametric_xyz()
        self.generate_o3d_pointcloud()
        self.generate_o3d_mesh() 
        self.crop_bifurcation_mesh()
        self.truncate_continuation_outlet()


    def generate_parametric_xyz(self, n_circ = 100, n_streamline = 120, 
                                n_car_range = 80, n_car_rad = 8, scalar = 1.1):
        inlet_extension = 0.05*self.L_d
        outlet_extension = 1.05*self.L_d
        carina_parameters = [self.R_p, self.R_d, self.R_o, self.R_i, self.iota_b]
        phi_s_pos_min = np.arctan(-1*inlet_extension/self.R_o)
        eta = np.tan(self.iota_gamma - self.iota_b)
        phi_s_pos_max = self.iota_b + np.arctan(eta + (outlet_extension/self.R_o))
        phi_s_pos = np.linspace(phi_s_pos_min, phi_s_pos_max, n_streamline)
        phi_c = np.linspace(-1*np.pi/2, 3*np.pi/2, n_circ)
        phi_s_carina_max = self.iota_gamma
        phi_s_carina = np.linspace(0, phi_s_carina_max, n_car_range)
        
        xyz_pos = np.zeros([n_circ*n_streamline, 3])
        xyz_pos_normals = np.zeros([n_circ*n_streamline, 3])
        
        car_pos = np.zeros([n_car_range*n_car_rad, 3])
        car_pos_normals =  np.zeros([n_car_range*n_car_rad, 3])
        
        car_neg = np.zeros([n_car_range*n_car_rad, 3])
        car_neg_normals = np.zeros([n_car_range*n_car_rad, 3])
                
        iteration = 0
        for i in range(len(phi_s_pos)):
            for j in range(len(phi_c)):
                xyz_pos[iteration,:], xyz_pos_normals[iteration,:] = r_positive(phi_s_pos[i], phi_c[j], 
                                                phi_c_pos_min, phi_c_pos_max, 
                                                self.parameters, scalar)
                iteration += 1
        
        for m in range(len(phi_s_carina)):
            _, omega, _, _ = rho_omega_theta_parameters(phi_s_carina[m], self.parameters, carina_parameters, scalar)
            psi_range = np.linspace(-1*omega, omega, n_car_rad)
            for n in range(len(psi_range)):
                car_pos[n*n_car_range + m, :], car_pos_normals[n*n_car_range + m, :] = carina_positive(
                    phi_s_carina[m], psi_range[n], self.parameters, carina_parameters, scalar)
                car_neg[n*n_car_range + m, :], car_neg_normals[n*n_car_range + m, :] = carina_negative(
                    phi_s_carina[m], psi_range[n], self.parameters, carina_parameters, scalar)

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
        
        xyz = np.concatenate((xyz_pos_trimmed, xyz_neg_trimmed, car_pos_trimmed, car_neg_trimmed))
        
        #Shift pointcloud to x: 0
        x_min = np.min(xyz[:,0])
        
        xyz[:,0] = xyz[:,0] - x_min - inlet_extension

        self.xyz = xyz
        return
    
    def generate_o3d_pointcloud(self):
        pointcloud = o3d.geometry.PointCloud()
        if len(self.xyz[0]) == 6:
            pointcloud.points = o3d.utility.Vector3dVector(self.xyz[:,:3])
            pointcloud.normals = o3d.utility.Vector3dVector(self.xyz[:,3:])
        else:
            pointcloud.points = o3d.utility.Vector3dVector(self.xyz)
        self.pointcloud = pointcloud
        return
    
    def generate_o3d_mesh(self, depth = 7, width = 0, scale = 1.1, 
                      linear_fit = False):
        poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(self.pointcloud, depth = depth, width = width,
                                                                                scale = scale, linear_fit = linear_fit)[0]
        self.mesh = poisson_mesh
        return
    
    def crop_bifurcation_mesh(self):
        cropped_mesh = clean_crop_x(self.mesh, min_x = 0, min_slice = True)
        #Crop Daughter ends
            # z pos crop
            #find x_max
        position_vec = bifurcation_vectors(self.parameters, False)
        r_pos, _ = position_vec.calculate_outlet_positions()
        normal_vec = bifurcation_vectors(self.parameters, True)
        n_pos, _ = normal_vec.calculate_outlet_normals()
        
        cos_theta = np.dot(n_pos, r_pos)/np.linalg.norm(r_pos)
        hypot_max = np.linalg.norm(r_pos)*cos_theta
        #rotate poisson cloud and do the first daughter mesh crop

        cropped_mesh.rotate((rotation_matrix_about_y(self.iota_b)), center = np.zeros(3))
        cropped_mesh = clean_crop_x(cropped_mesh, max_x = hypot_max, max_slice = True)
        
        cropped_mesh.rotate((rotation_matrix_about_y(-1*self.iota_b)), center = np.zeros(3))
        cropped_mesh.rotate((rotation_matrix_about_y(-1*self.iota_b)), center = np.zeros(3))
        cropped_mesh = clean_crop_x(cropped_mesh, max_x = hypot_max, max_slice = True)
        
        cropped_mesh.rotate((rotation_matrix_about_y(self.iota_b)), center = np.zeros(3))
        
        self.cropped_mesh = cropped_mesh

        return
    
    def truncate_continuation_outlet(self):
        position_vec = bifurcation_vectors(self.parameters, False)
        r_pos, _ = position_vec.calculate_outlet_positions()
        normal_vec = bifurcation_vectors(self.parameters, True)
        n_pos, _ = normal_vec.calculate_outlet_normals()
        cos_theta = np.dot(n_pos, r_pos)/np.linalg.norm(r_pos)
        hypot_max = np.linalg.norm(r_pos)*cos_theta
        truncation_distance = hypot_max*0.97
        truncated_mesh = self.cropped_mesh
        if self.n_cont_outlets == 1:
            truncated_mesh.rotate((rotation_matrix_about_y(self.iota_b)), center = np.zeros(3))
            truncated_mesh = clean_crop_x(truncated_mesh, max_x = truncation_distance, max_slice = True)
            truncated_mesh.rotate((rotation_matrix_about_y(-1*self.iota_b)), center = np.zeros(3))
        elif self.n_cont_outlets == 2:
            truncated_mesh.rotate((rotation_matrix_about_y(self.iota_b)), center = np.zeros(3))
            truncated_mesh = clean_crop_x(truncated_mesh, max_x = truncation_distance, max_slice = True)
            truncated_mesh.rotate((rotation_matrix_about_y(-1*self.iota_b)), center = np.zeros(3))
            truncated_mesh.rotate((rotation_matrix_about_y(-1*self.iota_b)), center = np.zeros(3))
            truncated_mesh = clean_crop_x(truncated_mesh, max_x = truncation_distance, max_slice = True)
            truncated_mesh.rotate((rotation_matrix_about_y(self.iota_b)), center = np.zeros(3))
        else:
            pass
        self.truncated_mesh = truncated_mesh
        return

    def get_free_edge_vertices(self, inlet_outlet_index):
        #inlet_outlet_index: 0 = inlet, 1 = positive_outlet, 2 = negative_outlet
        free_edges = self.truncated_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        bifurcation_free_vertex_index = np.unique(np.concatenate((numpy_free_edges[:,0], numpy_free_edges[:,1])))
        bifurcation_free_vertices = np.asarray(self.truncated_mesh.vertices)[bifurcation_free_vertex_index]
        free_vertices = []
        free_vertex_index = []
        if inlet_outlet_index == 0:
            #Get Inlet free edge vertices
            for ind, val in enumerate(bifurcation_free_vertices):
                if val[0] < np.abs(self.L_d/10):
                    free_vertices.append(val)
                    #free_vertex_index.append(bifurcation_free_vertex_index[ind])
                else:
                    pass
        elif inlet_outlet_index == 1:
            #Get positive z outlet free edge vertices
            for ind, val in enumerate(bifurcation_free_vertices):
                if val[0] > np.abs(self.L_d/10) and val[2] > 0:
                    free_vertices.append(val)
                    #free_vertex_index.append(bifurcation_free_vertex_index[ind])
                else:
                    pass
        elif inlet_outlet_index == 2:
            #Get negative z outlet free edge vertices
            for ind, val in enumerate(bifurcation_free_vertices):
                if val[0] > np.abs(self.L_d/10) and val[2] < 0:
                    free_vertices.append(val)
                    #free_vertex_index.append(bifurcation_free_vertex_index[ind])
                else:
                    pass
        free_vertices = np.array(free_vertices)
        #free_vertex_index = np.array(free_vertex_index)
        
        return free_vertices
        


    
    

