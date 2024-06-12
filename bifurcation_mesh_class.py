from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_free_vertices_class import bifurcation_free_vertices
import numpy as np
import open3d as o3d
from utility_functions_bifurcation import *
from stylianou_geometry_equations_legacy import *


class bifurcation_mesh:

    def __init__(self, parameters):
        """
        A class for the generation of the triangle mesh for a SYMMETRIC bifurcation unit
        using a Poisson mesh reconstruction. Generates high density xyz point cloud 
        using bifurcation unit parameters from [1] and [2], converts to open3d pointcloud,
        generates triangle mesh from Poisson mesh reconstruction, crops excess triangles
        from Poisson reconstruction using clean crop, and truncates or closes daughter
        outlets according to number of continuation outlets for constructing a larger tree.

        [1] T. Heistracher and W. Hofmann, “Physiologically realistic models of bronchial airway 
        bifurcations,” Journal of Aerosol Science, vol. 26, no. 3, pp. 497–509, Apr. 1995, 
        doi: 10.1016/0021-8502(94)00113-D.

        [2] F. S. Stylianou, J. Sznitman, and S. C. Kassinos, “Direct numerical simulation of particle laden 
        flow in a human airway bifurcation model,” International Journal of Heat and Fluid Flow, vol. 61, pp. 
        677–710, Oct. 2016, doi: 10.1016/j.ijheatfluidflow.2016.07.013.

        Args:
            parameters (list): bifurcation unit parameters, see github README for more information
        
        """
        self.parameters = parameters
        self.parent_radius = parameters[0]
        #for the purpose of generating a bifurcation unit, the parent length is always set to 0 to
        #prevent doubling the length of a single generation. Then, since the poisson mesh method
        #rounds the edges slightly, the radius is slightly increased before cropping
        self.parent_length = 0.05*parameters[3]
        self.daughter_radius = parameters[2]
        #Similar to above, the parent length is slightly increased for cropping
        self.daughter_length = parameters[3] * 1.05
        self.outer_radius = parameters[4]
        self.branching_angle = parameters[5]
        self.carina_angle = parameters[6]
        self.n_cont_outlets = parameters[7]
        self.generate_parametric_xyz()
        self.generate_o3d_pointcloud()
        self.generate_o3d_mesh() 
        self.crop_bifurcation_mesh()
        self.truncate_continuation_outlets()
        self.cap_terminal_outlets()
        self.cap_initial_inlet()


    def generate_parametric_xyz(self, n_circ = 100, n_streamline = 120, 
                                n_car_range = 80, n_car_rad = 8, scalar = 1.1):
        """
        Generates high density numpy xyz pointcloud array using parametric equations 
        from [2]. To modify pointcloud density, parameters above can be manipulated.

        Args:
            n_circ (int, optional): number of points around the circumference of the 
                bifurcation unit. 
            n_streamline (int, optional): number of points along length of bifurcation
                unit.
            n_car_range (int, optional ): number of points along the arc of the carina
                connecting the daughter outlet pathways.
            n_car_rad (int, optional): number of points along the radius of the carina 
                connecting the daughter outlet pathways.
            sclar (float, optional): uhhh I honestly don't remember. This
                is used to calculate the normal vectors along the bifurcation unit,
                which might not even be used. Either way, just don't touch this unless you
                want to mess with a lot of other code
        """
        #Since the Poisson algorithm tends to round the edges of the meshes at the 
        #inlet or outlet, these parameters extend the cylinders at each end. The 
        #mesh is still cropped at the original place
        inlet_extension = 0.05*self.L_d
        outlet_extension = 1.05*self.L_d
        #See [2] or README for description of these variables. Basically creating the
        #parameterization variables for construction of the pointcloud
        carina_parameters = [self.R_p, self.R_d, self.R_o, self.R_i, self.branching_angle]
        phi_s_pos_min = np.arctan(-1*inlet_extension/self.R_o)
        eta = np.tan(self.iota_gamma - self.branching_angle)
        phi_s_pos_max = self.branching_angle + np.arctan(eta + (outlet_extension/self.R_o))
        phi_s_pos = np.linspace(phi_s_pos_min, phi_s_pos_max, n_streamline)
        phi_c = np.linspace(-1*np.pi/2, 3*np.pi/2, n_circ)
        phi_s_carina_max = self.iota_gamma
        phi_s_carina = np.linspace(0, phi_s_carina_max, n_car_range)
        #Initialization of pointcloud and normals numpy arrays
        xyz_pos = np.zeros([n_circ*n_streamline, 3])
        xyz_pos_normals = np.zeros([n_circ*n_streamline, 3])
        car_pos = np.zeros([n_car_range*n_car_rad, 3])
        car_pos_normals =  np.zeros([n_car_range*n_car_rad, 3])
        car_neg = np.zeros([n_car_range*n_car_rad, 3])
        car_neg_normals = np.zeros([n_car_range*n_car_rad, 3])
        #Populate pointcloud and normals arrays of the bifurcation unit        
        iteration = 0
        for i in range(len(phi_s_pos)):
            for j in range(len(phi_c)):
                xyz_pos[iteration,:], xyz_pos_normals[iteration,:] = r_positive(phi_s_pos[i], phi_c[j], 
                                                phi_c_pos_min, phi_c_pos_max, 
                                                self.parameters, scalar)
                iteration += 1
        #Populate pointcloud and normals array of the carina
        for m in range(len(phi_s_carina)):
            _, omega, _, _ = rho_omega_theta_parameters(phi_s_carina[m], self.parameters, carina_parameters, scalar)
            psi_range = np.linspace(-1*omega, omega, n_car_rad)
            for n in range(len(psi_range)):
                car_pos[n*n_car_range + m, :], car_pos_normals[n*n_car_range + m, :] = carina_positive(
                    phi_s_carina[m], psi_range[n], self.parameters, carina_parameters, scalar)
                car_neg[n*n_car_range + m, :], car_neg_normals[n*n_car_range + m, :] = carina_negative(
                    phi_s_carina[m], psi_range[n], self.parameters, carina_parameters, scalar)
        #combine pointclouds and normals into single array
        xyz_pc_positive = np.concatenate((xyz_pos, xyz_pos_normals), axis = 1)
        carina_pc_positive = np.concatenate((car_pos, car_pos_normals), axis = 1)
        carina_pc_negative = np.concatenate((car_neg, car_neg_normals), axis = 1)
        #remove NAN values from bifurcation pathway array. These are expected.
        xyz_pos_trimmed = xyz_pc_positive[~np.isnan(xyz_pc_positive).any(axis=1),:] 
        xyz_neg_trimmed = np.zeros(np.shape(xyz_pos_trimmed))
        #Copy and reflect the positive xyz values across the xy plane. As of now,
        #this script only supports SYMMETRIC bifurcation units... assymetry is a
        #problem for another intern
        xyz_neg_trimmed[:,0] = xyz_pos_trimmed[:,0]
        xyz_neg_trimmed[:,1] = xyz_pos_trimmed[:,1]
        xyz_neg_trimmed[:,2] = -1*xyz_pos_trimmed[:,2]
        xyz_neg_trimmed[:,3] = xyz_pos_trimmed[:,3]
        xyz_neg_trimmed[:,4] = xyz_pos_trimmed[:,4]
        xyz_neg_trimmed[:,5] = -1*xyz_pos_trimmed[:,5]
        #Remove NAN values from carina xyz pointcloud
        car_pos_trimmed = carina_pc_positive[~np.isnan(carina_pc_positive).any(axis=1),:]
        car_neg_trimmed = carina_pc_negative[~np.isnan(carina_pc_negative).any(axis=1),:]
        #combine all xyz pointclouds and normal vectors to single array
        xyz = np.concatenate((xyz_pos_trimmed, xyz_neg_trimmed, car_pos_trimmed, car_neg_trimmed))
        #Shift pointcloud to x: 0.
        #This places the center of the unit inlet at the cartesian origin, and is necessary for
        #all future mesh manipulations
        x_min = np.min(xyz[:,0])
        xyz[:,0] = xyz[:,0] - x_min - inlet_extension
        #set class global variable to be passed into next function

        
        self.xyz = xyz        
        return
    
    def generate_o3d_pointcloud(self):
        """
        Convert numpy arrays into open3d pointcloud. See open3d documentation
        for more info
        """
        #Initialize open3d pointcloud object
        pointcloud = o3d.geometry.PointCloud()
        if len(self.xyz[0]) == 6:
            #If normals are given, set pointcloud normals
            pointcloud.points = o3d.utility.Vector3dVector(self.xyz[:,:3])
            pointcloud.normals = o3d.utility.Vector3dVector(self.xyz[:,3:])
        else:
            #Else only set pointcloud
            pointcloud.points = o3d.utility.Vector3dVector(self.xyz)
        self.pointcloud = pointcloud
        return
    
    def generate_o3d_mesh(self, depth = 7, width = 0, scale = 1.1, 
                      linear_fit = False):
        """"
        Generate triangle mesh from open3d pointcloud using Poisson mesh reconstruction algorithm.
        See open3d documentation for more information.

        Args:
            depth (int, optional): Poisson mesh parameter. See open3d documentation.
            width (int, optional): Poisson mesh parameter. See open3d documentation.
            scale (float, optional): Poisson mesh parameter. See open3d documentation.
            linear_fit (bool, optional): Poisson mesh parameter. See open3d documentation.
        """
        poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(self.pointcloud, depth = depth, width = width,
                                                                                scale = scale, linear_fit = linear_fit)[0]
        self.mesh = poisson_mesh
        return
    
    def crop_bifurcation_mesh(self):
        """
        Method to crop bifurcation unit mesh at inlet and both outlets, normal to the 
        flow plane. Uses a clean crop function in file:utility_functions_bifurcation.py
        adapted from something found on stack overflow. Poisson mesh returns a closed 
        (watertight) mesh with "bubbles" around opening. This method removes those bubbles
        and leaves the inlets and outlets open and undeformed.
        """
        #Crop mesh at inlet
        cropped_mesh = clean_crop_x(self.mesh, min_x = 0, min_slice = True)
        #Crop mesh at outlets
        #get (positive) outlet position and normal orientation.
        vectors = bifurcation_vectors(self.parameters)
        r_pos = vectors.r_pos
        n_pos = vectors.n_pos
        #Get crop distance of the outlets
        cos_theta = np.dot(n_pos, r_pos)/np.linalg.norm(r_pos)
        hypot_max = np.linalg.norm(r_pos)*cos_theta
        #rotate mesh and do the first outlet mesh crop. The clean crop
        #function ONLY works along the principle x-axis, so the entire mesh must be
        #rotated such that the outlet flow plane is normal to the x-axis.
        cropped_mesh.rotate((rotation_matrix_about_y(self.branching_angle)), center = np.zeros(3))
        cropped_mesh = clean_crop_x(cropped_mesh, max_x = hypot_max, max_slice = True)
        #rotate mesh (twice) and crop the second outlet mesh. 
        cropped_mesh.rotate((rotation_matrix_about_y(-1*self.branching_angle)), center = np.zeros(3))
        cropped_mesh.rotate((rotation_matrix_about_y(-1*self.branching_angle)), center = np.zeros(3))
        cropped_mesh = clean_crop_x(cropped_mesh, max_x = hypot_max, max_slice = True)
        #return mesh to original orientation
        cropped_mesh.rotate((rotation_matrix_about_y(self.branching_angle)), center = np.zeros(3))
        self.cropped_mesh = cropped_mesh
        return
    
    def truncate_continuation_outlets(self):
        """
        Method to truncate the continuation outlet(s). If the outlet will connect
        and continue to another bifurcation unit (or the acinus model), the outlet pathway
        needs to be truncated to 0.97 of the original length, to allow for the junction
        mesh to be created at a later step. This truncated distance is returned as part of 
        the junction mesh, so the final geometry dimensions are NOT affected.
        """
        #Get outlet locations and plane normals
        vectors = bifurcation_vectors(self.parameters)
        r_pos = vectors.r_pos
        n_pos = vectors.n_pos
        #get crop distance (0.97 of outlet pathway length)
        cos_theta = np.dot(n_pos, r_pos)/np.linalg.norm(r_pos)
        hypot_max = np.linalg.norm(r_pos)*cos_theta
        truncation_distance = hypot_max*0.97
        truncated_mesh = self.cropped_mesh
        #Truncate outlets
        if self.n_cont_outlets == 1:
            #Truncate ONLY positive z outlet for continuation outlet. other outlet will be capped.
            #Crop method is same as above.
            truncated_mesh.rotate((rotation_matrix_about_y(self.branching_angle)), center = np.zeros(3))
            truncated_mesh = clean_crop_x(truncated_mesh, max_x = truncation_distance, max_slice = True)
            truncated_mesh.rotate((rotation_matrix_about_y(-1*self.branching_angle)), center = np.zeros(3))
        elif self.n_cont_outlets == 2:
            #Truncate both outlets for continuation.
            truncated_mesh.rotate((rotation_matrix_about_y(self.branching_angle)), center = np.zeros(3))
            truncated_mesh = clean_crop_x(truncated_mesh, max_x = truncation_distance, max_slice = True)
            truncated_mesh.rotate((rotation_matrix_about_y(-1*self.branching_angle)), center = np.zeros(3))
            truncated_mesh.rotate((rotation_matrix_about_y(-1*self.branching_angle)), center = np.zeros(3))
            truncated_mesh = clean_crop_x(truncated_mesh, max_x = truncation_distance, max_slice = True)
            truncated_mesh.rotate((rotation_matrix_about_y(self.branching_angle)), center = np.zeros(3))
        else:
            #Truncate neither outlet, both will be capped off.
            pass
        self.truncated_mesh = truncated_mesh
        return self.truncated_mesh

        
    def cap_terminal_outlets(self):
        """
        Method for capping terminal outlets. If the outlet will NOT connect to another
        unit or acinus, close the outlet (untruncated) to finish with watertight mesh.
        """
        #Get vectors and free vertices of the mesh
        vectors = bifurcation_vectors(self.parameters)
        free_vertices = bifurcation_free_vertices(self.parameters, self.truncated_mesh)
        if self.n_cont_outlets == 1:
            #cap ONLY the negative-z outlet
            r_neg = vectors.r_neg
            vertices = np.asarray(self.truncated_mesh.vertices)
            vertices = np.vstack([vertices, r_neg])
            #capping is done by getting free edges on the outlet and connecting
            #each of the vertices in an edge to a new vertex at the center of the outlet. 
            #The resulting triangle is appended to the mesh
            negative_outlet_free_edges = free_vertices.negative_outlet_free_edges
            new_triangles = np.zeros((negative_outlet_free_edges.shape[0], 3))
            new_triangles[:,0] = negative_outlet_free_edges[:,0]
            new_triangles[:,1] = negative_outlet_free_edges[:,1]
            new_triangles[:,2] = np.shape(vertices)[0] - 1
            # Append the new triangles to the mesh
            triangles = np.asarray(self.truncated_mesh.triangles)
            triangles = np.vstack([triangles, new_triangles])
            #Create new open3d mesh object from triangles and vertices. Compute normals
            #and orient mesh (this is necessary for mesh analysis, though I'm not exactly
            #sure what is does)
            capped_mesh = o3d.geometry.TriangleMesh()
            capped_mesh.vertices = o3d.utility.Vector3dVector(np.array(vertices))
            capped_mesh.triangles = o3d.utility.Vector3iVector(np.array(triangles))
            capped_mesh.compute_triangle_normals()
            capped_mesh.orient_triangles()
            self.capped_mesh = capped_mesh
        elif self.n_cont_outlets == 0:
            #cap both outlets, method same as above
            r_pos = vectors.r_pos
            r_neg = vectors.r_neg
            vertices = np.asarray(self.truncated_mesh.vertices)
            vertices = np.vstack([vertices, r_pos, r_neg])
            positive_outlet_free_edges = free_vertices.positive_outlet_free_edges
            positive_new_triangles = np.zeros((positive_outlet_free_edges.shape[0], 3))
            positive_new_triangles[:,0] = positive_outlet_free_edges[:,0]
            positive_new_triangles[:,1] = positive_outlet_free_edges[:,1]
            positive_new_triangles[:,2] = np.shape(vertices)[0] - 2
            negative_outlet_free_edges = free_vertices.negative_outlet_free_edges
            negative_new_triangles = np.zeros((negative_outlet_free_edges.shape[0], 3))
            negative_new_triangles[:,0] = negative_outlet_free_edges[:,0]
            negative_new_triangles[:,1] = negative_outlet_free_edges[:,1]
            negative_new_triangles[:,2] = np.shape(vertices)[0] - 1
            triangles = np.asarray(self.truncated_mesh.triangles)
            triangles = np.vstack([triangles, positive_new_triangles, negative_new_triangles])
            capped_mesh = o3d.geometry.TriangleMesh()
            capped_mesh.vertices = o3d.utility.Vector3dVector(np.array(vertices))
            capped_mesh.triangles = o3d.utility.Vector3iVector(np.array(triangles))
            capped_mesh.compute_triangle_normals()
            capped_mesh.orient_triangles()
            self.capped_mesh = capped_mesh
        else:
            #cap neither outlet (they should already both be truncated, then)
            capped_mesh = self.truncated_mesh
            capped_mesh.compute_triangle_normals()
            capped_mesh.orient_triangles()
            self.capped_mesh = capped_mesh
        return
    
    def cap_initial_inlet(self):
        """
        Method for capping the inlet of the initial bifurcation unit mesh for an overall
        watertight mesh. Method is essentially the same as above
        """
        free_vertices = bifurcation_free_vertices(self.parameters, self.truncated_mesh)
        vertices = np.asarray(self.capped_mesh.vertices)
        vertices = np.vstack([vertices, np.zeros(3)])
        inlet_free_edges = free_vertices.inlet_free_edges
        new_triangles = np.zeros((inlet_free_edges.shape[0], 3))
        new_triangles[:,0] = inlet_free_edges[:,0]
        new_triangles[:,1] = inlet_free_edges[:,1]
        new_triangles[:,2] = np.shape(vertices)[0] - 1
        triangles = np.asarray(self.capped_mesh.triangles)
        triangles = np.vstack([triangles, new_triangles])
        initial_capped_mesh = o3d.geometry.TriangleMesh()
        initial_capped_mesh.vertices = o3d.utility.Vector3dVector(np.array(vertices))
        initial_capped_mesh.triangles = o3d.utility.Vector3iVector(np.array(triangles))
        initial_capped_mesh.compute_triangle_normals()
        initial_capped_mesh.orient_triangles()
        self.initial_capped_mesh = initial_capped_mesh
    
    def perform_mesh_analysis(self):
        """
        Method that returns analysis of the mesh using open3d defined functions.
        Returns:

        """
        analysis_mesh = self.initial_capped_mesh
        if self.n_cont_outlets != 0:
            analysis = [None, None]
        else:
            volume = analysis_mesh.get_volume()
            surface_area = analysis_mesh.get_surface_area()
            analysis = [volume, surface_area]
        return analysis
    