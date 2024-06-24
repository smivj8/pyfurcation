import numpy as np
import open3d as o3d
from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_free_vertices_class import bifurcation_free_vertices
import utility_functions_bifurcation as ufb
import stylianou_geometry_equations as sge
import pointcloud_construction_functions as pcf


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

        [2] F. S. Stylianou, J. Sznitman, and S. C. Kassinos, “Direct numerical simulation of
        particle laden flow in a human airway bifurcation model,” International Journal of Heat and 
        Fluid Flow, vol. 61, pp. 677–710, Oct. 2016, doi: 10.1016/j.ijheatfluidflow.2016.07.013.

        Args:
            parameters (list): bifurcation unit parameters, see github README for more information
        
        """
        self.parameters = parameters
        self.parent_radius = parameters[0]
        #for the purpose of generating a bifurcation unit, the parent length is always set to 0 to
        #prevent doubling the length of a single generation. Then, since the poisson mesh method
        #rounds the edges slightly, the radius is slightly increased before cropping
        self.parent_length = parameters[3] * 0.05
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


    def generate_parametric_xyz(self, n_streamline = 120, n_circ = 120,
                                n_carina_range = 60, n_carina_rad = 8):
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
        """
        pipe_section_pcd, pipe_normals = \
            pcf.construct_pipe_section_pointcloud(self.parameters, n_streamline, n_circ)
        carina_pcd, car_normals = \
            pcf.construct_carina_pointcloud(self.parameters, n_carina_range, n_carina_rad)

        self.bifurcation_unit_xyz = np.vstack([pipe_section_pcd, carina_pcd])
        self.bifurcation_unit_normals = np.vstack([pipe_normals, car_normals])

    def generate_o3d_pointcloud(self):
        """
        Convert numpy arrays into open3d pointcloud. See open3d documentation
        for more info
        """
        #Initialize open3d pointcloud object
        pointcloud = o3d.geometry.PointCloud()
        pointcloud.points = o3d.utility.Vector3dVector(self.bifurcation_unit_xyz)
        pointcloud.normals = o3d.utility.Vector3dVector(self.bifurcation_unit_normals)
        self.pointcloud = pointcloud

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

    def crop_bifurcation_mesh(self):
        """
        Method to crop bifurcation unit mesh at inlet and both outlets, normal to the 
        flow plane. Uses a clean crop function in file:utility_functions_bifurcation.py
        adapted from something found on stack overflow. Poisson mesh returns a closed 
        (watertight) mesh with "bubbles" around opening. This method removes those bubbles
        and leaves the inlets and outlets open and undeformed.
        """
        #Crop mesh at inlet
        cropped_mesh = ufb.clean_crop_x(self.mesh, min_x = 0, min_slice = True)
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
        cropped_mesh.rotate((ufb.rotation_matrix_about_y(self.branching_angle)), center = np.zeros(3))
        cropped_mesh = ufb.clean_crop_x(cropped_mesh, max_x = hypot_max, max_slice = True)
        #rotate mesh (twice) and crop the second outlet mesh. 
        cropped_mesh.rotate((ufb.rotation_matrix_about_y(-1*self.branching_angle)), center = np.zeros(3))
        cropped_mesh.rotate((ufb.rotation_matrix_about_y(-1*self.branching_angle)), center = np.zeros(3))
        cropped_mesh = ufb.clean_crop_x(cropped_mesh, max_x = hypot_max, max_slice = True)
        #return mesh to original orientation
        cropped_mesh.rotate((ufb.rotation_matrix_about_y(self.branching_angle)), center = np.zeros(3))
        self.cropped_mesh = cropped_mesh
    
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
            truncated_mesh.rotate((ufb.rotation_matrix_about_y(self.branching_angle)), center = np.zeros(3))
            truncated_mesh = ufb.clean_crop_x(truncated_mesh, max_x = truncation_distance, max_slice = True)
            truncated_mesh.rotate((ufb.rotation_matrix_about_y(-1*self.branching_angle)), center = np.zeros(3))
        elif self.n_cont_outlets == 2:
            #Truncate both outlets for continuation.
            truncated_mesh.rotate((ufb.rotation_matrix_about_y(self.branching_angle)), center = np.zeros(3))
            truncated_mesh = ufb.clean_crop_x(truncated_mesh, max_x = truncation_distance, max_slice = True)
            truncated_mesh.rotate((ufb.rotation_matrix_about_y(-1*self.branching_angle)), center = np.zeros(3))
            truncated_mesh.rotate((ufb.rotation_matrix_about_y(-1*self.branching_angle)), center = np.zeros(3))
            truncated_mesh = ufb.clean_crop_x(truncated_mesh, max_x = truncation_distance, max_slice = True)
            truncated_mesh.rotate((ufb.rotation_matrix_about_y(self.branching_angle)), center = np.zeros(3))
        else:
            #Truncate neither outlet, both will be capped off.
            pass
        self.truncated_mesh = truncated_mesh
        
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
        elif analysis_mesh.is_watertight():
            volume = analysis_mesh.get_volume()
            surface_area = analysis_mesh.get_surface_area()
            analysis = [volume, surface_area]
        else:
            analysis = [None, None]
        return analysis
    