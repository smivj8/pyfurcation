import numpy as np
import open3d as o3d
class bifurcation_free_vertices:

    def __init__(self, parameters, vertices_mesh):
        """
        A class for the determination and manipulation of the free vertices and
        edges of a bifurcation unit mesh. This class requires both the mesh parameters
        and the mesh itself (the vertices are not necessarily dependent on the 
        parameters). 

        Args:
            parameters (list): bifurcation unit parameters, see github README for more 
                information

            truncated_mesh (open3d.mesh object): bifurcation unit mesh, with any required
                continuation outlets truncated. 
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
        self.parameters = parameters
        self.vertices_mesh = vertices_mesh
        self.get_inlet_free_vertices()
        self.get_inlet_free_edges()
        self.get_positive_outlet_free_vertices()
        self.get_positive_outlet_free_edges()
        self.get_negative_outlet_free_vertices()
        self.get_negative_outlet_free_edges()
        

        #self.truncated_mesh = bifurcation_mesh(self.parameters).truncated_mesh


    def get_inlet_free_vertices(self):
        """
        Method for getting the free vertices of the inlet of the bifurcation unit 
        mesh.
        """
        free_edges = self.vertices_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        bifurcation_free_vertex_index = np.unique(np.concatenate((numpy_free_edges[:,0], numpy_free_edges[:,1])))
        bifurcation_vertices = np.asarray(self.vertices_mesh.vertices)[bifurcation_free_vertex_index]
        inlet_free_vertices = []
        for ind, val in enumerate(bifurcation_vertices):
            if val[0] < np.abs(self.daughter_length/10):
                inlet_free_vertices.append(val)
                #free_vertex_index.append(bifurcation_free_vertex_index[ind])
            else:
                pass
        inlet_free_vertices = np.array(inlet_free_vertices)
        self.inlet_free_vertices = inlet_free_vertices
        return
    
    def get_inlet_free_edges(self):
        """
        Method for getting the free esges of the inlet of the bifurcation unit 
        mesh.
        """
        free_edges = self.vertices_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        mesh_vertices = np.asarray(self.vertices_mesh.vertices)
        inlet_free_edges = []
        for ind, val in enumerate(numpy_free_edges):
            if mesh_vertices[numpy_free_edges[ind, 0], 0] < np.abs(self.daughter_length/10):
                inlet_free_edges.append(val)
            else:
                pass
        inlet_free_edges = np.array(inlet_free_edges)
        self.inlet_free_edges = inlet_free_edges

    
    def get_positive_outlet_free_vertices(self):
        """
        Method for getting the free vertices of the positive z-axis outlet
        of the bifurcation unit mesh.
        """
        free_edges = self.vertices_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        bifurcation_free_vertex_index = np.unique(np.concatenate((numpy_free_edges[:,0], numpy_free_edges[:,1])))
        bifurcation_vertices = np.asarray(self.vertices_mesh.vertices)[bifurcation_free_vertex_index]
        positive_outlet_free_vertices = []
        for ind, val in enumerate(bifurcation_vertices):
            if val[0] > np.abs(self.daughter_length/10) and val[2] > 0:
                positive_outlet_free_vertices.append(val)
                #free_vertex_index.append(bifurcation_free_vertex_index[ind])
            else:
                pass
        positive_outlet_free_vertices = np.array(positive_outlet_free_vertices)
        self.positive_outlet_free_vertices = positive_outlet_free_vertices
        return
    
    def get_positive_outlet_free_edges(self):
        """
        Method for getting the free edges of the positive z-axis outlet
        of the bifurcation unit mesh.
        """
        free_edges = self.vertices_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        mesh_vertices = np.asarray(self.vertices_mesh.vertices)
        positive_outlet_free_edges = []
        for ind, val in enumerate(numpy_free_edges):
            bool_1 = mesh_vertices[numpy_free_edges[ind, 0], 0] > np.abs(self.daughter_length/10)
            bool_2 = mesh_vertices[numpy_free_edges[ind, 0], 2] > 0
            if bool_1 and bool_2:
                positive_outlet_free_edges.append(val)
            else:
                pass
        positive_outlet_free_edges = np.array(positive_outlet_free_edges)
        self.positive_outlet_free_edges = positive_outlet_free_edges
    
    def get_negative_outlet_free_vertices(self):
        """
        Method for getting the free vertices of the negative z-axis outlet
        of the bifurcation unit mesh.
        """
        free_edges = self.vertices_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        bifurcation_free_vertex_index = np.unique(np.concatenate((numpy_free_edges[:,0], numpy_free_edges[:,1])))
        bifurcation_vertices = np.asarray(self.vertices_mesh.vertices)[bifurcation_free_vertex_index]
        negative_outlet_free_vertices = []
        for ind, val in enumerate(bifurcation_vertices):
            if val[0] > np.abs(self.daughter_length/10) and val[2] < 0:
                negative_outlet_free_vertices.append(val)
            else:
                pass
        negative_outlet_free_vertices = np.array(negative_outlet_free_vertices)
        self.negative_outlet_free_vertices = negative_outlet_free_vertices
        return
    
    def get_negative_outlet_free_edges(self):
        """
        Method for getting the free edges of the negative z-axis outlet
        of the bifurcation unit mesh.
        """
        free_edges = self.vertices_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        mesh_vertices = np.asarray(self.vertices_mesh.vertices)
        negative_outlet_free_edges = []
        for ind, val in enumerate(numpy_free_edges):
            bool_1 = mesh_vertices[numpy_free_edges[ind, 0], 0] > np.abs(self.daughter_length/10)
            bool_2 = mesh_vertices[numpy_free_edges[ind, 0], 2] < 0
            if bool_1 and bool_2:
                negative_outlet_free_edges.append(val)
            else:
                pass
        negative_outlet_free_edges = np.array(negative_outlet_free_edges)
        self.negative_outlet_free_edges = negative_outlet_free_edges
    
    def rotate_vertices(self, matrix):
        """
        Rotates free vertices about origin according to a rotation matrix
        """
        for ind, vertex in enumerate(self.inlet_free_vertices):
            self.inlet_free_vertices[ind] = matrix @ self.inlet_free_vertices[ind]
        if self.n_cont_outlets == 2:
            for ind, vertex in enumerate(self.positive_outlet_free_vertices):
                self.positive_outlet_free_vertices[ind] = matrix @ self.positive_outlet_free_vertices[ind]
            for ind, vertex in enumerate(self.negative_outlet_free_vertices):
                self.negative_outlet_free_vertices[ind] = matrix @ self.negative_outlet_free_vertices[ind]
        elif self.n_cont_outlets == 1:
            for ind, vertex in enumerate(self.positive_outlet_free_vertices):
                self.positive_outlet_free_vertices[ind] = matrix @ self.positive_outlet_free_vertices[ind]
        else:
            pass

    def translate_vertices(self, translation_vector):
        """
        Translates free vertices according to a translation vector.
        """
        self.inlet_free_vertices = self.inlet_free_vertices + translation_vector
        if self.n_cont_outlets == 2:
            self.positive_outlet_free_vertices = self.positive_outlet_free_vertices + translation_vector
            self.negative_outlet_free_vertices = self.negative_outlet_free_vertices + translation_vector
        elif self.n_cont_outlets == 1:
            self.positive_outlet_free_vertices = self.positive_outlet_free_vertices + translation_vector
            self.negative_outlet_free_vertices = self.negative_outlet_free_vertices
        else:
            self.positive_outlet_free_vertices = self.positive_outlet_free_vertices
            self.negative_outlet_free_vertices = self.negative_outlet_free_vertices
