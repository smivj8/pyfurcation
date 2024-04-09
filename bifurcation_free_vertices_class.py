import numpy as np
import open3d as o3d
from utility_functions_bifurcation import *
from stylianou_geometry_equations import *

class bifurcation_free_vertices:

    def __init__(self, parameters, truncated_mesh):
        self.R_p = parameters[0]; self.L_p = parameters[1]; self.R_d = parameters[2]
        self.L_d = parameters[3]; self.R_o = parameters[4]; self.R_i = parameters[5]
        self.iota_b = parameters[6]; self.delta_alpha = parameters[7]
        self.iota_gamma = parameters[8]; self.n_cont_outlets = parameters[9]
        self.parameters = parameters
        self.truncated_mesh = truncated_mesh
        self.get_inlet_free_vertices()
        self.get_inlet_free_edges()
        self.get_positive_outlet_free_vertices()
        self.get_positive_outlet_free_edges()
        self.get_negative_outlet_free_vertices()
        self.get_negative_outlet_free_edges()
        

        #self.truncated_mesh = bifurcation_mesh(self.parameters).truncated_mesh


    def get_inlet_free_vertices(self):
        free_edges = self.truncated_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        bifurcation_free_vertex_index = np.unique(np.concatenate((numpy_free_edges[:,0], numpy_free_edges[:,1])))
        bifurcation_vertices = np.asarray(self.truncated_mesh.vertices)[bifurcation_free_vertex_index]
        inlet_free_vertices = []
        for ind, val in enumerate(bifurcation_vertices):
            if val[0] < np.abs(self.L_d/10):
                inlet_free_vertices.append(val)
                #free_vertex_index.append(bifurcation_free_vertex_index[ind])
            else:
                pass
        inlet_free_vertices = np.array(inlet_free_vertices)
        self.inlet_free_vertices = inlet_free_vertices
        return
    
    def get_inlet_free_edges(self):
        free_edges = self.truncated_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        mesh_vertices = np.asarray(self.truncated_mesh.vertices)
        inlet_free_edges = []
        for ind, val in enumerate(numpy_free_edges):
            if mesh_vertices[numpy_free_edges[ind, 0], 0] < np.abs(self.L_d/10):
                inlet_free_edges.append(val)
            else:
                pass
        inlet_free_edges = np.array(inlet_free_edges)
        self.inlet_free_edges = inlet_free_edges

    
    def get_positive_outlet_free_vertices(self):
        free_edges = self.truncated_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        bifurcation_free_vertex_index = np.unique(np.concatenate((numpy_free_edges[:,0], numpy_free_edges[:,1])))
        bifurcation_vertices = np.asarray(self.truncated_mesh.vertices)[bifurcation_free_vertex_index]
        positive_outlet_free_vertices = []
        for ind, val in enumerate(bifurcation_vertices):
            if val[0] > np.abs(self.L_d/10) and val[2] > 0:
                positive_outlet_free_vertices.append(val)
                #free_vertex_index.append(bifurcation_free_vertex_index[ind])
            else:
                pass
        positive_outlet_free_vertices = np.array(positive_outlet_free_vertices)
        self.positive_outlet_free_vertices = positive_outlet_free_vertices
        return
    
    def get_positive_outlet_free_edges(self):
        free_edges = self.truncated_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        mesh_vertices = np.asarray(self.truncated_mesh.vertices)
        positive_outlet_free_edges = []
        for ind, val in enumerate(numpy_free_edges):
            bool_1 = mesh_vertices[numpy_free_edges[ind, 0], 0] > np.abs(self.L_d/10)
            bool_2 = mesh_vertices[numpy_free_edges[ind, 0], 2] > 0
            if bool_1 and bool_2:
                positive_outlet_free_edges.append(val)
            else:
                pass
        positive_outlet_free_edges = np.array(positive_outlet_free_edges)
        self.positive_outlet_free_edges = positive_outlet_free_edges
    
    def get_negative_outlet_free_vertices(self):
        free_edges = self.truncated_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        bifurcation_free_vertex_index = np.unique(np.concatenate((numpy_free_edges[:,0], numpy_free_edges[:,1])))
        bifurcation_vertices = np.asarray(self.truncated_mesh.vertices)[bifurcation_free_vertex_index]
        negative_outlet_free_vertices = []
        for ind, val in enumerate(bifurcation_vertices):
            if val[0] > np.abs(self.L_d/10) and val[2] < 0:
                negative_outlet_free_vertices.append(val)
                #free_vertex_index.append(bifurcation_free_vertex_index[ind])
            else:
                pass
        negative_outlet_free_vertices = np.array(negative_outlet_free_vertices)
        self.negative_outlet_free_vertices = negative_outlet_free_vertices
        return
    
    def get_negative_outlet_free_edges(self):
        free_edges = self.truncated_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        mesh_vertices = np.asarray(self.truncated_mesh.vertices)
        negative_outlet_free_edges = []
        for ind, val in enumerate(numpy_free_edges):
            bool_1 = mesh_vertices[numpy_free_edges[ind, 0], 0] > np.abs(self.L_d/10)
            bool_2 = mesh_vertices[numpy_free_edges[ind, 0], 2] < 0
            if bool_1 and bool_2:
                negative_outlet_free_edges.append(val)
            else:
                pass
        negative_outlet_free_edges = np.array(negative_outlet_free_edges)
        self.negative_outlet_free_edges = negative_outlet_free_edges
    
    def rotate_vertices(self, matrix):
        for ind, vertex in enumerate(self.inlet_free_vertices):
            self.inlet_free_vertices[ind] = matrix @ self.inlet_free_vertices[ind]
        for ind, vertex in enumerate(self.positive_outlet_free_vertices):
            self.positive_outlet_free_vertices[ind] = matrix @ self.positive_outlet_free_vertices[ind]
        for ind, vertex in enumerate(self.negative_outlet_free_vertices):
            self.negative_outlet_free_vertices[ind] = matrix @ self.negative_outlet_free_vertices[ind]
        

    def translate_vertices(self, translation_vector):
        # for ind, vertex in enumerate(self.inlet_free_vertices):
        #     self.inlet_free_vertices[ind] = self.inlet_free_vertices[ind] + translation_vector
        # for ind, vertex in enumerate(self.positive_outlet_free_vertices):
        #     self.positive_outlet_free_vertices[ind] = self.positive_outlet_free_vertices[ind] + translation_vector
        # for ind, vertex in enumerate(self.negative_outlet_free_vertices):
        #     self.negative_outlet_free_vertices[ind] = self.negative_outlet_free_vertices[ind] + translation_vector
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
