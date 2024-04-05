from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_mesh_class import bifurcation_mesh
import numpy as np
import open3d as o3d
from utility_functions_bifurcation import *
from stylianou_geometry_equations import *

class bifurcation_free_vertices:

    def __init__(self, parameters, truncated_mesh, inlet_outlet_index_in):
        self.R_p = parameters[0]; self.L_p = parameters[1]; self.R_d = parameters[2]
        self.L_d = parameters[3]; self.R_o = parameters[4]; self.R_i = parameters[5]
        self.iota_b = parameters[6]; self.delta_alpha = parameters[7]
        self.iota_gamma = parameters[8]; self.n_cont_outlets = parameters[9]
        self.parameters = parameters
        self.inlet_outlet_index = inlet_outlet_index_in
        self.truncated_mesh = truncated_mesh
        # inlet_outlet_index:
        #    0 -> inlet_vertices; 1 ->positive_outlet_vertices; 2 -> negative_outlet_vertices
        if self.inlet_outlet_index == 0:
            self.get_inlet_free_vertices()
        elif self.inlet_outlet_index == 1:
            self.get_positive_outlet_free_vertices()
        else:
            self.get_negative_outlet_free_vertices()
        

        #self.truncated_mesh = bifurcation_mesh(self.parameters).truncated_mesh


    def get_inlet_free_vertices(self):
        free_edges = self.truncated_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        bifurcation_free_vertex_index = np.unique(np.concatenate((numpy_free_edges[:,0], numpy_free_edges[:,1])))
        bifurcation_vertices = np.asarray(self.truncated_mesh.vertices)[bifurcation_free_vertex_index]
        free_vertices = []
        for ind, val in enumerate(bifurcation_vertices):
                if val[0] < np.abs(self.L_d/10):
                    free_vertices.append(val)
                    #free_vertex_index.append(bifurcation_free_vertex_index[ind])
                else:
                    pass
        free_vertices = np.array(free_vertices)
        self.inlet_free_vertices = free_vertices

        return self.inlet_free_vertices
    
    def get_positive_outlet_free_vertices(self):
        free_edges = self.truncated_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        bifurcation_free_vertex_index = np.unique(np.concatenate((numpy_free_edges[:,0], numpy_free_edges[:,1])))
        bifurcation_vertices = np.asarray(self.truncated_mesh.vertices)[bifurcation_free_vertex_index]
        free_vertices = []
        for ind, val in enumerate(bifurcation_vertices):
            if val[0] > np.abs(self.L_d/10) and val[2] > 0:
                free_vertices.append(val)
                #free_vertex_index.append(bifurcation_free_vertex_index[ind])
            else:
                pass
        free_vertices = np.array(free_vertices)
        self.positive_outlet_free_vertices = free_vertices

        return self.positive_outlet_free_vertices
    
    def get_negative_outlet_free_vertices(self):
        free_edges = self.truncated_mesh.get_non_manifold_edges(allow_boundary_edges = False)
        numpy_free_edges = np.asarray(free_edges)
        bifurcation_free_vertex_index = np.unique(np.concatenate((numpy_free_edges[:,0], numpy_free_edges[:,1])))
        bifurcation_vertices = np.asarray(self.truncated_mesh.vertices)[bifurcation_free_vertex_index]
        free_vertices = []
        for ind, val in enumerate(bifurcation_vertices):
            if val[0] > np.abs(self.L_d/10) and val[2] < 0:
                free_vertices.append(val)
                #free_vertex_index.append(bifurcation_free_vertex_index[ind])
            else:
                pass
        free_vertices = np.array(free_vertices)
        self.negative_outlet_free_vertices = free_vertices

        return self.negative_outlet_free_vertices
    
    def rotate_vertices(self, matrix):
        if self.inlet_outlet_index == 0:
            self.inlet_free_vertices = matrix @ self.inlet_free_vertices
            return self.inlet_free_vertices
        elif self.inlet_outlet_index == 1:
            self.positive_outlet_free_vertices = matrix @ self.positive_outlet_free_vertices
            return self.positive_outlet_free_vertices
        else:
            self.negative_outlet_free_vertices = matrix @ self.negative_outlet_free_vertices
            return self.negative_outlet_free_vertices
        

    def translate_vertices(self, translation_vector):
        if self.inlet_outlet_index == 0:
            self.inlet_free_vertices = self.inlet_free_vertices + translation_vector
            return self.inlet_free_vertices
        elif self.inlet_outlet_index == 1:
            self.positive_outlet_free_vertices = self.positive_outlet_free_vertices + translation_vector
            return self.positive_outlet_free_vertices
        else:
            self.negative_outlet_free_vertices = self.negative_outlet_free_vertices + translation_vector
            return self.negative_outlet_free_vertices
        