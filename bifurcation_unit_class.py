import numpy as np
import open3d as o3d
from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_mesh_class import bifurcation_mesh
from bifurcation_free_vertices_class import bifurcation_free_vertices

class bifurcation_unit:

    def __init__(self, parameters):
        self.R_p = parameters[0]; self.L_p = parameters[1]; self.R_d = parameters[2]
        self.L_d = parameters[3]; self.R_o = parameters[4]; self.R_i = parameters[5]
        self.iota_b = parameters[6]; self.delta_alpha = parameters[7]
        self.iota_gamma = parameters[8]; self.n_cont_outlets = parameters[9]
        self.parameters = parameters
        self.generate_bifurcation_unit()
        self.generate_initial_bifurcation_unit()
    
    def generate_bifurcation_unit(self):
        #Create Mesh
        mesh_class = bifurcation_mesh(self.parameters)
        self.bifurcation_unit_mesh = mesh_class.capped_mesh
        #Get Vectors 
        vector_class = bifurcation_vectors(self.parameters)
        self.r_pos = vector_class.r_pos
        self.r_neg = vector_class.r_neg
        self.n_pos = vector_class.n_pos
        self.n_neg = vector_class.n_neg
        #Get free vertices
        vertices_class = bifurcation_free_vertices(self.parameters, self.bifurcation_unit_mesh)
        self.inlet_free_vertices = vertices_class.inlet_free_vertices
        self.positive_outlet_free_vertices = vertices_class.positive_outlet_free_vertices
        self.negative_outlet_free_vertices = vertices_class.negative_outlet_free_vertices
        #Get free edges
        self.inlet_free_edges = vertices_class.inlet_free_edges
        self.positive_outlet_free_edges = vertices_class.positive_outlet_free_edges
        self.negative_outlet_free_edges = vertices_class.negative_outlet_free_edges

    def generate_initial_bifurcation_unit(self):
        #Create Mesh
        mesh_class = bifurcation_mesh(self.parameters)
        self.initial_bifurcation_unit_mesh = mesh_class.initial_capped_mesh
        #Get Vectors 
        vector_class = bifurcation_vectors(self.parameters)
        self.r_pos = vector_class.r_pos
        self.r_neg = vector_class.r_neg
        self.n_pos = vector_class.n_pos
        self.n_neg = vector_class.n_neg
        #Get free vertices
        vertices_class = bifurcation_free_vertices(self.parameters, self.bifurcation_unit_mesh)
        self.inlet_free_vertices = vertices_class.inlet_free_vertices
        self.positive_outlet_free_vertices = vertices_class.positive_outlet_free_vertices
        self.negative_outlet_free_vertices = vertices_class.negative_outlet_free_vertices
        #Get free edges
        self.inlet_free_edges = vertices_class.inlet_free_edges
        self.positive_outlet_free_edges = vertices_class.positive_outlet_free_edges
        self.negative_outlet_free_edges = vertices_class.negative_outlet_free_edges

    def rotate_bifurcation_unit(self, matrix):
        #Rotate Mesh
        self.bifurcation_unit_mesh.rotate(matrix, center = np.array(np.zeros(3)))
        #Rotate Vectors
        rotation_vector_class = bifurcation_vectors(self.parameters)
        rotation_vector_class.rotate_vector(matrix)
        self.r_pos = rotation_vector_class.r_pos
        self.r_neg = rotation_vector_class.r_neg
        self.n_pos = rotation_vector_class.n_pos
        self.n_neg = rotation_vector_class.n_neg
        #Rotate Vertices 
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

        return
    
    def translate_bifurcation_unit(self, translation_vector):
        #Translate Mesh
        self.bifurcation_unit_mesh.translate(translation_vector)
        #Translate Position Vectors
        translation_vector_class = bifurcation_vectors(self.parameters)
        translation_vector_class.translate_vector(translation_vector)
        self.r_pos = translation_vector_class.r_pos
        self.r_neg = translation_vector_class.r_neg
        #Translate Vertices
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

