import numpy as np
import open3d as o3d
from bifurcation_vectors_class import bifurcation_vectors
from bifurcation_mesh_class import bifurcation_mesh
from bifurcation_free_vertices_class import bifurcation_free_vertices

class bifurcation_unit:

    def __init__(self, parameters, initial_bool):
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
        if not initial_bool:
            self.generate_bifurcation_unit()
        else:
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
        self.bifurcation_unit_mesh = mesh_class.initial_capped_mesh
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

