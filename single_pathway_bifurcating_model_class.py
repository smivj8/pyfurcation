import numpy as np
from bifurcation_unit_class import bifurcation_unit
from bifurcation_vectors_class import bifurcation_vectors
from utility_functions_bifurcation import *


class single_pathway_bifurcating_model:

    def __init__(self, tree_parameters, axial_rotation_angles):
        self.tree_parameters = tree_parameters
        self.n_generations = len(self.tree_parameters)
        self.axial_rotation_angles = axial_rotation_angles
        self.generate_single_pathway_tree_skeleton()
        self.generate_single_pathway_bifurcating_tree()
        

    def generate_single_pathway_tree_skeleton(self):
        tree_outlet_positions = np.zeros((self.n_generations, 3))
        tree_outlet_normals = np.zeros((self.n_generations, 3))
        vectors = bifurcation_vectors(self.tree_parameters[0])
        tree_outlet_positions[0] = vectors.r_pos
        tree_outlet_normals[0] = vectors.n_pos
        for j in range(1, self.n_generations):
            vectors = bifurcation_vectors(self.tree_parameters[j])
            R1 = unit_vector_match_rotation(np.array([1,0,0]), tree_outlet_normals[j-1])
            R2 = rotation_about_vector_u(tree_outlet_normals[j-1], self.axial_rotation_angles[j])
            #Should these be together or separate, hmmm....
            vectors.rotate_vector(R2@R1)
            vectors.translate_vector(tree_outlet_positions[j-1])
            #update new position and normal
            tree_outlet_positions[j] = vectors.r_pos
            tree_outlet_normals[j] = vectors.n_pos
        self.tree_outlet_positions = tree_outlet_positions
        self.tree_outlet_normals = tree_outlet_normals
        

    def generate_single_pathway_bifurcating_tree(self):
        #Initialize mesh and outlet vertices
        #TODO: Figure out how to implement the capping of the initial mesh
        bifurcation_unit_gen_0 = bifurcation_unit(self.tree_parameters[0])
        tree_mesh = bifurcation_unit_gen_0.initial_bifurcation_unit_mesh
        cont_outlet_vertices = bifurcation_unit_gen_0.positive_outlet_free_vertices
        #Generate rest of tree
        for k in range(1, self.n_generations):
            #initialize bifurcation unit and transform
            bifurcation_unit_gen_k = bifurcation_unit(self.tree_parameters[k])
            R1 = unit_vector_match_rotation(np.array([1,0,0]), self.tree_outlet_normals[k-1])
            R2 = rotation_about_vector_u(self.tree_outlet_normals[k-1], self.axial_rotation_angles[k])
            bifurcation_unit_gen_k.rotate_bifurcation_unit(R2@R1)
            bifurcation_unit_gen_k.translate_bifurcation_unit(self.tree_outlet_positions[k-1])
            #get inlet vertices and combine 
            unit_inlet_vertices = bifurcation_unit_gen_k.inlet_free_vertices
            #TODO Everything BUT the junction meshing is now working. It works for the first junction but not afterwards...
            #   NOTE: It is NOT fixed by setting the axial rotations to zero, so the issue is somewhere else, in the ordering
            #         vertices function maybe?    
            junction_mesh = create_junction_mesh(cont_outlet_vertices, unit_inlet_vertices, self.tree_outlet_normals[k-1], self.tree_outlet_positions[k-1])
            tree_mesh += (bifurcation_unit_gen_k.bifurcation_unit_mesh + junction_mesh)
            cont_outlet_vertices = bifurcation_unit_gen_k.positive_outlet_free_vertices
        self.tree_mesh = tree_mesh
        self.cont_outlet_vertices = cont_outlet_vertices
