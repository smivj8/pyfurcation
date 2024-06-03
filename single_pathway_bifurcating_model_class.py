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
            vectors.rotate_vector(R1)
            vectors.rotate_vector(R2)
            #Next test: vectors.rotate_vector(R2 @ R1)
            vectors.translate_vector(tree_outlet_positions[j-1])
            #update new position and normal
            tree_outlet_positions = vectors.r_pos
            tree_outlet_normals = vectors.n_pos
        self.tree_outlet_positions = tree_outlet_positions
        self.tree_outlet_normals = tree_outlet_normals

    def generate_single_pathway_bifurcating_tree(self):
        bifurcation_unit_gen_0 = bifurcation_unit(self.tree_parameters[0])
        tree_mesh = bifurcation_unit_gen_0.bifurcation_unit_mesh
        


    def append_bifurcation_unit_to_tree(self, unit_parameters, rotation_angle):
        #Create bifurcation unit generation_k
        bifurcation_unit_gen_k = bifurcation_unit(unit_parameters)
        #Transform unit gen k to correct position and orientation (at continuation outlet)
        R1 = unit_vector_match_rotation(np.array([1,0,0]), self.cont_outlet_orientation)
        R2 = rotation_about_vector_u(self.cont_outlet_orientation, rotation_angle)
        bifurcation_unit_gen_k.rotate_bifurcation_unit(R2 @ R1)
        # bifurcation_unit_gen_k.rotate_bifurcation_unit(R2)
        bifurcation_unit_gen_k.translate_bifurcation_unit(self.cont_outlet_position)
        #Create junction mesh
        junction_mesh = create_junction_mesh(self.cont_outlet_vertices, bifurcation_unit_gen_k.inlet_free_vertices, self.cont_outlet_orientation)
        #update mesh
        self.tree_mesh += (junction_mesh + bifurcation_unit_gen_k.bifurcation_unit_mesh)
        vectors = bifurcation_vectors(unit_parameters)
        r_pos_unit = vectors.r_pos
        n_pos_unit = vectors.n_pos
        #update position and orientation vectors
        self.cont_outlet_position = (R2 @ (R1 @ r_pos_unit)) + self.cont_outlet_position
        self.cont_outlet_orientation = bifurcation_unit_gen_k.n_pos
        #update vertices
        self.cont_outlet_vertices = bifurcation_unit_gen_k.positive_outlet_free_vertices

    
    # def generate_single_pathway_bifurcating_tree(self):
    #     bifurcation_unit_gen_0 = bifurcation_unit(self.tree_parameters[0])
    #     self.tree_mesh = bifurcation_unit_gen_0.bifurcation_unit_mesh
    #     self.cont_outlet_position = bifurcation_unit_gen_0.r_pos
    #     self.cont_outlet_orientation = bifurcation_unit_gen_0.n_pos
    #     self.cont_outlet_vertices = bifurcation_unit_gen_0.positive_outlet_free_vertices
    #     for k in range(1, self.n_generations):
    #         #self.append_bifurcation_unit_to_tree(self.tree_parameters[k], self.axial_rotation_angles[k])
    #         #Create bifurcation unit generation_k
    #         bifurcation_unit_gen_k = bifurcation_unit(self.tree_parameters[k])
    #         #Transform unit gen k to correct position and orientation (at continuation outlet)
    #         R1 = unit_vector_match_rotation(np.array([1,0,0]), self.cont_outlet_orientation)
    #         R2 = rotation_about_vector_u(self.cont_outlet_orientation, self.axial_rotation_angles[k])
    #         bifurcation_unit_gen_k.rotate_bifurcation_unit(R2 @ R1)
    #         # bifurcation_unit_gen_k.rotate_bifurcation_unit(R2)
    #         bifurcation_unit_gen_k.translate_bifurcation_unit(self.cont_outlet_position)
    #         #Create junction mesh
    #         junction_mesh = create_junction_mesh(self.cont_outlet_vertices, bifurcation_unit_gen_k.inlet_free_vertices, self.cont_outlet_orientation)
    #         #update mesh
    #         self.tree_mesh += (junction_mesh + bifurcation_unit_gen_k.bifurcation_unit_mesh)
    #         #update position and orientation vectors
    #         self.cont_outlet_position = bifurcation_unit_gen_k.r_pos
    #         self.cont_outlet_orientation = bifurcation_unit_gen_k.n_pos
    #         #update vertices
    #         self.cont_outlet_vertices = bifurcation_unit_gen_k.positive_outlet_free_vertices

    