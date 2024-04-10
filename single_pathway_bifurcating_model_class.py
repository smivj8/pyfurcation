import numpy as np
from bifurcation_unit_class import bifurcation_unit
from utility_functions_bifurcation import *


class single_pathway_bifurcating_model:

    def __init__(self, tree_parameters, axial_rotation_angles):
        self.tree_parameters = tree_parameters
        self.n_generations = len(self.tree_parameters)
        self.axial_rotation_angles = axial_rotation_angles
        self.generate_single_pathway_bifurcating_tree()

    def generate_single_pathway_bifurcating_tree(self):
        bifurcation_unit_gen_0 = bifurcation_unit(self.tree_parameters[0])
        self.tree_mesh = bifurcation_unit_gen_0.bifurcation_unit_mesh
        self.cont_outlet_position = bifurcation_unit_gen_0.r_pos
        self.cont_outlet_orientation = bifurcation_unit_gen_0.n_pos
        self.cont_outlet_vertices = bifurcation_unit_gen_0.positive_outlet_free_vertices
        for k in range(1, self.n_generations):
            self.append_bifurcation_unit_to_tree(self.tree_parameters[k], self.axial_rotation_angles[k])

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
        #update position and orientation vectors
        self.cont_outlet_position = bifurcation_unit_gen_k.r_pos
        self.cont_outlet_orientation = bifurcation_unit_gen_k.n_pos
        #update vertices
        self.cont_outlet_vertices = bifurcation_unit_gen_k.positive_outlet_free_vertices

    


    