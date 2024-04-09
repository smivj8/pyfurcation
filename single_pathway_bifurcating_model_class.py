import numpy as np
from bifurcation_unit_class import bifurcation_unit


class single_pathway_bifurcating_model:

    def __init__(self, tree_parameters):
        self.tree_parameters = tree_parameters
        self.n_generations = len(self.tree_parameters)
    
    def initialize_tree(self):
        bifurcation_unit_gen_0 = bifurcation_unit(self.tree_parameters[0])
        self.tree_mesh = bifurcation_unit_gen_0.bifurcation_unit_mesh
        self.cont_outlet_position = bifurcation_unit_gen_0.r_pos
        self.cont_outlet_orientation = bifurcation_unit_gen_0.n_pos
        self.cont_outlet_vertices = bifurcation_unit_gen_0.positive_outlet_free_vertices