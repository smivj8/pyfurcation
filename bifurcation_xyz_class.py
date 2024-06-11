import numpy as np
from stylianou_geometry_equations_new import *

class bifurcation_xyz:

    def __init__(self, parameters):
        self.parameters = parameters
        self.parent_radius = parameters[0]

