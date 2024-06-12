"""
A set of functions for defining the xyz pointcloud of a bifurcation unit with the 
modified carina according to [1] and [2]. All equations are from Appendix A - 
Mathematical Description of the Bifurcation Geometry of [2], and are labeled
A.# in accordance to the source.

Construction of bifurcation unit is dependent on a list of parameters. The parameters
argument completely defines the dimensions of a bifurcation unit. The parameters argument is
a list structured as follows: 

parameters = [parent_radius, parent_length, daughter_radius, daughter_length, 
              outer_radius, branching_angle, carina_angle, n_cont_outlets]

See [2] for a full description of the parameter descriptions.

For the purpose of creating a multiple generation branch, the parent length is set
to parent_length = 0 by default. This is because the parent length of a bifurcation
unit in generation k equals the daughter length of a bifurcation unit in generation 
k-1.

n_cont_outlets is the number of daughter outlets (0, 1, or 2) that will continue
on and attach to a successive bifurcation unit, and is not used in the creation of
the bifurcation unit pointcloud. As such, it is not accessed in this script.

[1] T. Heistracher and W. Hofmann, “Physiologically realistic models of bronchial airway 
bifurcations,” Journal of Aerosol Science, vol. 26, no. 3, pp. 497–509, Apr. 1995, 
doi: 10.1016/0021-8502(94)00113-D.

[2] F. S. Stylianou, J. Sznitman, and S. C. Kassinos, “Direct numerical simulation of particle
laden flow in a human airway bifurcation model,” International Journal of Heat and Fluid Flow, 
vol. 61, pp. 677–710, Oct. 2016, doi: 10.1016/j.ijheatfluidflow.2016.07.013.
"""
import numpy as np

####################################################################################################
#                                                                                                  #
#               Functions describing the xyz coordinates of the parent, merging,                   #
#                       and daughter pipe sections of the bifurcation unit.                        #
#                                                                                                  #
####################################################################################################

def transition_radius(phi_s, parameters):
    """
    Since the parent pipe section and daughter pipe section 
    have different radii, a smooth transition between the two is defined
    by the function relation given by Equation (A.10), pg. 700.

    Args:
        phi_s (float): streamline parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    parent_radius = parameters[0]
    daughter_radius = parameters[2]
    branching_angle = parameters[5]
    #Define sigmoidal function
    def sigmoidal(t_val, t_min, t_max):
        """
        A polynomial sigmoidal function to define a smooth transition between
        two parameteric values. A general definition of the equation
        defined in Equation (A.11), pg. 700.

        Args:
            t_val (float): the parametric value

            t_min (float): the minimum value of the parameter

            t_max (float): the maximum value of the parameter
        """
        if t_val < t_min:
            value = 0
        elif t_val >= t_min and phi_s <= t_max:
            value = 3 * ((t_val/(t_max - t_min))**2) - 2 * ((t_val/(t_max - t_min))**3)
        else:
            value = 1
        return value
    #calculate transition radius length at phi_s
    transition_radius_value = parent_radius - (parent_radius - daughter_radius) * sigmoidal(phi_s, 0, branching_angle)
    return transition_radius_value

def get_parent_pipe(phi_s, phi_c, parameters):
    """
    Construct the xyz position of the parent pipe section at parameteric
    location (phi_s, phi_c) for given parameters. 
    Equations (A.1), (A.2), (A.3), pg. 699

    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    parent_radius = parameters[0]
    parent_length = parameters[1]
    outer_radius = parameters[4]
    #Calculate x,y,z coordinates, pack into numpy array
    x_p = parent_length + outer_radius * np.tan(phi_s) #(A.1)
    y_p = parent_radius * np.cos(phi_c) #(A.2)
    z_p = parent_radius * np.sin(phi_c) #(A.3)
    xyz_p = np.array([x_p, y_p, z_p])
    return xyz_p

def get_merging_pipe_positive_z(phi_s, phi_c, parameters):
    """
    Construct the xyz position of the merging pipe section in the positive-z 
    range at parameteric location (phi_s, phi_c) for given parameters. 
    Positive-z variation of Equations (A.4), (A.5), (A.6), pg. 699

    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    parent_length = parameters[1]
    outer_radius = parameters[4]
    #Calculate xyz coordinates and pack into numpy array
    val = transition_radius(phi_s, parameters) * np.sin(phi_c)
    x_m_pos = parent_length + (outer_radius * np.sin(phi_s)) - (np.sin(phi_s) * val) #(A.4)
    y_m_pos = transition_radius(phi_s, parameters) * np.cos(phi_c) #(A.5)
    z_m_pos = outer_radius * (1-np.cos(phi_s)) + np.cos(phi_s) * val #(A.6)
    xyz_m_pos = np.array([x_m_pos, y_m_pos, z_m_pos])
    return xyz_m_pos

def get_merging_pipe_negative_z(phi_s, phi_c, parameters):
    """
    Construct the xyz position of the merging pipe section in the negative-z 
    range at parameteric location (phi_s, phi_c) for given parameters. 
    Negative-z variation of Equations (A.4), (A.5), (A.6), pg. 699

    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    parent_length = parameters[1]
    outer_radius = parameters[4]
    #Calculate xyz coordinates and pack into numpy array
    val = transition_radius(phi_s, parameters) * np.sin(phi_c)
    x_m_neg = parent_length + (outer_radius * np.sin(phi_s)) + (np.sin(phi_s) * val) #(A.4)
    y_m_neg = transition_radius(phi_s, parameters) * np.cos(phi_c) #(A.5)
    z_m_neg = -1 * outer_radius * (1-np.cos(phi_s)) + np.cos(phi_s) * val #(A.6)
    xyz_m_neg = np.array([x_m_neg, y_m_neg, z_m_neg])
    return xyz_m_neg

def get_daughter_pipe_positive_z(phi_s, phi_c, parameters):
    """
    Construct the xyz position of the daughter pipe section in the positive-z 
    range at parameteric location (phi_s, phi_c) for given parameters. 
    Positive-z variation of Equations (A.12), (A.13), (A.14), pg. 699

    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    outer_radius = parameters[4]
    branching_angle = parameters[5]
    #Calculate xyz coordinates and pack into numpy array
    xyz_m_pos = get_merging_pipe_positive_z(branching_angle, phi_c, parameters) #get xyz_m_pos coords
    zeta = outer_radius * np.tan(phi_s - branching_angle) #cartesian transform of phi_s
    x_d_pos = xyz_m_pos[0] + zeta*np.cos(branching_angle) #(A.12)
    y_d_pos = xyz_m_pos[1] #(A.13)
    z_d_pos = xyz_m_pos[2] + zeta*np.sin(branching_angle) #(A.14)
    xyz_d_pos = np.array([x_d_pos, y_d_pos, z_d_pos])
    return xyz_d_pos

def get_daughter_pipe_negative_z(phi_s, phi_c, parameters):
    """
    Construct the xyz position of the daughter pipe section in the negative-z 
    range at parameteric location (phi_s, phi_c) for given parameters. 
    Negative-z variation of Equations (A.12), (A.13), (A.14), pg. 699

    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    outer_radius = parameters[4]
    branching_angle = parameters[5]
    #Calculate xyz coordinates and pack into numpy array
    xyz_m_neg = get_merging_pipe_negative_z(branching_angle, phi_c, parameters) #get xyz_m_neg coords
    zeta = outer_radius * np.tan(phi_s - branching_angle) #cartesian transform of phi_s
    x_d_neg = xyz_m_neg[0] + zeta*np.cos(branching_angle) #(A.12)
    y_d_neg = xyz_m_neg[1] #(A.13)
    z_d_neg = xyz_m_neg[2] - zeta*np.sin(branching_angle) #(A.14)
    xyz_d_neg = np.array([x_d_neg, y_d_neg, z_d_neg])
    return xyz_d_neg

####################################################################################################
#                                                                                                  #
#          Helper Functions for determining the range of circumferential angles (phi_c)            #
#              for a given location along the streamline parameterization (phi_s)                  #
#                                                                                                  #
####################################################################################################

def get_position_alpha(parameters):
    """
    A method for getting the position vector of location ALPHA in Figure A.3, according to equation 
    (A.16). Equation (A.16) references equation (A.15) which, although defined in this code later,
    is simply the piecewise combinate of equations (A.1-6, 12-14) described above.

    Args:
        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    carina_angle = parameters[6]
    #calculcate position A
    position_alpha = get_parent_pipe(0, carina_angle/2, parameters) #(A.16)
    return position_alpha

def get_angle_phi_gamma(parameters):
    """
    Method for getting the angle phi_s corresponding to location GAMMA 
    (denoted PHI_GAMMA) in Figures A.2 and A.3. Determined using the relations found
    at the bottom of pg. 700 and top of pg. 701. The relation is not given an equation
    number. It is assumed that R_c = 0.4 * R_d, per Table 1 of [2]. I did the derivation
    of this by hand, so just ask me there are any questions

    Args:
        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    daughter_radius = parameters[2]
    outer_radius = parameters[4]
    branching_angle = parameters[5]
    # Calculate PHI_GAMMA
    value_1 = (1 + ((1.4 * daughter_radius)/outer_radius)) * (1/np.tan(branching_angle))
    value_2 = 1 / np.sin(branching_angle)
    delta = np.arctan(value_1 - value_2)
    phi_gamma = delta + branching_angle
    return phi_gamma

def get_position_gamma(parameters):
    """
    A method for getting the position vector of location GAMMA in Figures A.2 and A.3, 
    according to equation (A.17). Equation (A.17) references equation (A.15) which, 
    although defined in this code later, is simply the piecewise combination of 
    equations (A.1-6, 12-14) described above.

    Args:
        parameters (list): parameters of bifurcation unit, see docstring
    """
    phi_gamma = get_angle_phi_gamma(parameters)
    position_gamma = get_daughter_pipe_positive_z(phi_gamma, -1*np.pi/2, parameters) #(A.17)
    return position_gamma

def get_alpha_gamma_xz_projection(parameters):
    """
    A method for getting the constants a and b from the projection of the 
    ALPHA-GAMMA line onto the xz plane in equation (A.18).

    Args:
        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Get the positional values of ALPHA and GAMMA
    [x_alpha, _, z_alpha] = get_position_alpha(parameters)
    [x_gamma, _, z_gamma] = get_position_gamma(parameters)
    #calculate the constants a and b per equation (A.18)
    a = (z_gamma - z_alpha)/(x_gamma - x_alpha)
    b = z_alpha - ((z_gamma - z_alpha)/(x_gamma - x_alpha)) * x_alpha
    return a, b

####################################################################################################
#                                                                                                  #
#             Functions describing the maximum and minimum circumferential angles                  #
#                in the positive- and negative-z domain at a give angle phi_s                      #
#                                                                                                  #
####################################################################################################

def get_phi_c_pos_min(phi_s, parameters):
    """
    Calculate the minimum circumferential angle in the positive-z domain (phi_c_pos_min) at
    angle phi_s along the streamline parametric axis according to equation (A.19). Valid for 
    any value phi_s in the domain. 

    Args:
        phi_s (float): streamline parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    parent_length = parameters[1]
    daughter_radius = parameters[2]
    outer_radius = parameters[4]
    branching_angle = parameters[5]
    # Get constants a and b, get phi_GAMMA
    a, b = get_alpha_gamma_xz_projection(parameters)
    phi_gamma = get_angle_phi_gamma(parameters)
    #define piecewise function (A.19)
    if phi_s < 0:
        phi_c_pos_min = -np.pi/2
    elif phi_s >= 0 and phi_s <= branching_angle:
        numerator = b + (a * (parent_length + outer_radius * np.sin(phi_s))) - (outer_radius * (1-np.cos(phi_s)))
        denominator = (np.cos(phi_s) + a * np.sin(phi_s)) * transition_radius(phi_s, parameters)
        phi_c_pos_min = np.arcsin(numerator/denominator)
    elif phi_s > branching_angle and phi_s <= phi_gamma:
        numerator_1 = b + (a * (parent_length + outer_radius * np.sin(branching_angle))) - (outer_radius * (1-np.cos(branching_angle)))
        numerator_2 = (a * np.cos(branching_angle) - np.sin(branching_angle)) * outer_radius * np.tan(phi_s - branching_angle)
        denominator = (np.cos(branching_angle) + a * np.sin(branching_angle)) * daughter_radius
        phi_c_pos_min = np.arcsin((numerator_1 + numerator_2)/(denominator))
    elif phi_s > phi_gamma:
        phi_c_pos_min = -np.pi/2
    else:
        #Sanity check/Error log in case of debugging
        print("Error: Could not evaluate phi_c_pos_min given values. Printing Error log... \n")
        print(f"phi_s value: {phi_s}")
        print(f"PHI_GAMMA value: {phi_gamma}")
        print(f"a, b value: {a}, {b}")
        phi_c_pos_min = None
        assert False
    return phi_c_pos_min

def get_phi_c_pos_max(phi_c_pos_min):
    """
    Calculate phi_c_pos_max from phi_c_pos_min according to Equation (A.7).

    Args:
        phi_s (float): streamline parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    phi_c_pos_max = np.pi - phi_c_pos_min #(A.7)
    return phi_c_pos_max

def get_phi_c_neg_min(phi_c_pos_min):
    """
    Calculate phi_c_neg_min from phi_c_pos_min according to Equation (A.8).

    Args:
        phi_c_cos_min (float) 
    """
    phi_c_neg_min = phi_c_pos_min + np.pi #(A.8)
    return phi_c_neg_min

def get_phi_c_neg_max(phi_c_pos_min):
    """
    Calculate phi_c_pos_max from phi_c_pos_min according to Equation (A.9).

    Args:
        phi_s (float): streamline parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    phi_c_pos_max = np.pi - phi_c_pos_min
    phi_c_neg_max = phi_c_pos_max + np.pi #(A.9)
    return phi_c_neg_max

####################################################################################################
#                                                                                                  #
#                Functions describing the maximum and minimum streamline angles                    #
#                            given the bifurcation unit parameters.                                #
#                                                                                                  #
####################################################################################################

def get_phi_s_min(parameters):
    """
    Method for calculating the minimum value of the streamline angle parameteric 
    variable phi_s based on the parameters. Method is adapted from the ranges of transformed 
    cartesian parameteric variable xi introduced in equation (A.1). The derivation of these
    ranges was done analytically by hand (by me).

    Args:
        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    parent_length = parameters[1]
    outer_radius = parameters[4]
    #calculate phi_s_min
    phi_s_min = np.arctan(-1*parent_length/outer_radius)
    assert phi_s_min <= 0 #sanity check that phi_s_min is within the expected range
    return phi_s_min

def get_phi_s_max(parameters):
    """
    Method for calculating the maximum value of the streamline angle parameteric 
    variable phi_s based on the parameters. Method is adapted from the ranges of transformed 
    cartesian parameteric variable zeta introduced in equation (A.12). The derivation of these
    ranges was done analytically by hand (by me).

    Args:
        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    daughter_length = parameters[3]
    outer_radius = parameters[4]
    branching_angle = parameters[5]
    #calculate phi_s_max
    phi_gamma = get_angle_phi_gamma(parameters)
    delta = np.tan(phi_gamma - branching_angle)
    phi_s_max = branching_angle + np.arctan((daughter_length/outer_radius) + delta)
    assert phi_s_max >= branching_angle #sanity check that phi_s_max is in expected range
    return phi_s_max

####################################################################################################
#                                                                                                  #
#             Functions calculating the unit normal vectors of the parent, merging,                #
#                     and daughter pipe sections of the bifurcation unit.                          #
#                                                                                                  #
####################################################################################################

def get_parent_pipe_normal(phi_s, phi_c, parameters):
    """
    Construct the unit normal vector of the parent pipe section at parameteric
    location (phi_s, phi_c) for given parameters.

    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Get the xyz value
    xyz_p = get_parent_pipe(phi_s, phi_c, parameters)
    #calculate unit normal vectors to the pointcloud.
    scalar = 1.1
    parent_dummy = np.array([xyz_p[0], scalar*xyz_p[1], scalar*xyz_p[2]])
    #Create vector and normalize
    normals_p = (parent_dummy - xyz_p) * (1/np.linalg.norm((parent_dummy - xyz_p)))
    return normals_p

def get_merging_pipe_positive_z_normal(phi_s, phi_c, parameters):
    """
    Construct the normal vector of the merging pipe section in the positive-z 
    range at parameteric location (phi_s, phi_c) for given parameters.

    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    parent_length = parameters[1]
    outer_radius = parameters[4]
    #Calculate xyz coordinates
    xyz_m_pos = get_merging_pipe_positive_z(phi_s, phi_c, parameters)
    #Calculate unit vectors normal to the pointcloud
    val = transition_radius(phi_s, parameters) * np.sin(phi_c)
    scalar = 1.1
    dummy_x_m_pos = parent_length + (outer_radius * np.sin(phi_s)) - (np.sin(phi_s) * val * scalar)
    dummy_y_m_pos = transition_radius(phi_s, parameters) * np.cos(phi_c) * scalar
    dummy_z_m_pos = outer_radius * (1-np.cos(phi_s)) + (np.cos(phi_s) * val * scalar)
    merge_dummy_pos = np.array([dummy_x_m_pos, dummy_y_m_pos, dummy_z_m_pos])
    #create vector and normalize
    normals_m_pos = (merge_dummy_pos - xyz_m_pos) * (1/np.linalg.norm(merge_dummy_pos - xyz_m_pos))
    return normals_m_pos

def get_merging_pipe_negative_z_normal(phi_s, phi_c, parameters):
    """
    Construct the normal vector of the merging pipe section in the negative-z 
    range at parameteric location (phi_s, phi_c) for given parameters.

    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    parent_length = parameters[1]
    outer_radius = parameters[4]
    #Calculate xyz coordinates and pack into numpy array
    xyz_m_neg = get_merging_pipe_negative_z(phi_s, phi_c, parameters)
    #Calculate unit vectors normal to the pointcloud
    val = transition_radius(phi_s, parameters) * np.sin(phi_c)
    scalar = 1.1
    dummy_x_m_neg = parent_length + (outer_radius * np.sin(phi_s)) + (np.sin(phi_s) * val * scalar)
    dummy_y_m_neg = transition_radius(phi_s, parameters) * np.cos(phi_c) * scalar
    dummy_z_m_neg = -1 * outer_radius * (1-np.cos(phi_s)) + (np.cos(phi_s) * val * scalar)
    merge_dummy_neg = np.array([dummy_x_m_neg, dummy_y_m_neg, dummy_z_m_neg])
    #create vector and normalize
    normals_m_neg = (merge_dummy_neg - xyz_m_neg) * (1/np.linalg.norm(merge_dummy_neg - xyz_m_neg))
    return normals_m_neg

def get_daughter_pipe_positive_z_normal(phi_s, phi_c, parameters):
    """
    Construct the normal vector of the daughter pipe section in the positive-z 
    range at parameteric location (phi_s, phi_c) for given parameters.

    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    outer_radius = parameters[4]
    branching_angle = parameters[5]
    #Calculate xyz coordinates
    xyz_d_pos = get_daughter_pipe_positive_z(phi_s, phi_c, parameters)
    #Calculate unit vectors normal to the pointcloud
    scalar = 1.1
    zeta = outer_radius * np.tan(phi_s - branching_angle) #cartesian transform of phi_s
    dummy_parameters = parameters
    dummy_parameters[0] = parameters[0] * scalar
    dummy_parameters[2] = parameters[2] * scalar
    xyz_m_dummy_pos, _ = get_merging_pipe_positive_z(branching_angle, phi_c, dummy_parameters)
    dummy_x_d_pos = xyz_m_dummy_pos[0] + zeta*np.cos(branching_angle)
    dummy_y_d_pos = xyz_m_dummy_pos[1]
    dummy_z_d_pos = xyz_m_dummy_pos[2] + zeta*np.sin(branching_angle)
    daughter_dummy_pos = np.array([dummy_x_d_pos, dummy_y_d_pos, dummy_z_d_pos])
    #create vector and normalize
    normals_d_pos = (daughter_dummy_pos - xyz_d_pos) * (1/np.linalg.norm(daughter_dummy_pos - xyz_d_pos))
    return normals_d_pos

def get_daughter_pipe_negative_z_normal(phi_s, phi_c, parameters):
    """
    Construct the normal vector of the daughter pipe section in the negative-z 
    range at parameteric location (phi_s, phi_c) for given parameters.

    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    outer_radius = parameters[4]
    branching_angle = parameters[5]
    #Calculate xyz coordinates
    xyz_d_neg = get_daughter_pipe_negative_z(phi_s, phi_c, parameters)
    #Calculate unit vectors normal to the pointcloud
    zeta = outer_radius * np.tan(phi_s - branching_angle) #cartesian transform of phi_s
    scalar = 1.1
    dummy_parameters = parameters
    dummy_parameters[0] = parameters[0] * scalar
    dummy_parameters[2] = parameters[2] * scalar
    xyz_m_dummy_neg, _ = get_merging_pipe_negative_z(branching_angle, phi_c, dummy_parameters)
    dummy_x_d_neg = xyz_m_dummy_neg[0] + zeta*np.cos(branching_angle)
    dummy_y_d_neg = xyz_m_dummy_neg[1]
    dummy_z_d_neg = xyz_m_dummy_neg[2] - zeta*np.sin(branching_angle)
    daughter_dummy_neg = np.array([dummy_x_d_neg, dummy_y_d_neg, dummy_z_d_neg])
    #create vector and normalize
    normals_d_neg = (daughter_dummy_neg - xyz_d_neg) * (1/np.linalg.norm(daughter_dummy_neg - xyz_d_neg))
    return normals_d_neg

####################################################################################################
#                                                                                                  #
#            Functions for defining the positive and negative position and normal                  #
#               vectors of the bifurcation unit from the individual dections                       #
#                                                                                                  #
####################################################################################################

def get_r_parent(phi_s, phi_c, parameters):
    """
    Method combining equations (A.1-6, 12-14) together piecewise into equation (A.15) to build the
    pointcloud. This method only accepts values of phi_c and phi_s that are in the necessary range.
    This method builds the parent section of the bifurcation unit.
    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Calculate minimum and maximum angle values
    phi_s_min = get_phi_s_min(parameters)
    phi_c_pos_min = get_phi_c_pos_min(phi_s, parameters)
    phi_c_pos_max = get_phi_c_pos_max(phi_c_pos_min)
    #Check here that phi_c is in range, instead of in every condition
    assert (phi_c >= phi_c_pos_min and phi_c <= phi_c_pos_max)
    #construct the piecewise equation (A.15), ensurung phi_s is in range
    if phi_s >= phi_s_min and phi_s <= 0:
        r_parent = get_parent_pipe(phi_s, phi_c, parameters)
    else:
        print("Error: Could not evaluate r_positive given phi_s. Printing Error log... \n")
        print(f"phi_s value: {phi_s}")
        r_parent = None
        assert False

    return r_parent

def get_r_positive(phi_s, phi_c, parameters):
    """
    Method combining equations (A.1-6, 12-14) together piecewise into equation (A.15) to build the
    pointcloud. This method only accepts values of phi_c and phi_s that are in the necessary range.
    This method builds the positive-z domain of the bifurcation unit.
    
    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    branching_angle = parameters[5]
    #Calculate minimum and maximum angle values
    phi_s_max = get_phi_s_max(parameters)
    phi_c_pos_min = get_phi_c_pos_min(phi_s, parameters)
    phi_c_pos_max = get_phi_c_pos_max(phi_c_pos_min)
    #Check here that phi_c is in range, instead of in every condition
    assert (phi_c >= phi_c_pos_min and phi_c <= phi_c_pos_max)
    #construct the piecewise equation (A.15), ensurung phi_s is in range
    if phi_s >= 0 and phi_s <= branching_angle:
        r_positive = get_merging_pipe_positive_z(phi_s, phi_c, parameters)
    elif phi_s > branching_angle and phi_s <= phi_s_max:
        r_positive = get_daughter_pipe_positive_z(phi_s, phi_c, parameters)
    else:
        print("Error: Could not evaluate r_positive given phi_s. Printing Error log... \n")
        print(f"phi_s value: {phi_s}")
        r_positive = None
        assert False

    return r_positive

def get_r_negative(phi_s, phi_c, parameters):
    """
    Method combining equations (A.1-6, 12-14) together piecewise into equation (A.15) to build the
    pointcloud. This method only accepts values of phi_c and phi_s that are in the necessary range.
    This method builds the negative-z domain of the bifurcation unit.
    
    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    branching_angle = parameters[5]
    #Calculate minimum and maximum angle values
    phi_s_max = get_phi_s_max(parameters)
    phi_c_pos_min = get_phi_c_pos_min(phi_s, parameters)
    phi_c_neg_min = get_phi_c_neg_min(phi_c_pos_min)
    phi_c_neg_max = get_phi_c_neg_max(phi_c_pos_min)
    #Check here that phi_c is in range, instead of in every condition
    assert (phi_c >= phi_c_neg_min and phi_c <= phi_c_neg_max)
    #construct the piecewise equation (A.15), ensurung phi_s is in range
    if phi_s >= 0 and phi_s <= branching_angle:
        r_negative = get_merging_pipe_negative_z(phi_s, phi_c, parameters)
    elif phi_s > branching_angle and phi_s <= phi_s_max:
        r_negative = get_daughter_pipe_negative_z(phi_s, phi_c, parameters)
    else:
        print("Error: Could not evaluate r_positive given phi_s. Printing Error log... \n")
        print(f"phi_s value: {phi_s}")
        r_negative = None
        assert False

    return r_negative

#TODO get parent and branching normals????

####################################################################################################
#                                                                                                  #
#                 Functions to prep for creating the carina geometry pointcloud                    #
#                                                                                                  #
####################################################################################################

def get_original_carina_unit_vectors(phi_s, phi_c, parameters):
    """
    Method for calculating the unit vectors tangential (along the streamline and circumference) and 
    normal to a constant transitional radius at phi_s, using equation (A.20). These vectors form a 
    right hand system
    
    Args:
        phi_s (float): streamline parameteric angle, in radians

        phi_c (float): circumferential parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    branching_angle = parameters[5]
    #Calculate the vectors according to equation (A.20)
    phi_s_tilde = np.max([0, np.min([phi_s, branching_angle])])
    #Unit vector tangential to the circumference
    e_circ = np.array([-1*np.sin(phi_s_tilde) * np.cos(phi_c), #(A.20)
                       -1*np.sin(phi_c),
                       np.cos(phi_s_tilde) * np.cos(phi_c)])
    #unit vector tangential to a constant radiusstreamline
    e_stream_prime = np.array([np.cos(phi_s_tilde), #(A.20)
                               0,
                               np.sin(phi_s_tilde)])
    #unit vector normal to a constant radius streamline
    normal_prime = np.array([-1*np.sin(phi_s_tilde) * np.sin(phi_c), #(A.20)
                             np.cos(phi_c),
                             np.cos(phi_s_tilde) * np.sin(phi_c)])
    
    return e_circ, e_stream_prime, normal_prime

def get_upsilon(phi_s, parameters):
    """
    Method for calculating the angle upsilon at phi_s based on equation (A.21).
    This angle is used to transform the unit vectors above to the correct
    orientation and location to construct the carina.
    
    Args:
        phi_s (float): streamline parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Unpack parameters list for any required parameters in this function
    parent_radius = parameters[0]
    daughter_radius = parameters[2]
    outer_radius = parameters[4]
    branching_angle = parameters[5]
    #Calculate upsilon
    val1 = (branching_angle - phi_s)/branching_angle
    val2 = (parent_radius - daughter_radius)/(outer_radius * branching_angle)
    upsilon = np.arctan(-6 * (phi_s/branching_angle)*val1*val2) #(A.21)
    return upsilon

def transform_carina_unit_vectors(phi_s, phi_c, parameters):
    """
    Method for transforming the carina unit vectors e_streamline_prime and 
    normal_prime from equation (A.20) to the correct orientation and location
    using matrix multiplication by Q defined in Equation (A.22), according to 
    Equations (A.23) and (A.24).
    
    Args:
        phi_s (float): streamline parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Get upsilon and the original unit vectors
    upsilon = get_upsilon(phi_s, parameters)
    e_circ, e_stream_prime, normal_prime = get_original_carina_unit_vectors(phi_s, phi_c, parameters)
    #Define the transformation matrix Q 
    A = np.array([[0,            -1*e_circ[2],  e_circ[1]   ],         #(A.22)
                  [e_circ[2],    0,             -1*e_circ[0]],
                  [-1*e_circ[1], e_circ[0],     0           ]])
    Q = np.identity(3) + np.sin(upsilon)*A + (1-np.cos(upsilon)) * A@A #(A.22)
    #Transform the original streamline tangent and normal vectors
    e_stream = Q @ e_stream_prime   #(A.23)
    unit_normal = Q @ normal_prime  #(A.24)
    return e_stream, unit_normal

def get_carina_position_and_normal_vectors(phi_s, parameters):
    """
    Method for calculating the position vector of location LAMBDA, normal
    vector of location LAMBDA, position vector of location KAPPA, position
    vector of location MU, and position vector of location NU in Figure A.3,
    using Equations (A.25), (A.26), (A.27), (A.28), and (A.29). These 
    values are necessary to construct the carina.
    
    Args:
        phi_s (float): streamline parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Find phi_c_pos_min at the given phi
    phi_c_pos_min = get_phi_c_pos_min(phi_s, parameters)
    #get position LAMBDA
    position_lambda = get_r_positive(phi_s, phi_c_pos_min, parameters) #(A.25)
    #get normal vector at position LAMBDA
    _, normal_lambda = transform_carina_unit_vectors(phi_s, phi_c_pos_min, parameters) #(A.26)
    #get position KAPPA
    position_kappa = position_lambda - ((position_lambda[2]/normal_lambda[2])*normal_lambda) #(A.27)
    #get position MU
    position_mu = position_lambda - ((position_lambda[1]/normal_lambda[1]) * normal_lambda) #(A.28)
    #get position NU
    position_nu = np.array([position_mu[0], 0, 0]) #(A.29)
    return position_lambda, position_kappa, position_mu, position_nu, normal_lambda

def get_local_carina_radius_and_angles(phi_s, parameters):
    """
    Method for calculating the value of the carina radius (rho) and the
    range of the carina angle (omega) as shown in Figure A.3, using Equations
    (A.30) and (A.31). Also calculates the values of sin(theta) and cos(theta),
    angles relating the positions of KAPPA and NU, using equations (A.32) and (A.33), 
    to reduce the number of times the position vectors are calculated.
    These values are necessary to construct the carina.
    
    Args:
        phi_s (float): streamline parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring

    Returns: 
        carina_radius (float): local radius of the carina at phi_s

        carina_angle (float): local range of the carina angle

        sin_theta (float): value necessary for the construction of the carina

        cos_theta (float): value necesarry for the construction of the carina
    """
    #Get necessary position vectors
    position_lambda, position_kappa, position_mu, position_nu, _ = \
                                    get_carina_position_and_normal_vectors(phi_s, parameters)
    #Calculate the carina radius
    carina_radius = np.linalg.norm((position_kappa - position_lambda)) #(A.30)
    #Calculate the carina angle
    carina_angle = np.arctan(position_mu[2]/(np.linalg.norm((position_kappa - position_nu))))#(A.31)
    #Calculate sin(theta) and cos(theta)
    sin_theta = (position_kappa[0] - position_nu[0])/(np.linalg.norm((position_kappa - position_nu))) #(A.32)
    cos_theta = (position_kappa[1])/(np.linalg.norm((position_kappa - position_nu))) #(A.33)
    return carina_radius, carina_angle, sin_theta, cos_theta


####################################################################################################
#                                                                                                  #
#            Functions for defining the positive and negative position vectors and                 #
#                      the unit normal vectors of the carina pointcloud                            #
#                                                                                                  #
####################################################################################################

def get_carina_position(phi_s, psi, parameters):
    """
    Construct the xyz position of the bifurcation carina in the positive-y 
    range (denoted up) at parameteric location (phi_s, psi) for given parameters. 
    Equations (A.34), (A.35), (A.36)

    Args:
        phi_s (float): streamline parameteric angle, in radians

        psi (float): carina parametric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #Get necessary position vectors, carina radius, and theta angles
    _, position_kappa, _, _, _ = get_carina_position_and_normal_vectors(phi_s, parameters)
    carina_radius, carina_angle, sin_theta, cos_theta = get_local_carina_radius_and_angles(phi_s, parameters)
    #Check that phi_s and psi are within range
    phi_gamma = get_angle_phi_gamma(parameters)
    assert (phi_s >= 0 and phi_s <= phi_gamma)
    assert (psi >= -1*carina_angle and psi <= carina_angle)
    #calculate the xyz position of the carina point
    x_carina_up = position_kappa[0] - (carina_radius * np.cos(psi) * sin_theta)
    y_carina_up = position_kappa[1] - (carina_radius * np.cos(psi) * cos_theta)
    z_carina_up = carina_radius * np.sin(psi)
    r_carina_up = np.array([x_carina_up, y_carina_up, z_carina_up])
    r_carina_down = np.array([x_carina_up, -1*y_carina_up, z_carina_up])
    return r_carina_up, r_carina_down


def get_carina_normal(phi_s, psi, parameters):
    """
    Construct the unit normal vector of the positive-y carina at parameteric
    location (phi_s, psi) for given parameters.

    Args:
        phi_s (float): streamline parameteric angle, in radians

        psi (float): carina parameteric angle, in radians

        parameters (list): parameters of bifurcation unit, see docstring
    """
    #get position KAPPA and normal LAMBDA vectors
    _, position_kappa, _, _, normal_lambda = get_carina_position_and_normal_vectors(phi_s, parameters)
    #get r_carina_up
    r_carina_up, _ = get_carina_position(phi_s, psi, parameters)
    normal_carina_up = -1*(position_kappa - r_carina_up) / np.linalg.norm((position_kappa - r_carina_up))*np.sign(normal_lambda[2])
    normal_carina_down = np.array([normal_carina_up[0], -1*normal_carina_up[1], normal_carina_up[2]])
    return normal_carina_up, normal_carina_down
