"""
Functions for the pre-processing of the general bifurcation tree. Includes reading in the .csv file,
ensuring that the number of bifurcation units per generation matches the number of available/open
continuation outlets from the previous generation, and that all parameters are within acceptable 
ranges (TBD). Will also include functions for initializing the position and orientation of each unit
in the tree.
"""
import csv
import numpy as np

def read_csv_parameters_to_numpy_array(tree_parameterscsv_filename: str):
    """
    Read the tree parameters .csv file (see Github README for information) into the function to 
    be converted into a numpy array for use in the bifurcation tree generation. 

    Args:
        tree_parameters_csv_filename (str): the filename (and path, if necessary) of the tree
        parameter .csv file

    Returns:
        tree_parameters: numpy array of the tree parameters
    """
    # Read .csv file into python as list
    with open(tree_parameterscsv_filename) as fp:
        reader = csv.reader(fp, delimiter=",", quotechar='"')
        next(reader, None)  # skip the headers
        data_read = [row for row in reader]
    #convert data to numpy array
    tree_parameters = np.array(data_read).astype(np.float64)
    return tree_parameters

#From here, there are a number of small functions to get statistics of the bifurcation tree
#and to check their validity.

def check_number_of_columns(tree_parameters):
    """
    Check that the number of columns in the tree parameters numpy array is the 
    expected 9 
    """
    n_params = np.shape(tree_parameters)[1]
    if n_params != 9:
        raise AssertionError(f"Incorrect format of tree parameters. Has {n_params}" + \
                              "columns, should have 9 columns") 
    else:
        return

def check_initial_generation_index(tree_parameters):
    """
    Check that the initial generation index of the of the tree parameters is 0. If not,
    corrects the indexing error and returns the revised tree parameters list.
    """
    initial_generation_index = tree_parameters[0,0]
    if int(initial_generation_index) != 0:
        print(f"Incorrect Initial Generation Index! Is: {initial_generation_index}. " + \
              "Should be 0. \nCorrecting tree Parameters List...")
        tree_parameters[:,0] = int(tree_parameters[:,0] - initial_generation_index)
        return tree_parameters
    else:
        return tree_parameters


