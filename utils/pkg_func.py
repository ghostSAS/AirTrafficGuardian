import numpy as np
import pandas as pd

from scipy.special import gamma, factorial, comb

import matplotlib.pyplot as plt
from matplotlib import animation
import mpl_toolkits.mplot3d.axes3d as p3

import time


############################ Helper functions ###################################

def timer(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} takes {end_time - start_time:.4f} seconds")
        return result
    return wrapper

def combinations_of_2(l):
    for i, j in zip(*np.triu_indices(len(l), 1)):
        yield l[i], l[j]
        

def is_type(input, inputtype):
    """
    Checks if the input is of a given type, raises an error otherwise

    :param input: Any input
    :type  input: Any
    :param inputtype: Type class
    :type  inputtype: TypeA

    :raises TypeError: When inputs type does not match the inputtype argument

    :return: The function input
    :rtype: TypeA
    """
    if isinstance(input, inputtype):
        return input
    else:
        raise TypeError("Input must be of type "+str(inputtype))

def is_positive(input):
    """
    Checks if the a number is greater than 0

    :param input: A number
    :type  input: <numbers.Real

    :raises ValueError: When number is less than or equal to 0

    :return: input
    :rtype: TypeA
    """
    if input > 0:
        return input
    else:
        raise ValueError("Input must be positive")

def is_valid_spline_type(input):
    """
    Checks that the input is a string of defining an appropriate
    spline type

    :param input: A string, either "polynomial" or "bezier" 
    :type  input: string

    :raises TypeError: If the input is not a string
    :raises ValueError: If the input is not "polynomial" or "bezier" 

    :return: input
    :rtype: string
    """
    input = is_type(input, str).lower()
    valid_splines = ["poly", "bezier"]
    if input in valid_splines:
        return input
    else:
        raise ValueError("Input must be either polynomial or bezier")

def is_valid_projection_type(input):
    """
    Checks that the input is a string of defining an appropriate
    projection type

    :param input: A string, either "1d" or "3d" 
    :type  input: string

    :raises TypeError: If the input is not a string
    :raises ValueError: If the input is not "1d" or "3d" 

    :return: input
    :rtype: string
    """
    input = is_type(input, str).lower()
    valid_splines = ["3d", "1d"]
    if input in valid_splines:
        return input
    else:
        raise ValueError("Input must be either 1D or 3D")


def is_valid_solver_type(input):
    """
    Checks that the input is a string of defining an appropriate
    solver type

    :param input: A string, either "kkt" or "cpo" 
    :type  input: string

    :raises TypeError: If the input is not a string
    :raises ValueError: If the input is not "polynomial" or "bezier" 

    :return: input
    :rtype: string
    """
    input = is_type(input, str).lower()
    valid_splines = ["kkt", "qp"]
    if input in valid_splines:
        return input
    else:
        raise ValueError("Input must be either kkt or qp")
    
    
def convert2_np_array(mat):
    assert(isinstance(mat, str))
    
    mat = mat.replace('[','').replace(']','')
    mat = " ".join(mat.split())
    
    rows = mat.split(';')
    row_dig = []
    for row in rows:
        row = row.strip()

        nums = row.split(' ')        
        num_dig = [float(num) for num in nums]
        row_dig.append(num_dig)
    
    return np.array(row_dig)
        

