import numpy as np
from scipy.special import gamma, factorial, comb

import matplotlib.pyplot as plt
from matplotlib import animation
import mpl_toolkits.mplot3d.axes3d as p3

import time

def timer(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} takes {end_time - start_time:.4f} seconds")
        return result
    return wrapper
