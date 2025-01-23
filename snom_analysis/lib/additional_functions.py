"""Some additional functions for the snom_analysis package."""

import numpy as np

def set_nan_to_zero(data) -> np.array:
    xres = len(data[0])
    yres = len(data)
    for y in range(yres):
        for x in range(xres):
            if str(data[y][x]) == 'nan':
                data[y][x] = 0
    return data       

# needed for the realign function
def gauss_function(x, A, mu, sigma, offset):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + offset

def get_largest_abs(val1, val2):
    if abs(val1) > abs(val2): return abs(val1)
    else: return abs(val2)