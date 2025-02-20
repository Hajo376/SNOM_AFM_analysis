"""Some additional functions for the snom_analysis package."""

##############################################################################
# Copyright (C) 2020-2025 Hans-Joachim Schill

# This file is part of snom_analysis.

# snom_analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# snom_analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with snom_analysis.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

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

def calculate_colorbar_size(fig, ax, colorbar_size):
    """This function converts a colorbar size in % of the fig width to a colorbar size in % of the axis width."""
    # size of the figure in inches
    fig_width = fig.get_figwidth()
    # size of the axis in inches
    ax_width = ax.get_position().width * fig_width
    # calculate the size of the colorbar in percent
    # it should always be x % of the figure width
    return colorbar_size * fig_width / ax_width 