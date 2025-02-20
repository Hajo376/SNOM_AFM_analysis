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

# import numpy as np
# from matplotlib import pyplot as plt



def horizontal_profile(array):
    xres = len(array[0])
    yres = len(array)
    print(f'xres: {xres}')
    print(f'yres: {yres}')
    profile = []
    for x in range(xres):
        mean = 0
        for y in range(yres):
            mean += array[y][x]/yres
        profile.append(mean)
    return profile
