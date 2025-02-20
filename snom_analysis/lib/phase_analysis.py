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

def flatten_phase_profile(profile:list, phase_orientation:int=1) -> list:
    """flattened_profile = []
    previous_element = profile[0]
    offset = 0
    for element in profile:
        if (phase_orientation == 1) and (element - previous_element < -np.pi):
            offset += 2*np.pi
            # print('increased offset')
        elif (phase_orientation == -1) and (element - previous_element > np.pi):
            offset -= 2*np.pi
            # print('reduced offset')
        flattened_profile.append(element + offset)
        previous_element = element

    return flattened_profile"""
    # new alternative, use built-in numpy function unwrap
    return np.unwrap(profile)


def get_smallest_difference(value1, value2):
    # make shure value1 is smaller than value2
    orientation = 1
    if value1 > value2:
        copy = value1
        value1 = value2
        value2 = copy
        # orientation = -1
    difference = abs(value1 - value2)
    if difference > np.pi:
        difference = np.pi*2 - difference
    return difference*orientation

def get_difference(value1, value2):
    difference = value1-value2
    if difference < 0:
        difference += 2*np.pi
    return difference

def get_difference_2(value1, value2):
    return value1-value2
    
def get_profile_difference(profile1:list, profile2:list) -> list:
    # difference = [abs(profile1[i] - profile2[i]) for i in range(len(profile1))]
    # difference = [profile1[i] - profile2[i] for i in range(len(profile1))]
    # difference = [abs(profile1[i] - profile2[i]) if abs(profile1[i] - profile2[i])< np.pi else 2*np.pi - abs(profile1[i] - profile2[i]) for i in range(len(profile1))]
    difference = []
    for i in range(len(profile1)):
        # difference.append(get_smallest_difference(profile1[i], profile2[i]))
        difference.append(get_difference(profile1[i], profile2[i]))

    return difference
    # pass

def get_profile_difference_2(profile1:list, profile2:list) -> list:
    difference = []
    for i in range(len(profile1)):
        difference.append(get_difference_2(profile1[i], profile2[i]))
    return difference

def get_modeindex_from_linearfunction(slope, pixelsize, wavelength=1600):
    # wavelength and pixelsize in nm
    period = np.pi*2/slope*pixelsize
    mode_index = wavelength/period
    return mode_index
