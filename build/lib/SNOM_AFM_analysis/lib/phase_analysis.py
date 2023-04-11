import numpy as np

def Flatten_Phase_Profile(profile:list, phase_orientation:int=1) -> list:
    flattened_profile = []
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

    return flattened_profile

def Get_Smallest_Difference(value1, value2):
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

def Get_Difference(value1, value2):
    difference = value1-value2
    if difference < 0:
        difference += 2*np.pi
    return difference

def Get_Difference_2(value1, value2):
    return value1-value2
    
def Get_Profile_Difference(profile1:list, profile2:list) -> list:
    # difference = [abs(profile1[i] - profile2[i]) for i in range(len(profile1))]
    # difference = [profile1[i] - profile2[i] for i in range(len(profile1))]
    # difference = [abs(profile1[i] - profile2[i]) if abs(profile1[i] - profile2[i])< np.pi else 2*np.pi - abs(profile1[i] - profile2[i]) for i in range(len(profile1))]
    difference = []
    for i in range(len(profile1)):
        # difference.append(Get_Smallest_Difference(profile1[i], profile2[i]))
        difference.append(Get_Difference(profile1[i], profile2[i]))

    return difference
    # pass

def Get_Profile_Difference_2(profile1:list, profile2:list) -> list:
    difference = []
    for i in range(len(profile1)):
        difference.append(Get_Difference_2(profile1[i], profile2[i]))
    return difference

def Get_Modeindex_from_LinearFunction(slope, pixelsize, wavelength=1600):
    # wavelength and pixelsize in nm
    period = np.pi*2/slope*pixelsize
    mode_index = wavelength/period
    return mode_index
