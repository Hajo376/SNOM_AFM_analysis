# import numpy as np
# from matplotlib import pyplot as plt



def Horizontal_Profile(array):
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
