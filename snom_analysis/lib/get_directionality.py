'''This script is ment to extract the directionality from a measurement of a chiral coupler.'''

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
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import signal
import pathlib
this_files_path = pathlib.Path(__file__).parent.absolute()
import os
import tkinter as tk
from tkinter import filedialog

from .snom_colormaps import *

class ChiralCoupler():
    def __init__(self, height_data, data_array, channels, orientation='vertical'):
        self.height_data = height_data
        self.data_array = data_array
        self.channels = channels
        # self.integ_width = integ_width
        self.orientation = orientation

        self.wg_positions = None
        # if 'A' in channel:
        #     self.display_channel = 'A'
        # elif 'P' in channel:
        #     self.display_channel = 'P'
        # self.display_channel = 
        # data used for displaying with the fit to the waveguides, typically the height profile is shown as this is the basis for the fits
        self.display_data = self.height_data
        self.Fit_Wgs()

    def Fit_Wgs(self) -> None:
        '''This fits the data and returns a list with both positions of the waveguides.
        Like: [row, col_wg1, col_wg2]
        '''
        if self.orientation == 'vertical':
            width = len(self.height_data[0])
            p0 = [100, (width)/4, (width)/4*3, 5, 0]
            bounds = ([0, 0, width/2, 0, -1000], [1000, width/2, width, width/2, 1000])
            self.list_of_coefficients = []
            print('starting the fitting procedure')
            for row in range(len(self.height_data)):
                coeff, var_matrix = curve_fit(Double_Gauss, range(0, width), self.height_data[row], p0=p0, bounds=bounds)
                self.list_of_coefficients.append(coeff)
                p0 = coeff #set the starting parameters for the next fit
                # print(coeff)
            mean_sigma = [self.list_of_coefficients[i][3] for i in range(len(self.list_of_coefficients))]
            mean_sigma = np.mean(mean_sigma)
            self.integ_width = int(mean_sigma*2)
            self.Get_Integration_Width()

    def Get_Integration_Width(self) -> None:
        self.Plot_Data_With_Fit_And_Bounds(self.height_data, 'Z')
        user_input = input(f'Do you like the extraction width? Current width is {self.integ_width} y/n: ')
        if user_input == 'n':
            self.integ_width = int(input('Please enter the new extraction width as an integer: '))
            self.Get_Integration_Width() # call recursively until user is happy with integ_width
        elif user_input == 'y':
            pass
        else:
            print('Please enter n or y!')

    def Plot_Data_With_Fit_And_Bounds(self, data, channel) -> None:
        fig, axis = plt.subplots()    
        # fig.set_figheight(self.figsizey)
        # fig.set_figwidth(self.figsizex) 
        
        yres = len(data)

        if 'Z' in channel:
            cmap = 'gray'
            label = 'height [nm]'
            axis.plot([self.list_of_coefficients[i][1] for i in range(yres)], range(yres), color='red')
            axis.plot([self.list_of_coefficients[i][2] for i in range(yres)], range(yres), color='red')
        elif 'A' in channel:
            cmap = SNOM_amplitude
            label = 'amplitude [a.u.]'
        elif 'P' in channel:
            cmap = SNOM_phase
            label = 'phase'
        img = axis.pcolormesh(data, cmap=cmap)
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(label, rotation=270)
        axis.set_title('finding the wgs')
        axis.axis('scaled')
        # left wg extraction bounds
        axis.plot([self.list_of_coefficients[i][1]-self.integ_width-1 for i in range(yres)], range(yres), color='white')
        axis.plot([self.list_of_coefficients[i][1]+self.integ_width+1 for i in range(yres)], range(yres), color='white')
        # right wg extraction bounds
        axis.plot([self.list_of_coefficients[i][2]-self.integ_width-1 for i in range(yres)], range(yres), color='white')
        axis.plot([self.list_of_coefficients[i][2]+self.integ_width+1 for i in range(yres)], range(yres), color='white')
        axis.invert_yaxis()
        plt.show()

    def Extract_WG_Data(self, channel):
        wg_positions = [[i, self.list_of_coefficients[i][1], self.list_of_coefficients[i][2]] for i in range(len(self.height_data))]
        # print('wg_positions', wg_positions)
        # print('integ_width: ', integ_width)
        self.left_wg_data = []
        self.right_wg_data = []
        for row in range(len(wg_positions)):
            wg_left = round(wg_positions[row][1])
            wg_right = round(wg_positions[row][2])
            self.left_wg_data.append(self.data_array[self.channels.index(channel)][row][int(wg_left-self.integ_width):int(wg_left+self.integ_width)])
            self.right_wg_data.append(self.data_array[self.channels.index(channel)][row][int(wg_right-self.integ_width):int(wg_right+self.integ_width)])
    
    def Data_Mean(self) -> list:
        mean_data = []
        for channel in self.channels:
            self.Extract_WG_Data(channel)
            left_mean_array = []
            right_mean_array = []
            for row in range(len(self.data_array[self.channels.index(channel)])):
                sum_left = 0
                sum_right = 0
                amp_left = self.left_wg_data[row]
                amp_right = self.right_wg_data[row]
                for j in range(2*self.integ_width):
                    sum_left+=amp_left[j]
                    sum_right+=amp_right[j]
                left_mean_array.append(sum_left/2*self.integ_width)
                right_mean_array.append(sum_right/2*self.integ_width)
            mean_data.append([channel, left_mean_array, right_mean_array])
        return mean_data
    
    def Export_Mean_Data(self, filename):
        mean_data = self.Data_Mean()
        for i in range(len(mean_data)):
            channel = mean_data[i][0]
            self.Export_Data(channel, filename, mean_data[i][1], mean_data[i][2])


    def Integrate_Amplitude(self, channel):
        if channel in self.channels:
            self.Plot_Data_With_Fit_And_Bounds(self.data_array[self.channels.index(channel)], channel)
        else:
            print(f'This channel ({channel}) is not in memory!')
            exit()
        
        self.Extract_WG_Data(channel)
        self.amp_left_integ_array = []
        self.amp_right_integ_array = []
        for row in range(len(self.data_array[self.channels.index(channel)])):
            sum_left = 0
            sum_right = 0
            amp_left = self.left_wg_data[row]
            amp_right = self.right_wg_data[row]
            for j in range(2*self.integ_width):
                sum_left+=amp_left[j]
                sum_right+=amp_right[j]
            # print('sum_left:', sum_left)
            self.amp_left_integ_array.append(sum_left)
            self.amp_right_integ_array.append(sum_right)
        # print('len(data[0])', len(data[0]))
        # print('len(self.amp_left_integ_array)', len(self.amp_left_integ_array))
        plt.plot(range(len(self.amp_left_integ_array)), self.amp_left_integ_array, label='left wg')
        plt.plot(range(len(self.amp_right_integ_array)), self.amp_right_integ_array, label='right wg')
        plt.title(f'{channel}')
        plt.legend()
        plt.show()

    def Display_Phase_Profile(self, channel):
        if channel in self.channels:
            self.Plot_Data_With_Fit_And_Bounds(self.data_array[self.channels.index(channel)], channel)
        else:
            print(f'This channel ({channel}) is not in memory!')
            exit()
        self.Extract_WG_Data(channel)
        self.phase_left_integ_array = []
        self.phase_right_integ_array = []
        for row in range(len(self.data_array[self.channels.index(channel)])):
            sum_left = 0
            sum_right = 0
            phase_left = self.left_wg_data[row]
            phase_right = self.right_wg_data[row]
            N = 2*self.integ_width
            for j in range(N):
                sum_left+=phase_left[j]
                sum_right+=phase_right[j]
            # print('sum_left:', sum_left)
            self.phase_left_integ_array.append(sum_left/(N))
            self.phase_right_integ_array.append(sum_right/(N))
        # print('len(data[0])', len(data[0]))
        # print('len(phase_left_integ_array)', len(phase_left_integ_array))
        # plt.plot(range(len(phase_left_integ_array)), phase_left_integ_array, label='left wg')
        plt.plot(range(len(self.phase_right_integ_array)), self.phase_right_integ_array, label='right wg')
        '''
        coeff = self.Fit_Sawtooth(phase_right_integ_array)
        # x = range(len(phase_left_integ_array))
        x = np.linspace(0, 1, 500)
        plt.plot(x*len(phase_left_integ_array), Sawtooth_v2(x, coeff[0], coeff[1]), label='Fit to right wg')
        '''
        flattened_profile = Flatten_Profile(self.phase_right_integ_array)
        coeff = Fit_Linear_Function(flattened_profile)
        plt.plot(range(len(self.phase_right_integ_array)), flattened_profile, label='flattened phase')
        plt.plot(range(len(self.phase_right_integ_array)), Linear_Function(np.arange(len(self.phase_right_integ_array)), coeff[0], coeff[1]), label='fit: flattened phase')
        plt.title(f'{channel}')
        plt.legend()
        plt.show()

    # def __Export_File(self, filepath, filename, data):
    #     data = open(filename + filename, 'w')


    def Export_Data(self, channel, filename, left_data=None, right_data=None):
        filepath = os.path.normpath(os.path.join(this_files_path, '../chiralcouplers_analysis/extracted_data'))
        if 'A' in channel:
            if left_data==None:
                left_data = self.amp_left_integ_array
            if right_data==None:
                right_data = self.amp_right_integ_array
            filename_extension = f'{channel}_mean.txt'
            with open(filepath + '/' + filename + filename_extension, 'w') as datafile:
                datafile.write('#pixel\tamp leftwg\tamp rightwg\n')
                for i in range(len(self.amp_left_integ_array)):
                    datafile.write(f'{i}\t{round(self.amp_left_integ_array[i], 3)}\t{round(self.amp_right_integ_array[i], 3)}\n')
            # with open(filepath + '/' + f'right_wg_integ_{channel}.txt', 'w') as datafile:
            #     datafile.write('#pixel\tamp')
            #     for i in range(len(self.amp_right_integ_array)):
            #         datafile.write(f'{i}\t{self.amp_right_integ_array[i]}\n')
        elif 'P' in channel:
            if left_data==None:
                left_data = self.phase_left_integ_array
            if right_data==None:
                right_data = self.phase_right_integ_array
            filename_extension = f'{channel}_mean.txt'
            with open(filepath + '/' + filename + filename_extension, 'w') as datafile:
                datafile.write('#pixel\tphase leftwg\tphase rightwg')
                for i in range(len(self.phase_left_integ_array)):
                    datafile.write(f'{i}\t{round(self.phase_left_integ_array[i], 3)}\t{round(self.phase_right_integ_array[i], 3)}\n')
            # with open(filepath + '/' + f'right_wg_integ_{channel}.txt', 'w') as datafile:
            #     datafile.write('#pixel\tamp')
            #     for i in range(len(self.phase_right_integ_array)):
            #         datafile.write(f'{i}\t{self.amp_right_integ_array[i]}\n')

    def Display_Phase_Difference(self, channel):
        if channel in self.channels:
            self.Plot_Data_With_Fit_And_Bounds(self.data_array[self.channels.index(channel)], channel)
        else:
            print(f'This channel ({channel}) is not in memory!')
            exit()
        self.Extract_WG_Data(channel)
        phase_left_integ_array = []
        phase_right_integ_array = []
        phase_difference = []
        for row in range(len(self.left_wg_data)):
            sum_left = 0
            sum_right = 0
            phase_left = self.left_wg_data[row]
            phase_right = self.right_wg_data[row]
            N = 2*self.integ_width
            for j in range(N):
                sum_left+=phase_left[j]
                sum_right+=phase_right[j]
            # print('sum_left:', sum_left)
            phase_left_integ_array.append(sum_left/(N))
            phase_right_integ_array.append(sum_right/(N))
            # phase_difference.append(Get_Smallest_Difference(sum_left/(N), sum_right/(N)))
            phase_difference.append(Get_Phase_Difference_V2(sum_left/(N), sum_right/(N)))

        # extend the phase range by adding copy shifted down by 2 pi
        extended_phase_difference = []
        for phase in phase_difference:
            extended_phase_difference.append(phase)
            if phase >= np.pi:
                extended_phase_difference.append(phase- 2*np.pi)
        plt.plot(range(len(phase_difference)), phase_difference, 'x', label='phase difference')
        plt.title(f'{channel}')
        plt.legend()
        plt.show()
        
        # plt.plot(range(len(extended_phase_difference)), extended_phase_difference, 'x', label='extended phase difference')
        # plt.title(f'{channel}')
        # plt.legend()
        # plt.show()

    def Fit_Sawtooth(self, data):
        x = np.linspace(0, 1, len(data))
        n_periods=8
        period = len(data)/n_periods
        p0 = [n_periods, 0]
        bounds = ([1, -100], [100, 100])
        coeff, var_matrix = curve_fit(Sawtooth_v2, x, data, p0=p0, bounds=bounds)
        print('Fit of Sawtooth successful!')
        print(coeff)
        return coeff
        # x = range(len(data))
        # p0 = [2*np.pi/len(data), 8, 0]
        # bounds = ([1*np.pi/len(data), 1, -2*np.pi/len(data)], [4*np.pi/len(data), 100, 2*np.pi/len(data)])
        # coeff, var_matrix = curve_fit(Sawtooth, x, data, p0=p0, bounds=bounds)
        # print('Fit of Sawtooth successful!')
        # print(coeff)
        # return coeff
        # plt.plot(x, Sawtooth(x, coeff[0], coeff[1], coeff[2]))
        # #plot a sawtooth
        # period = 2 * np.pi /len(phase_left_integ_array)
        # # period = len(phase_left_integ_array)/8
        # n_periods = 8
        # offset = 17
        # phase_orientation = 1
        # t = np.linspace(0, len(phase_left_integ_array), 1000)
        # # t = np.arange(len(phase_left_integ_array))
        # plt.plot(t, np.pi*(1+signal.sawtooth(period * n_periods * (t-offset), phase_orientation)))

def Flatten_Profile(data):
    threshold = np.pi*0.8
    flattened_profile = []
    phase_shift = 0
    for i in range(len(data)):
        if data[i] < (data[i-1]-threshold):
            phase_shift += np.pi*2
        flattened_profile.append(data[i]+phase_shift)
    return flattened_profile

def Fit_Decaying_Amplitude():
    with_oscillation = False
    scaling = 4
    pixelsize = 50/scaling # in nm
    root = tk.Tk()
    root.withdraw()
    filepath = filedialog.askopenfilename(initialdir=os.path.normpath(os.path.join(this_files_path, '../chiralcouplers_analysis/extracted_data')))
    with open(filepath, 'r') as f:
        x = np.genfromtxt(f, delimiter='\t', skip_header=1, usecols=0)
    with open(filepath, 'r') as f:
        leftwg = np.genfromtxt(f, delimiter='\t', skip_header=1, usecols=1)
    with open(filepath, 'r') as f:
        rightwg = np.genfromtxt(f, delimiter='\t', skip_header=1, usecols=2)
        # data = np.loadtxt(f, delimiter='\t', skip_header=1)
    # print('x: ', x)
    # print('leftwg: ', leftwg)
    plt.plot(x*pixelsize*0.001,leftwg, label='leftwg', color='tab:blue') # convert pixels to distance in µm
    plt.plot(x*pixelsize*0.001,rightwg, label='rightwg', color='tab:orange')
    if with_oscillation == True:
        p0 = [0.2, 50, 0.1, 0.29, 0, 0.04]
        # bounds = ([0, -100], [100, 100])
        coeff_leftwg, var_matrix = curve_fit(Damped_Exponenial_Function, x, leftwg, p0=p0)
        coeff_rightwg, var_matrix = curve_fit(Damped_Exponenial_Function, x, leftwg, p0=p0)
        # plt.plot(x, Damped_Oscillation_Function(x, 0.1, 0.29, 0, 0.04)+Exponenial_Function(x, 0.2, 0.01, 0))
        # plt.plot(x, Damped_Exponenial_Function(x, 0.2, 50, 0.1, 0.29, 0, 0.04), label='initial_guess')
        plt.plot(x, Damped_Exponenial_Function(x, coeff_leftwg[0], coeff_leftwg[1], coeff_leftwg[2], coeff_leftwg[3], coeff_leftwg[4], coeff_leftwg[5], ), label='fit_leftwg', color='tab:blue')
        plt.plot(x, Damped_Exponenial_Function(x, coeff_rightwg[0], coeff_rightwg[1], coeff_rightwg[2], coeff_rightwg[3], coeff_rightwg[4], coeff_rightwg[5], ), label='fit_rightwg', color='tab:orange')
        print('oscillation period leftwg = ', 2*np.pi/(np.sqrt(1-coeff_leftwg[5]**2 )*coeff_leftwg[3])*pixelsize, ' nm')
        print('oscillation period rightwg = ', 2*np.pi/(np.sqrt(1-coeff_rightwg[5]**2 )*coeff_rightwg[3])*pixelsize, ' nm')
    else:
        p0 = [1, 50*scaling]
        coeff_leftwg, var_matrix = curve_fit(Exponenial_Function, x, leftwg, p0=p0)
        coeff_rightwg, var_matrix = curve_fit(Exponenial_Function, x, rightwg, p0=p0)
        # plt.plot(x,leftwg, label='leftwg', color='tab:blue')
        # plt.plot(x,rightwg, label='rightwg', color='tab:orange')
        plt.plot(x*pixelsize*0.001, Exponenial_Function(x, coeff_leftwg[0], coeff_leftwg[1]), '--', label='fit_leftwg', color='tab:blue')
        plt.plot(x*pixelsize*0.001, Exponenial_Function(x, coeff_rightwg[0], coeff_rightwg[1]), '--', label='fit_rightwg', color='tab:orange')
    print('propagation length leftwg = ', coeff_leftwg[1]*pixelsize, ' nm') 
    print('propagation length rightwg = ', coeff_rightwg[1]*pixelsize, ' nm') 
    # print(coeff)
    plt.xlabel('$x \ [µm]$')
    plt.ylabel('$abs(E_z) \ [a.u.]$')
    plt.legend()
    plt.show()

def Damped_Exponenial_Function(x, A, d, B, omega_0, phase, gamma):
    return Exponenial_Function(x, A, d)+Damped_Oscillation_Function(x, B, omega_0, phase, gamma)

def Exponenial_Function(x, A, d):
    return A*np.exp(-(x)/(2*d)) # this is for the electric field, for intensity the /2 should be removed



def Damped_Oscillation_Function(t, A, omega_0, phase, gamma):
    '''Damped oscillation function
    
    Args:
        A [float]: amplitude of the oscillation
        period [float]: period of oscillation 
        phase [float]: initial phase of the oscillation
        gamma [float]: damping factor of the oscillation
    
    '''
    return A*np.exp(-gamma*omega_0*t)*np.sin(np.sqrt(1-gamma**2)*omega_0*t + phase)

def Fit_Linear_Function(data):
    p0 = [1, 0]
    bounds = ([0, -100], [100, 100])
    coeff, var_matrix = curve_fit(Linear_Function, np.arange(len(data)), data, p0=p0, bounds=bounds)
    return coeff

def Linear_Function(x, a, b):
    return a*x+b

def Double_Gauss(x, A, mu_1, mu_2, sigma, offset):
    return A*(np.exp(-(x-mu_1)**2/(2.*sigma**2)) + np.exp(-(x-mu_2)**2/(2.*sigma**2))) + offset

def Sawtooth_v2(x, n_periods, offset):
    phase_orientation = 1
    return np.pi*(1+signal.sawtooth(2*np.pi*n_periods * (x-offset), phase_orientation))

def Sawtooth(x, period, n_periods, offset):
    phase_orientation = 1
    return np.pi*(1+signal.sawtooth(period * n_periods * (x-offset), phase_orientation))


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

def Get_Phase_Difference_V2(value1, value2):
    sign = +1
    if value1 > value2:
        copy = value1
        value1 = value2
        value2 = copy
        # sign = -1
    diff1 = value2-value1
    diff2 = value1 + 2*np.pi-value2
    
    # if abs(diff1) <= abs(diff2):
    #     return sign*diff1
    # else:
    #     return sign*diff2
    return diff1




def main():
    Fit_Decaying_Amplitude()

if __name__ == '__main__':
    main()