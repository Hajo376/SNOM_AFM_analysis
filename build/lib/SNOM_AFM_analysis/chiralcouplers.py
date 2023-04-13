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
from datetime import datetime


from SNOM_AFM_analysis.lib.snom_colormaps import *
from SNOM_AFM_analysis.python_classes_snom import Open_Measurement
from SNOM_AFM_analysis.lib.phase_analysis import *
# from SNOM_AFM_analysis.lib.get_directionality import ChiralCoupler
from SNOM_AFM_analysis.lib.rectangle_selector import Select_Rectangle

class ChiralCouplers(Open_Measurement):
    """This class builds on the Open_Measurement base class for snom analysis.
    This class was initiated for a more detailed analysis of measurements of the chiral couplers.
    The basic idea is to find the waveguides by reducing the measurement and then fitting the waveguide positions.
    Once the positions of the waveguides are known we can get the amplitudes and phaseprofiles along the waveguides. 
    This will be used to calculate the directionality and optionally the phase difference.

    Availabe methods are:
        Find_Waveguides: Tries to find and fit the two waveguides, must be executed before any data can be extracted or used.
        Extract_WG_Data: This function will extract the data of the specified channel for the left and right waveguide.
            Requires Find_Waveguides to be executed beforehand.
        Extract_WG_Data_with_Errors: This function will extract the data of the specified channel for the left and right waveguide and the background in between.
            Requires Find_Waveguides to be executed beforehand.
        Export_All_Data: Exports the Data on the waveguides, averaged amplitude and averaged phase.
            Requires Find_Waveguides to be executed beforehand.
        Export_All_Data_single_phase_line: Exports the Data on the waveguides, averaged amplitude and phase in the center.
            Requires Find_Waveguides to be executed beforehand.
        Export_All_Data_single_phase_line_with_Errors: Exports the Data on the waveguides, averaged amplitude and phase in the center plus background in between.
            Requires Find_Waveguides to be executed beforehand.
        Fit_Sawtooth:

    """
    def __init__(self, directory_name:str, channels:list=None, title:str=None, autoscale:bool=True, orientation:str='vertical'):
        """Create a ChiralCouplers object. Inherits from Open_Measurement class and can therefore use all functions from that class.
        Only has additional functionality to analyse chiral couplers measurements or simulations.

        Args:
            directory_name (str): measurement folder, will be given to Open_Measurement class
            channels (list, optional): channels to get the profile data from, will be given to Open_Measurement class. Defaults to None.
            title (_type_, optional): title to be displayed, will be given to Open_Measurement class. Defaults to None.
            autoscale (bool, optional): autoscale data? Will be given to Open_Measurement class. Defaults to True.
            orientation (str, optional): orientation of the waveguides, 'vertical' or 'horizontal'. Defaults to 'vertical'.
        """
        self.orientation = orientation
        self.removed_lines = []
        super().__init__(directory_name, channels, title, autoscale)
    
    def Find_Waveguides(self, channels:list=None):
        self._Initialize_Data(channels)
        
        # first select the area to retieve the data from, should be a bit away from the coupler to ignore direct laser interactions
        # but not too far that the signal vanishes
        # print(self.height_channel)
        if self.height_channel in self.channels:
            height_data = self.all_data[self.channels.index(self.height_channel)]
        else:
            height_data = self._Load_Data([self.height_channel])[0][0]
        XRes = len(height_data[0])
        YRes = len(height_data)
        selection = Select_Rectangle(height_data, self.height_channel)
        self._Write_to_Logfile('Select_Rectangle', selection)
        # print('selection: ', selection)
        # set all data outside the selection to zero and apply auto cut feature to cut the data to specified area
        mask_array = np.zeros((YRes, XRes))
        for y in range(YRes):
            if (y > selection[0][1]) and (y < selection[1][1]):
                for x in range(XRes):
                    if (x > selection[0][0]) and (x < selection[1][0]):
                        mask_array[y][x] = 1
                      
        height_data = np.multiply(height_data, mask_array)
        height_data = self._Auto_Cut_Data(height_data)
        
        for channel in self.channels:
            self.all_data[self.channels.index(channel)] = np.multiply(self.all_data[self.channels.index(channel)], mask_array)
            self.all_data[self.channels.index(channel)] = self._Auto_Cut_Data(self.all_data[self.channels.index(channel)])
            if channel in self.amp_channels:
                amp_channel = channel
            elif (channel in self.phase_channels) or (channel in self.corrected_phase_channels):
                phase_channel = channel
        
        self.height_data = height_data
        self.data_array = self.all_data
        self.wg_positions = None
        # data used for displaying with the fit to the waveguides, typically the height profile is shown as this is the basis for the fits
        self.display_data = self.height_data
        # finally start the fitting procedure
        self._Fit_Wgs()
    
    def _Fit_Wgs(self, wg_height=65, max_shift=0.05) -> None:
        '''This fits the height channel with two gaussians to find the two waveguide positions for each row.
        The coefficients will be stored in an instance variable list_of_coefficients. The list will also be saved as .txt.
        The coefficients will be used later to extract the data from the measurements on the waveguides.
        '''
        # max shift is to limit jumping of the center point
        if self.orientation == 'vertical':
            width = len(self.height_data[0])
            max_shift = max_shift*width
            p0 = [wg_height, (width)/4, (width)/4*3, 5, 0]
            bounds = ([wg_height-15, 0, width/2, 2, -1000], [wg_height+15, width/2, width, 25, 1000])
            self.list_of_coefficients = []
            print('starting the fitting procedure')
            for row in range(len(self.height_data)):
                # check if row in removed areas
                if row in self.removed_lines:
                    self.list_of_coefficients.append(p0)
                else:
                    coeff, var_matrix = curve_fit(Double_Gauss, range(0, width), self.height_data[row], p0=p0, bounds=bounds)
                    # coeff = [A, mu_1, mu_2, sigma, offset]
                    self.list_of_coefficients.append(coeff)
                    p0 = coeff #set the starting parameters for the next fit
                    # bounds = ([wg_height-15, 0, width/2, 0, -1000], [wg_height+15, width/2, width, width/2, 1000])
                    # bounds = ([wg_height-15, p0[1]-max_shift, p0[2]-max_shift, 0, -1000], [wg_height+15, p0[1]+max_shift, p0[2]+max_shift, width/2, 1000])
                    # print(coeff)
            if len(self.removed_lines) != 0:
                # interpolate missing lines by linear function between first fit before and afterwards
                self._Fill_Removed_Lines()

            mean_sigma = [self.list_of_coefficients[i][3] for i in range(len(self.list_of_coefficients))]
            mean_sigma = np.mean(mean_sigma)
            self.integ_width = int(mean_sigma*2)
            self._Get_Fit_Feedback()
            self._Get_Integration_Width()
        # save list of coefficients
        # filepath = os.path.normpath(os.path.join(this_files_path, '../chiralcouplers_analysis/extracted_data/fit_coefficients'))# old and too specific
        filepath = os.path.normpath(os.path.join(self.directory_name, '../extracted_data/fit_coefficients'))
        file = open(filepath + '/' + self.filename + '.txt', 'w')
        now = datetime.now()
        current_datetime = now.strftime("%d/%m/%Y %H:%M:%S")
        file.write(current_datetime + '\n' + 'filename = ' + self.directory_name.split('/')[-1] + '\n')
        file.write('#\tA\tmu_1\tmu_2\tsigma\toffset\n')
        for i in range(len(self.list_of_coefficients)):
            file.write(f'{i}\t{round(self.list_of_coefficients[i][0], 3)}\t{round(self.list_of_coefficients[i][1], 3)}\t{round(self.list_of_coefficients[i][2], 3)}\t{round(self.list_of_coefficients[i][3], 3)}\t{round(self.list_of_coefficients[i][4], 3)}\n')
        file.close()

    def _Fill_Removed_Lines(self):
        jumps = self._Find_Jumps_in_Array(self.removed_lines)
        if len(jumps)==0:
            jumps = [len(self.removed_lines)]
        print('jumps: ', jumps)
        start = 0
        for jump in jumps:
            # get previous positions:
            left_prev = self.list_of_coefficients[self.removed_lines[start]-1][1]
            right_prev = self.list_of_coefficients[self.removed_lines[start]-1][2]
            print('left_prev=', left_prev)
            # get position after jump:
            left_after = self.list_of_coefficients[self.removed_lines[jump-1]+1][1]
            right_after = self.list_of_coefficients[self.removed_lines[jump-1]+1][2]
            print('left_after=', left_after)
            #determine interpolation function:
            # f(x)=mx+n
            m_left = (left_after-left_prev)/(jump-1 - start)
            m_right = (right_after-right_prev)/(jump-1 - start)
            m_mean = (m_left+m_right)/2
            print('m_mean: ', m_mean)
            print('start: ', start)
            print('jump: ', jump)

            def function_left(x):
                return m_mean*x + left_prev - m_mean*start
            
            def function_right(x):
                return m_mean*x + right_prev - m_mean*start

            for line in range(start, jump):
                left_mu = function_left(line)
                right_mu = function_right(line)
                self.list_of_coefficients[self.removed_lines[line]] = [0, left_mu, right_mu, 0, 0]
                print('Filled line: ', self.list_of_coefficients[self.removed_lines[line]])
            start = jump

    def _Find_Jumps_in_Array(self, array):
        '''finds jumps in ordered list of integers'''
        jumps = []
        last_element = array[0]
        for i in range(len(array)):
            if last_element < array[i]-1:
                jumps.append(i)
            last_element = array[i]
        jumps.append(len(array)-1) # append also the last index
        return jumps

    def _Select_Lines_to_Remove(self):
        # self.removed_areas = []
        self.removed_lines = []
        user_input = 'y'
        while user_input == 'y':
            selection = Select_Rectangle(self.height_data, self.height_channel)
            for i in range(selection[0][1], selection[1][1]):
                if i not in self.removed_lines:
                    self.removed_lines.append(i)
            user_input = input('Do you want to select an other area to remove? enter y/n: ')
        print('The removed areas are: ', self.removed_lines)

    def _Get_Fit_Feedback(self):
        self._Plot_Data_With_Fit_And_Bounds(self.height_data, 'Z')
        user_input = input('Are you happy with the current fit? y/n: ')
        if user_input == 'y':
            pass
        elif user_input == 'n':
            user_input = input('Do you want to retry the fit from the beginning? enter 0:\nDo you want to remove areas from the fit? enter 1: ')
            if user_input == '0':
                self.Find_Waveguides(self.channels)
                #operations like scaling and blurring will be lost...
            elif user_input == '1':
                self._Select_Lines_to_Remove()
                self._Fit_Wgs()
            else:
                print('Please enter n or y!')
                self._Get_Fit_Feedback()

        else:
            print('Please enter n or y!')
            self._Get_Fit_Feedback()

    def _Get_Integration_Width(self) -> None:
        self._Plot_Data_With_Fit_And_Bounds(self.height_data, 'Z')
        user_input = input(f'Do you like the extraction width? Current width is {self.integ_width} y/n: ')
        if user_input == 'n':
            self.integ_width = int(input('Please enter the new extraction width as an integer: '))
            self._Get_Integration_Width() # call recursively until user is happy with integ_width
        elif user_input == 'y':
            self._Write_to_Logfile('integ_width', self.integ_width)
        else:
            print('Please enter n or y!')
            self._Get_Integration_Width()

    def _Plot_Data_With_Fit_And_Bounds(self, data:np.array, channel:str) -> None:
        """This method will display the data and additionally show the center of the gaussian fits plus the chosen export width.

        Args:
            data (np.array): the data to display
            channel (str): specifies the colormap
        """
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

    def Extract_WG_Data(self, channel:str):
        """This function will extract the data of the specified channel for the left and right waveguide. 
        An instance variable containing both data sets is created (left_wg_data, right_wg_data).

        Args:
            channel (str): channel to extract the data from
        """
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
    
    def Extract_WG_Data_with_Errors(self, channel:str):
        """This function will extract the data of the specified channel for the left and right waveguide. 
        An instance variable containing both data sets is created (left_wg_data, right_wg_data).
        Additionally also the area in between the two waveguides is extracted as background_data.
        The background data will be used for an error estimation.

        Args:
            channel (str): channel to extract the data from
        """
        wg_positions = [[i, self.list_of_coefficients[i][1], self.list_of_coefficients[i][2]] for i in range(len(self.height_data))]
        # print('wg_positions', wg_positions)
        # print('integ_width: ', integ_width)
        self.left_wg_data = []
        self.right_wg_data = []
        self.background_data = []
        for row in range(len(wg_positions)):
            wg_left = round(wg_positions[row][1])
            wg_right = round(wg_positions[row][2])
            self.left_wg_data.append(self.data_array[self.channels.index(channel)][row][int(wg_left-self.integ_width):int(wg_left+self.integ_width)])
            self.right_wg_data.append(self.data_array[self.channels.index(channel)][row][int(wg_right-self.integ_width):int(wg_right+self.integ_width)])
            self.background_data.append(self.data_array[self.channels.index(channel)][row][int(wg_left+self.integ_width*2):int(wg_right-self.integ_width*2)])

    def _Data_Mean(self) -> None:
        self.mean_data = []
        # print('in _Data_Mean, print channels: ', self.channels)
        for channel in self.channels:
            if 'Z' not in channel:
                self.Extract_WG_Data(channel)
                left_mean_array = []
                right_mean_array = []
                if 'P' in channel:
                    self.left_center_phase = []
                    self.right_center_phase = []
                for row in range(len(self.data_array[self.channels.index(channel)])):
                    # create additional list only containing center phase data points
                    if 'P' in channel:
                        left_phase=self.left_wg_data[row][int(self.integ_width/2)]
                        right_phase=self.right_wg_data[row][int(self.integ_width/2)]
                        self.left_center_phase.append(left_phase)
                        self.right_center_phase.append(right_phase)
                    sum_left = 0
                    sum_right = 0
                    amp_left = self.left_wg_data[row]
                    amp_right = self.right_wg_data[row]
                    for j in range(2*self.integ_width):
                        sum_left+=amp_left[j]
                        sum_right+=amp_right[j]
                    left_mean_array.append(sum_left/(2*self.integ_width))
                    right_mean_array.append(sum_right/(2*self.integ_width))
                self.mean_data.append([channel, left_mean_array, right_mean_array])
    
    def _Data_Mean_with_Errors(self) -> None:
        """This function will create the mean values from the extracted amplitude data. The mean data is stored in a new instance variable mean_data.
        The phase data is just used to get the phase profile in the center of the waveguide, this leads to more reliable results, as averaging phase typically does not work out well.
        """
        self.mean_data = []
        # print('in _Data_Mean, print channels: ', self.channels)
        for channel in self.channels:
            if 'Z' not in channel:
                self.Extract_WG_Data_with_Errors(channel)
                left_mean_array = []
                right_mean_array = []
                background_mean_array = []
                if 'P' in channel:
                    self.left_center_phase = []
                    self.right_center_phase = []
                for row in range(len(self.data_array[self.channels.index(channel)])):
                    # create additional list only containing center phase data points
                    if 'P' in channel:
                        left_phase=self.left_wg_data[row][int(self.integ_width/2)]
                        right_phase=self.right_wg_data[row][int(self.integ_width/2)]
                        self.left_center_phase.append(left_phase)
                        self.right_center_phase.append(right_phase)
                    sum_left = 0
                    sum_right = 0
                    sum_background = 0
                    amp_left = self.left_wg_data[row]
                    amp_right = self.right_wg_data[row]
                    amp_background = self.background_data[row]
                    for j in range(2*self.integ_width):
                        sum_left+=amp_left[j]
                        sum_right+=amp_right[j]
                    for j in range(len(amp_background)):
                        sum_background+=amp_background[j]
                    left_mean_array.append(sum_left/(2*self.integ_width))
                    right_mean_array.append(sum_right/(2*self.integ_width))
                    background_mean_array.append(sum_background/len(amp_background))
                self.mean_data.append([channel, left_mean_array, right_mean_array, background_mean_array])

    def _Get_Data_Pairs(self) -> list:
        all_indices = []
        allowed_channel_types = ['A', 'P']
        data_pairs = []
        # print('Trying to find pairs')
        # for i in range(len(self.mean_data)):
            # print(self.mean_data[i][0])
        for i in range(len(self.mean_data)):
            channel = self.mean_data[i][0]
            # print('channel1: ', channel)
            channel.split()
            channel_type = channel[2]
            # print('channel type1: ', channel_type)
            index = channel[1] # basically the demodulation order of the channel, so 1 for O1A
            # print('index1: ', index)
            # channels must be sorted! amp comes first then phase
            if index not in all_indices and channel_type in allowed_channel_types and channel_type==allowed_channel_types[0]:
                for j in range(len(self.mean_data)):
                    channel = self.mean_data[j][0]
                    # print('channel2: ', channel)
                    channel.split()
                    new_channel_type = channel[2]
                    # print('channel type2: ', new_channel_type)
                    new_index = channel[1]
                    # print('index2: ', new_index)

                    if new_index == index and new_channel_type in allowed_channel_types and new_channel_type!=channel_type:
                        # print('found a pair!')
                        all_indices.append(index)
                        if len(self.mean_data[i]) == 4: # channel, left_mean_array, right_mean_array, background_mean_array
                            data_pairs.append([index , [self.mean_data[i][1], self.mean_data[i][2], self.mean_data[i][3]], [self.mean_data[j][1], self.mean_data[j][2]]]) # [index, [amp_leftwg, amp_rightwg, amp_background], [phase_leftwg, phase_rightwg]]
                        else:
                            data_pairs.append([index , [self.mean_data[i][1], self.mean_data[i][2]], [self.mean_data[j][1], self.mean_data[j][2]]]) # [index, [amp_leftwg, amp_rightwg], [phase_leftwg, phase_rightwg]]
        return data_pairs

    def Export_All_Data(self):
        self._Data_Mean() # create mean data
        data_pairs = self._Get_Data_Pairs() # get data pairs from the mean data list 
        for j in range(len(data_pairs)):
            filepath = os.path.normpath(os.path.join(this_files_path, '../chiralcouplers_analysis/extracted_data'))
            filename = f'{self.filename}_all_data_O{data_pairs[j][0]}.txt'
            filepath = filepath + '/' + filename
            with open(filepath, 'w') as datafile:
                datafile.write('#pixel\tamp leftwg\tamp rightwg\tphase leftwg\tphase rightwg\n')
                for i in range(len(data_pairs[j][1][0])):
                    datafile.write(f'{i}\t{round(data_pairs[j][1][0][i], 5)}\t{round(data_pairs[j][1][1][i], 5)}\t{round(data_pairs[j][2][0][i],3)}\t{round(data_pairs[j][2][1][i],3)}\n')
        print(f'Exported all data to {filepath}!')
    
    def Export_All_Data_single_phase_line(self, filepath=None):
        self._Data_Mean() # create mean data
        data_pairs = self._Get_Data_Pairs() # get data pairs from the mean data list 
        for j in range(len(data_pairs)):
            if filepath == None:
                filepath = os.path.normpath(os.path.join(this_files_path, '../chiralcouplers_analysis/extracted_data'))
            filename = f'{self.filename}_all_data_O{data_pairs[j][0]}_center_phase.txt'
            complete_filename = filepath + '/' + filename
            with open(complete_filename, 'w') as datafile:
                datafile.write('#pixel\tamp leftwg\tamp rightwg\tphase leftwg\tphase rightwg\n')
                for i in range(len(data_pairs[j][1][0])):
                    datafile.write(f'{i}\t{round(data_pairs[j][1][0][i], 5)}\t{round(data_pairs[j][1][1][i], 5)}\t{round(self.left_center_phase[i],3)}\t{round(self.right_center_phase[i],3)}\n')
        print(f'Exported all data to {complete_filename}!')

    def Export_All_Data_single_phase_line_with_Errors(self, filepath:str=None):
        """This method first extracts the waveguide data by selecting an area in the scan and then fitting the waveguide positions.
        Subsequently the extracted amplitude data is averaged, for the phase only the values in the middle of the waveguides are used.
        The averaged profiles and phase data are stored in data_pairs for each demodulation order. This is based on the channels specified earlier and requires that
        Find_Waveguides was executed prior to this method.

        Args:
            filepath (str, optional): where to save the extracted data? The filename is created automatically. Defaults to None.
        """
        self._Data_Mean_with_Errors() # create mean data
        data_pairs = self._Get_Data_Pairs() # get data pairs from the mean data list 
        for j in range(len(data_pairs)):
            if filepath == None:
                filepath = os.path.normpath(os.path.join(this_files_path, '../chiralcouplers_analysis/extracted_data'))
            filename = f'{self.filename}_all_data_O{data_pairs[j][0]}_center_phase_with_errors.txt'
            complete_filename = filepath + '/' + filename
            with open(complete_filename, 'w') as datafile:
                datafile.write('#pixel\tamp leftwg\tamp rightwg\tphase leftwg\tphase rightwg\tamp background\n')
                for i in range(len(data_pairs[j][1][0])):
                    datafile.write(f'{i}\t{round(data_pairs[j][1][0][i], 5)}\t{round(data_pairs[j][1][1][i], 5)}\t{round(self.left_center_phase[i],3)}\t{round(self.right_center_phase[i],3)}\t{round(data_pairs[j][1][2][i], 5)}\n')
        print(f'Exported all data to {complete_filename}!')


class LoadData():
    """This class is ment to load the profiles created with the ChiralCouplers class and comsol profiles.
    From the profiles the directionality can be calculated and compared.

    Available functions are:
        Load_Extracted_Data: Loads the extracted profiles for amp and phase in left and right wg and if with_errors==True also background.
        Cut_Data: Remove part of the profiles from end or beginning, only used for comsol data to remove oscillations in the beginning.
        Fit_Amplitudes: This fits the amplitude profiles and creates the propagation lenght and directionality. Subsequent load and fit calls will increment the lists.
        Plot_Propagation_Lengths: Displays the calculated propagation lengths.
        Plot_Directionality: Displays the directionalities.
        Export_Directionalities: Exports the directonalities to .txt file.
        Plot_Phase_Gradient: Displays the phase profiles in left and right wg in several ways like raw phase, corrected, phase difference and variation.
    """
    def __init__(self, initialdir:str=None, title:str=None, autoload:bool=True, chirality:int=+1, pixelsize:float=None, with_errors:bool=True):
        """This class is ment to load the profiles created with the ChiralCouplers class and comsol profiles.
        From the profiles the directionality can be calculated and compared.

        Args:
            initialdir (str, optional): initial directory, will be used for the file dialogues.
            title (str, optional): title, will be used for future plots.
            autoload (bool, optional): if set to True you will automatically be asked to select data files without having to call Load_Extracted_Data manually.
            chirality (int, optional): chirality of the measurement/simulation. Important for directionality calculation.
            pixelsize (float, optional): size of pixels in nm, important for plot labels, defaults to 50nm/scaling
            with_errors (bool, optional): do you want to include an error estimation based on the background in between the two waveguides? Only for snom measurements.
        """
        self.initialdir = initialdir
        self.chirality = chirality
        self.all_propagation_lengths = []
        self.all_directionalitys = []
        self.all_measured_directionalities = []
        self.all_mode_indices = []
        self.filename = None
        if pixelsize is None: # for the measurements
            scaling = 4
            pixelsize = 50/scaling # in nm 
        self.pixelsize = pixelsize
        self.with_errors = with_errors
        if autoload is True:
            self.Load_Extracted_Data()
        if title is not None:
            self.title = title

    def Load_Extracted_Data(self, filepath:str=None, initialdir:str=None) -> None:
        """Select the file to load the data from. File format: amp left, amp right, phase left, phase right and if with_errors==True also background, tab delimited, 2 lines header. 

        Args:
            filepath (str, optional): if not specified you will be prompted with a filedialoge. Defaults to None.
            initialdir (str, optional): sets the starting directory for the filedialoge. Defaults to None.
        """
        if initialdir is not None:
            self.initialdir = initialdir
        root = tk.Tk()
        root.withdraw()
        if filepath is None:
            if self.initialdir == None:
                filepath = filedialog.askopenfilename(initialdir=this_files_path)
            else:
                filepath = filedialog.askopenfilename(initialdir=os.path.normpath(os.path.join(this_files_path, self.initialdir)))
        self.filename = filepath.split('/')[-1]
        datafile = open(filepath, 'r')
        self.amplitude_left = np.genfromtxt(datafile, delimiter='\t', skip_header=2, usecols=1)
        datafile.close()
        datafile = open(filepath, 'r')
        self.amplitude_right = np.genfromtxt(datafile, delimiter='\t', skip_header=2, usecols=2)
        datafile.close()
        datafile = open(filepath, 'r')
        self.phase_left = np.genfromtxt(datafile, delimiter='\t', skip_header=2, usecols=3)
        datafile.close()
        datafile = open(filepath, 'r')
        self.phase_right = np.genfromtxt(datafile, delimiter='\t', skip_header=2, usecols=4)
        datafile.close()
        if self.with_errors == True:
            datafile = open(filepath, 'r')
            self.amp_background = np.genfromtxt(datafile, delimiter='\t', skip_header=2, usecols=5)
            datafile.close()
    
    def Cut_Data(self, length:float=4000, from_where='beginning'):
        """This method will cut the loaded profiles to a specific length.
        Only used for comsol data to remove the strong oscillations in the beginning.

        Args:
            length  (float, optional): length in nm to remove from the profiles. Careful, the definition of the pixelsize must be correct! Defaults to 4000.
            from_where (str, optional): from where to remove the specified length, can be either be 'beginning' to remove from the beginning or 'end' to remove from the end.
                Defaults to 'beginning'.
        """
        datapoints = len(self.amplitude_left)
        reduced_datapoinst = datapoints-int(length/self.pixelsize)
        amplitude_left = np.zeros(reduced_datapoinst)
        amplitude_right = np.zeros(reduced_datapoinst)
        phase_left = np.zeros(reduced_datapoinst)
        phase_right = np.zeros(reduced_datapoinst)
        count = 0
        for i in range(datapoints):
            if (from_where == 'beginning' and i>=length/self.pixelsize) or (from_where == 'end' and i<=(datapoints-length/self.pixelsize)):
                amplitude_left[count] = self.amplitude_left[i]
                amplitude_right[count] = self.amplitude_right[i]
                phase_left[count] = self.phase_left[i]
                phase_right[count] = self.phase_left[i]
                count += 1
        self.amplitude_left = amplitude_left
        self.amplitude_right = amplitude_right
        self.phase_left = phase_left
        self.phase_right = phase_right

    def Fit_Amplitudes(self, with_oscillation=False) -> None:
        """This function fits the amplitude data from the loaded profiles.
        The data is fitted with exponential function to extract the propagation length and the directionality.
        The directionalities are saved in the instance variable all_directionalitys. This contains the directionality and an error estimate.
        Repeated calls of Load_Extracted_Data and Fit_Amplitudes will create a list of directionalites, which can then be saved for a parameter variation.

        Args:
            with_oscillation (bool, optional): ties to fit the exponential additionally with a decaying oscillation. Defaults to False.
        """
        root = tk.Tk()
        root.withdraw()  
        x=np.arange(len(self.amplitude_left))      
        plt.plot(x*self.pixelsize*0.001,self.amplitude_left, label='leftwg', color='tab:blue') # convert pixels to distance in µm
        plt.plot(x*self.pixelsize*0.001,self.amplitude_right, label='rightwg', color='tab:orange')
        if with_oscillation == True:
            p0 = [0.2, 50, 0.1, 0.29, 0, 0.04]
            # bounds = ([0, -100], [100, 100])
            coeff_leftwg, var_matrix = curve_fit(Damped_Exponenial_Function, x, self.amplitude_left, p0=p0)
            coeff_rightwg, var_matrix = curve_fit(Damped_Exponenial_Function, x, self.amplitude_left, p0=p0)
            # plt.plot(x, Damped_Oscillation_Function(x, 0.1, 0.29, 0, 0.04)+Exponenial_Function(x, 0.2, 0.01, 0))
            # plt.plot(x, Damped_Exponenial_Function(x, 0.2, 50, 0.1, 0.29, 0, 0.04), label='initial_guess')
            plt.plot(x, Damped_Exponenial_Function(x, coeff_leftwg[0], coeff_leftwg[1], coeff_leftwg[2], coeff_leftwg[3], coeff_leftwg[4], coeff_leftwg[5], ), label='fit_leftwg', color='tab:blue')
            plt.plot(x, Damped_Exponenial_Function(x, coeff_rightwg[0], coeff_rightwg[1], coeff_rightwg[2], coeff_rightwg[3], coeff_rightwg[4], coeff_rightwg[5], ), label='fit_rightwg', color='tab:orange')
            print('oscillation period leftwg = ', 2*np.pi/(np.sqrt(1-coeff_leftwg[5]**2 )*coeff_leftwg[3])*self.pixelsize, ' nm')
            print('oscillation period rightwg = ', 2*np.pi/(np.sqrt(1-coeff_rightwg[5]**2 )*coeff_rightwg[3])*self.pixelsize, ' nm')
        else:
            # p0 = [1, 50*scaling, 0]
            p0 = [1, self.pixelsize, 0]
            coeff_leftwg, var_matrix = curve_fit(Exponenial_Function, x, self.amplitude_left, p0=p0)
            coeff_rightwg, var_matrix = curve_fit(Exponenial_Function, x, self.amplitude_right, p0=p0)
            print('offset leftwg: ', coeff_leftwg[2])
            print('offset rightwg: ', coeff_rightwg[2])
            # plt.plot(x,leftwg, label='leftwg', color='tab:blue')
            # plt.plot(x,rightwg, label='rightwg', color='tab:orange')
            plt.plot(x*self.pixelsize*0.001, Exponenial_Function(x, coeff_leftwg[0], coeff_leftwg[1], coeff_leftwg[2]), '--', label='fit_leftwg', color='tab:blue')
            plt.plot(x*self.pixelsize*0.001, Exponenial_Function(x, coeff_rightwg[0], coeff_rightwg[1], coeff_rightwg[2]), '--', label='fit_rightwg', color='tab:orange')
        print('propagation length leftwg = ', coeff_leftwg[1]*self.pixelsize, ' nm') 
        print('propagation length rightwg = ', coeff_rightwg[1]*self.pixelsize, ' nm') 
        self.all_propagation_lengths.append([coeff_leftwg[1]*self.pixelsize, coeff_rightwg[1]*self.pixelsize])
        # calculate area under exp. function from 0 to inf for left and right wg:
        area_left = 2*coeff_leftwg[0]*coeff_leftwg[1]
        area_right = 2*coeff_rightwg[0]*coeff_rightwg[1]
        if self.chirality == 1:
            directionality = (area_right**2)/(area_left**2 + area_right**2)
        elif self.chirality == -1:
            directionality = (area_left**2)/(area_left**2 + area_right**2)
        if self.with_errors == True:
            #sum all background amplitudes to get background area
            background = sum(self.amp_background) # for now the uncertainty for left and right is just the background sum
            if self.chirality == 1:
                directionality_error = Get_Directionality_Uncertainity(area_right, background, area_left, background)
            elif self.chirality == -1:
                directionality_error = Get_Directionality_Uncertainity(area_left, background, area_right, background)
        self.all_directionalitys.append([directionality, directionality_error])
        self._Get_Measured_Directionality()
        # print(coeff)
        plt.xlabel('$x \ [µm]$')
        plt.ylabel('$abs(E_z) \ [a.u.]$')
        plt.legend()
        plt.show()
    
    def _Get_Measured_Directionality(self):
        left_sum = np.sum(self.amplitude_left)
        right_sum = np.sum(self.amplitude_right)
        if self.chirality == 1:
            directionality = (right_sum**2)/(left_sum**2 + right_sum**2)
        elif self.chirality == -1:
            directionality = (left_sum**2)/(left_sum**2 + right_sum**2)
        if self.with_errors == True:
            #sum all background amplitudes to get background area
            background = sum(self.amp_background) # for now the uncertainty for left and right is just the background sum
            if self.chirality == 1:
                directionality_error = Get_Directionality_Uncertainity(right_sum, background, left_sum, background)
            elif self.chirality == -1:
                directionality_error = Get_Directionality_Uncertainity(left_sum, background, right_sum, background)
        self.all_measured_directionalities.append([directionality, directionality_error])

    def Plot_Propagation_Lengths(self, x_values:list, x_label:str) -> None:
        """This function will plot all propagation lengths in memory, only useful if multiple measurement profiles have been loaded and fitted.

        Args:
            x_values (list): a predefined list of x-values. This would for example be a list of radii when you want to display the data of a radius variation.
            x_label (_type_): label for the x-axis
        """
        
        propagation_lengths_leftwg = [element[0]*0.001 for element in self.all_propagation_lengths]
        propagation_lengths_rightwg = [element[1]*0.001 for element in self.all_propagation_lengths]
        
        plt.plot(x_values, propagation_lengths_leftwg, 'x', color='tab:blue', label='leftwg')
        plt.plot(x_values, propagation_lengths_rightwg, 'x', color='tab:orange', label='rightwg')
        plt.title(self.title)
        plt.xlabel(x_label)
        plt.ylabel('Propagation Length [µm]')
        plt.legend()
        plt.show()
        # fig, ax = plt.subplots()
    
    def Plot_Directionality(self, x_values, x_label) -> None:
        if self.with_errors == True:
            measured_directionalities = [x[0] for x in self.all_measured_directionalities]
            yerr = [x[1] for x in self.all_measured_directionalities]
            plt.errorbar(x_values, measured_directionalities, yerr=yerr, color='tab:blue', label='from data')
            fit_directionalities = [x[0] for x in self.all_directionalitys]
            yerr = [x[0] for x in self.all_directionalitys]
            plt.errorbar(x_values, fit_directionalities, yerr=yerr, color='tab:orange', label='from fit')
        else: 
            plt.plot(x_values, self.all_measured_directionalities, 'x', color='tab:blue', label='from data')
            plt.plot(x_values, self.all_directionalitys, 'x', color='tab:orange', label='from fit')
        plt.title(self.title)
        plt.xlabel(x_label)
        plt.ylabel('Directionality [a.u.]')
        plt.legend()
        plt.show()

    def Export_Directionalities(self, datapoints:list, filepath:str, variation_type:str='radius', data_type:str='snom'):
        """This method will export the created list of directionalities to a .txt file.

        Args:
            datapoints (list): list of the parameters which were varied in this measurement sequence
            filepath (str): complete filepath with '.txt' ending for the directionalities list.
            variation_type (str, optional): will be included in file header. Defaults to 'radius'.
            data_type (str, optional): either 'snom' or 'comsol'. Defaults to 'snom'. For 'snom' the directionalities from the exponential fit will be used,
                for 'comsol' the sum of the data will be used, since sometimes the field is so weak that the fit does not work properly.
        """
        file = open(filepath, 'w')
        if self.with_errors == True:
            file.write(f'{variation_type}\tfit_directionality\tuncertainty\tsum_directionality\tuncertainty\n')
            for i in range(len(datapoints)):
                if data_type == 'snom':
                    file.write(f'{datapoints[i]}\t{self.all_directionalitys[i][0]}\t{self.all_directionalitys[i][1]}\t{self.all_measured_directionalities[i][0]}\t{self.all_measured_directionalities[i][1]}\n') # export values from fitcomparison, because data is very noisy
                elif data_type == 'comsol':
                    file.write(f'{datapoints[i]}\t{self.all_measured_directionalities[i]}\n') # export amplitude comparison, because fit for low intensity wg does not work properly
            file.close()
        else:
            file.write(f'{variation_type}\tdirectionality\n')
            for i in range(len(datapoints)):
                if data_type == 'snom':
                    file.write(f'{datapoints[i]}\t{self.all_directionalitys[i]}\n') # export values from fitcomparison, because data is very noisy
                elif data_type == 'comsol':
                    file.write(f'{datapoints[i]}\t{self.all_measured_directionalities[i]}\n') # export amplitude comparison, because fit for low intensity wg does not work properly
            file.close()

    def Plot_Phase_Gradient(self):
        """Displays the phase in several ways like raw, corrected, flattened with fit, phase difference between left and right wg and varaiations.
        """
        x_values = np.arange(len(self.phase_left))
        plt.plot(x_values*self.pixelsize*0.001, self.phase_left, 'x', color='tab:blue', label='leftwg')
        plt.plot(x_values*self.pixelsize*0.001, self.phase_right, 'x', color='tab:orange', label='rightwg')
        plt.xlabel('X [µm]')
        plt.ylabel('Phase')
        plt.title(self.filename)
        plt.legend()
        plt.show()

        phase_left_cleaned = self._Remove_Artificial_Phase(self.phase_left)
        phase_right_cleaned = self._Remove_Artificial_Phase(self.phase_right)
        plt.plot(x_values*self.pixelsize*0.001, phase_left_cleaned, 'x', color='tab:blue', label='leftwg cleaned + flattened')
        plt.plot(x_values*self.pixelsize*0.001, phase_right_cleaned, 'x', color='tab:orange', label='rightwg cleaned + flattened')
        # plt.title(self.title)
        plt.xlabel('X [µm]')
        plt.ylabel('Phase')
        plt.legend()
        plt.show()
        difference = Get_Profile_Difference(phase_left_cleaned, phase_right_cleaned)
        plt.plot(x_values*self.pixelsize*0.001, difference, '.', label='Phase difference')
        plt.title(self.filename)
        plt.xlabel('X [µm]')
        plt.ylabel('Phase difference')
        mean_difference = np.mean(difference)
        if mean_difference < np.pi/2: 
            plt.ylim((-np.pi, np.pi))
        elif mean_difference >= np.pi/2:
            plt.ylim((0, np.pi*2))
        plt.legend()
        plt.show()
        # fit flattened phase profiles
        p0 = [1, 0]
        # bounds = ([0, -100], [100, 100])
        coeff_left, pcov_left = Fit_Linear_Function(phase_left_cleaned, p0)
        coeff_right, pcov_right = Fit_Linear_Function(phase_right_cleaned, p0)
        # perr_left = np.sqrt(np.diag(pcov_left))
        # perr_right = np.sqrt(np.diag(pcov_right))
        plt.plot(x_values*self.pixelsize*0.001, phase_left_cleaned, 'x', color='tab:blue', label='leftwg cleaned + flattened')
        plt.plot(x_values*self.pixelsize*0.001, phase_right_cleaned, 'x', color='tab:orange', label='rightwg cleaned + flattened')
        plt.plot(x_values*self.pixelsize*0.001, Linear_Function(x_values, coeff_left[0], coeff_left[1]), color='tab:blue', label='leftwg fit')
        plt.plot(x_values*self.pixelsize*0.001, Linear_Function(x_values, coeff_right[0], coeff_right[1]), color='tab:orange', label='rightwg fit')
        plt.title(self.filename)
        plt.xlabel('X [µm]')
        plt.ylabel('Phase')
        plt.legend()
        plt.show()
        phase_differece_fit = Linear_Function(x_values, coeff_left[0], coeff_left[1]) - Linear_Function(x_values, coeff_right[0], coeff_right[1])
        plt.plot(x_values*self.pixelsize*0.001, phase_differece_fit, '.', color='tab:blue', label='Phase difference from fit')
        # phase_differece_fit = coeff_left-coeff_right
        # plt.hlines(phase_differece_fit, color='tab:blue', label='Phase difference from fit')
        plt.plot(x_values*self.pixelsize*0.001, difference, '.', color='tab:orange', label='Phase difference')
        # fit difference with linear function and compare to phase_differece_fit
        coeff_difference, pcov_difference = Fit_Linear_Function(difference)
        # perr_difference = np.sqrt(np.diag(pcov_difference))
        plt.plot(x_values*self.pixelsize*0.001, Linear_Function(x_values, coeff_difference[0], coeff_difference[1]), color='tab:red', label='Phase difference fit')
        print('The mean phase difference is: ', np.mean(phase_differece_fit))

        plt.title(self.filename)
        plt.xlabel('X [µm]')
        plt.ylabel('Phase difference')
        mean_difference = np.mean(difference)
        if mean_difference < np.pi/2: 
            plt.ylim((-np.pi, np.pi))
        elif mean_difference >= np.pi/2:
            plt.ylim((0, np.pi*2))
        plt.legend()
        plt.show()

        #mode index:
        mode_index_left = Get_Modeindex_from_LinearFunction(coeff_left[0], self.pixelsize, wavelength=1600)
        mode_index_right = Get_Modeindex_from_LinearFunction(coeff_right[0], self.pixelsize, wavelength=1600)
        print('The mode index of the left wg is: ', mode_index_left)
        print('The mode index of the right wg is: ', mode_index_right)

        # plot phase gradient minus fitfunction to look for oscillations...
        plt.plot(x_values*self.pixelsize*0.001, phase_left_cleaned-Linear_Function(x_values, coeff_left[0], coeff_left[1]), color='tab:blue', label='left wg')
        plt.plot(x_values*self.pixelsize*0.001, phase_right_cleaned-Linear_Function(x_values, coeff_right[0], coeff_right[1]), color='tab:orange', label='right wg')
        plt.title(self.filename)
        plt.xlabel('X [µm]')
        plt.ylabel('Phase variation')
        plt.legend()
        plt.show()

    def _Remove_Artificial_Phase(self, phasearray, bounds=1.7):
        new_array = np.zeros((len(phasearray)))
        previous_element = phasearray[0]
        phaseshift = 0
        for i in range(len(phasearray)):
            element = phasearray[i]
            if previous_element >= np.pi*2-bounds:
                if element < bounds:
                    phaseshift += np.pi*2
                    previous_element = element
                    # print('in a')
                elif element >= previous_element:
                    previous_element = element
                    # print('in b')
                else:
                    element = np.pi*2# set all phase values on the falling flank to 2pi, as they should not exist 
                    # print('in c')
            else:
                previous_element = element
            new_array[i] = element + phaseshift
        return new_array

    

def Flatten_Profile(data):
    threshold = np.pi*0.8
    flattened_profile = []
    phase_shift = 0
    for i in range(len(data)):
        if data[i] < (data[i-1]-threshold):
            phase_shift += np.pi*2
        flattened_profile.append(data[i]+phase_shift)
    return flattened_profile

def Damped_Exponenial_Function(x, A, d, B, omega_0, phase, gamma):
    return Exponenial_Function(x, A, d)+Damped_Oscillation_Function(x, B, omega_0, phase, gamma)

def Exponenial_Function(x, A, d, offset):
    offset = 0#for now
    return A*np.exp(-(x)/(2*d))+offset # this is for the electric field, for intensity the /2 should be removed



def Damped_Oscillation_Function(t, A, omega_0, phase, gamma):
    '''Damped oscillation function
    
    Args:
        A [float]: amplitude of the oscillation
        period [float]: period of oscillation 
        phase [float]: initial phase of the oscillation
        gamma [float]: damping factor of the oscillation
    
    '''
    return A*np.exp(-gamma*omega_0*t)*np.sin(np.sqrt(1-gamma**2)*omega_0*t + phase)

def Fit_Linear_Function(data, p0=None, bounds=None):
    if p0 != None:
        if bounds!= None:
            coeff, var_matrix = curve_fit(Linear_Function, np.arange(len(data)), data, p0=p0, bounds=bounds)
        else:
            coeff, var_matrix = curve_fit(Linear_Function, np.arange(len(data)), data, p0=p0)
    else:
        coeff, var_matrix = curve_fit(Linear_Function, np.arange(len(data)), data)
    return coeff, var_matrix

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

def Get_Directionality_Uncertainity(high_amp, high_amp_uncert, low_amp, low_amp_uncert):
    return (2*high_amp)/((high_amp**2 + low_amp**2)**2)*np.sqrt((low_amp**2*high_amp_uncert)**2 + (high_amp*low_amp*low_amp_uncert)**2) # simply gaussian error propagation


# unused?
'''
    def Fit_Decaying_Amplitude(self):
        with_oscillation = False
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
        plt.plot(x*self.pixelsize*0.001,leftwg, label='leftwg', color='tab:blue') # convert pixels to distance in µm
        plt.plot(x*self.pixelsize*0.001,rightwg, label='rightwg', color='tab:orange')
        if with_oscillation == True:
            p0 = [0.2, 50, 0.1, 0.29, 0, 0.04]
            # bounds = ([0, -100], [100, 100])
            coeff_leftwg, var_matrix = curve_fit(Damped_Exponenial_Function, x, leftwg, p0=p0)
            coeff_rightwg, var_matrix = curve_fit(Damped_Exponenial_Function, x, leftwg, p0=p0)
            # plt.plot(x, Damped_Oscillation_Function(x, 0.1, 0.29, 0, 0.04)+Exponenial_Function(x, 0.2, 0.01, 0))
            # plt.plot(x, Damped_Exponenial_Function(x, 0.2, 50, 0.1, 0.29, 0, 0.04), label='initial_guess')
            plt.plot(x, Damped_Exponenial_Function(x, coeff_leftwg[0], coeff_leftwg[1], coeff_leftwg[2], coeff_leftwg[3], coeff_leftwg[4], coeff_leftwg[5], ), label='fit_leftwg', color='tab:blue')
            plt.plot(x, Damped_Exponenial_Function(x, coeff_rightwg[0], coeff_rightwg[1], coeff_rightwg[2], coeff_rightwg[3], coeff_rightwg[4], coeff_rightwg[5], ), label='fit_rightwg', color='tab:orange')
            print('oscillation period leftwg = ', 2*np.pi/(np.sqrt(1-coeff_leftwg[5]**2 )*coeff_leftwg[3])*self.pixelsize, ' nm')
            print('oscillation period rightwg = ', 2*np.pi/(np.sqrt(1-coeff_rightwg[5]**2 )*coeff_rightwg[3])*self.pixelsize, ' nm')
        else:
            p0 = [1, self.pixelsize]
            coeff_leftwg, var_matrix = curve_fit(Exponenial_Function, x, leftwg, p0=p0)
            coeff_rightwg, var_matrix = curve_fit(Exponenial_Function, x, rightwg, p0=p0)
            # plt.plot(x,leftwg, label='leftwg', color='tab:blue')
            # plt.plot(x,rightwg, label='rightwg', color='tab:orange')
            plt.plot(x*self.pixelsize*0.001, Exponenial_Function(x, coeff_leftwg[0], coeff_leftwg[1]), '--', label='fit_leftwg', color='tab:blue')
            plt.plot(x*self.pixelsize*0.001, Exponenial_Function(x, coeff_rightwg[0], coeff_rightwg[1]), '--', label='fit_rightwg', color='tab:orange')
        print('propagation length leftwg = ', coeff_leftwg[1]*self.pixelsize, ' nm') 
        print('propagation length rightwg = ', coeff_rightwg[1]*self.pixelsize, ' nm') 
        # print(coeff)
        plt.xlabel('$x \ [µm]$')
        plt.ylabel('$abs(E_z) \ [a.u.]$')
        plt.legend()
        plt.show()

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


    def Display_Phase_Profile(self, channel):
        if channel in self.channels:
            self._Plot_Data_With_Fit_And_Bounds(self.data_array[self.channels.index(channel)], channel)
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
        
        # coeff = self.Fit_Sawtooth(phase_right_integ_array)
        # x = range(len(phase_left_integ_array))
        # x = np.linspace(0, 1, 500)
        # plt.plot(x*len(phase_left_integ_array), Sawtooth_v2(x, coeff[0], coeff[1]), label='Fit to right wg')
        
        flattened_profile = Flatten_Profile(self.phase_right_integ_array)
        p0 = [1, 0]
        bounds = ([0, -100], [100, 100])
        coeff, pcov = Fit_Linear_Function(flattened_profile, p0, bounds)
        plt.plot(range(len(self.phase_right_integ_array)), flattened_profile, label='flattened phase')
        plt.plot(range(len(self.phase_right_integ_array)), Linear_Function(np.arange(len(self.phase_right_integ_array)), coeff[0], coeff[1]), label='fit: flattened phase')
        plt.title(f'{channel}')
        plt.legend()
        plt.show()

    

    def Integrate_Amplitude(self, channel):
        if channel in self.channels:
            self._Plot_Data_With_Fit_And_Bounds(self.data_array[self.channels.index(channel)], channel)
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

    def Export_Mean_Data(self, filename):
        self._Data_Mean()
        for i in range(len(self.mean_data)):
            channel = self.mean_data[i][0]
            self.Export_Data(channel, filename, self.mean_data[i][1], self.mean_data[i][2])
    
    def Export_Mean_Data_with_Errors(self, filename):
        self._Data_Mean_with_Errors()
        for i in range(len(self.channels)):
            channel = self.mean_data[i][0]
            self.Export_Data(channel, filename, self.mean_data[i][1], self.mean_data[i][2])

    def Load_Extracted_Data(self):
        root = tk.Tk()
        root.withdraw()
        filepath = filedialog.askopenfilename(initialdir=os.path.normpath(os.path.join(this_files_path, '../chiralcouplers_analysis/extracted_data')))
        datafile = open(filepath, 'r')
        self.amplitude_left = np.genfromtxt(datafile, delimiter='\t', skip_header=2, usecols=1)
        datafile.close()
        datafile = open(filepath, 'r')
        self.amplitude_right = np.genfromtxt(datafile, delimiter='\t', skip_header=2, usecols=2)
        datafile.close()
        datafile = open(filepath, 'r')
        self.phase_left = np.genfromtxt(datafile, delimiter='\t', skip_header=2, usecols=3)
        datafile.close()
        datafile = open(filepath, 'r')
        self.phase_right = np.genfromtxt(datafile, delimiter='\t', skip_header=2, usecols=4)
        datafile.close()

    # should be included as external functionality, plot profiles, 
    def Display_Phase_Difference(self, channel):
        if channel in self.channels:
            self._Plot_Data_With_Fit_And_Bounds(self.data_array[self.channels.index(channel)], channel)
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

'''