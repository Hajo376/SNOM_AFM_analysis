##########################################################################
# This code was created by Hans-Joachim Schill, University of Bonn, 2022 #
##########################################################################
from ast import Try
from operator import lt
from scipy.ndimage import gaussian_filter
from struct import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_point_clicker import clicker# used for getting coordinates from images
from matplotlib_scalebar.scalebar import ScaleBar # used for creating scale bars
from enum import Enum, auto


# from python_manipulation_classes.get_bounds import Get_Reduced_Data  # legacy method based on pyqtgraph
from .lib.rectangle_selector import Select_Rectangle
from .lib.phase_analysis import Flatten_Phase_Profile, Get_Profile_Difference
from .lib.snom_colormaps import *

# Definitions for plotting
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
'''
SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 20
'''
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

class Definitions(Enum):
    vertical = auto()
    horizontal = auto()

class File_Type(Enum):
    """Different file types may vary accross different platforms."""
    standard = auto()
    standard_new = auto()
    aachen_ascii = auto()
    aachen_gsf = auto()
    # for the parameters:
    html = auto()
    txt = auto()
    html_new = auto()# new software version creates slightly different html file

class Tag_Type(Enum):
    # Tag types:
    scan_type = auto()
    center_pos = auto()
    rotation = auto()
    scan_area = auto()
    pixel_area = auto()
    integration_time = auto()
    tip_frequency = auto()
    tip_amplitude = auto()
    tapping_amplitude = auto()

class File_Definitions:
    #standard file definitions
    file_type = File_Type.standard
    parmeters_type = File_Type.html

class Plot_Definitions:
    hide_ticks = True
    figsizex = 10
    figsizey = 5
    show_titles = True
    tight_layout = True
    colorbar_width = 10 # in percent, standard is 5 or 10
    hspace = 0.4 #standard is 0.4

class SimpleManipulation(File_Definitions, Plot_Definitions):
    
    # reorganize and put the following variables in the init function
    # ToDo
    all_subplots = []
    mask_array = []
    scalebar = [] # specifiy the channels to add the scalebar to, should contain the channel plus the scalbar object in a nested list
    # piezo motion distortion correction:
    upper_y_bound = None
    lower_y_bound = None
    align_points = None
    y_shifts = None
    scaling_factor = 1

    def __init__(self, directory_name, channels:list=None, title=None, autoscale=True) -> None:
        '''
        autoscale will automatically scale the data such that each pixel is quadratic in dimensions.
        '''
        self.directory_name = directory_name
        self.filename = directory_name.split('/')[-1]
        self.measurement_title = title # If a measurement_title is specified it will precede the automatically created title based on the channel dictionary
        self.autoscale = autoscale
        self.logfile_path = self._Initialize_Logfile()
        self._Initialize_File_Type()
        if channels == None: # the standard channels which will be used if no channels are specified
            channels = [self.preview_ampchannel, self.preview_phasechannel, self.height_channel]
        self.channels = channels
        self._Create_Tag_Dict()
        self._Create_Channels_Tag_Dict()
        self.XRes, self.YRes = self.measurement_tag_dict[Tag_Type.pixel_area]
        self.XReal, self.YReal = self.measurement_tag_dict[Tag_Type.scan_area] # in µm
        # Create a variable containing the data of the specified channels which can be varied later on with e.g. the gaussian_filter and subsequent methods
        # Aditionally a dictionary is created which contains the channel information
        # The dictionary and self.channels are not identical as self.channels only contains the raw channel information, whereas the dictionary is modified to contain
        # the information of which modifications like scaling, filtering ... have been applied.
        self.all_data, self.all_data_dict = self._Load_Data(self.channels)
        if autoscale == True:
            self.Quadratic_Pixels()

    def _Find_Filetype(self) -> None:
        '''This function aims at finding specific characteristics in the filename to idendify the filetype.
        For example the difference in File_Type.standard and File_Type.standard_new are an additional ' raw' at the end of the filename.'''
        try:
            f=open(f"{self.directory_name}/{self.filename} O1A.gsf","br")
        except:
            # filetype is at least not standard
            try:
                f=open(f"{self.directory_name}/{self.filename} O1A raw.gsf","br")
            except:
                print("The correct filetype could automatically be found. Please try again and specifiy the filetype.")
            else:
                f.close()
                File_Definitions.file_type = File_Type.standard_new
                File_Definitions.parmeters_type = File_Type.html_new
        else:
            f.close()
            File_Definitions.file_type = File_Type.standard
            File_Definitions.parmeters_type = File_Type.html

    def _Initialize_File_Type(self) -> None:
        self._Find_Filetype() # try to find the filetype automatically
        self.file_type = File_Definitions.file_type
        self.parameters_type = File_Definitions.parmeters_type
        if self.file_type == File_Type.standard or self.file_type == File_Type.standard_new:
            self.height_channel = "Z C"
            self.all_channels = ["O1A","O1P","O2A","O2P","O3A","O3P","O4A","O4P","O5A","O5P"]
            self.phase_channels = ['O1P','O2P','O3P','O4P','O5P']
            self.corrected_phase_channels = ['O1P_corrected','O2P_corrected','O3P_corrected','O4P_corrected','O5P_corrected']
            self.amp_channels = ['O1A','O2A','O3A','O4A','O5A']
            self.real_channels = ['O1R', 'O2R', 'O3R', 'O4R', 'O5R']
            self.preview_ampchannel = 'O2A'
            self.preview_phasechannel = 'O2P'
            self.height_indicator = 'Z'
            self.amp_indicator = 'A'
            self.phase_indicator = 'P'
        elif self.file_type == File_Type.aachen_ascii or self.file_type == File_Type.aachen_gsf:
            self.all_channels = ["O1-F-abs","O1-F-arg","O2-F-abs","O2-F-arg","O3-F-abs","O3-F-arg","O4-F-abs","O4-F-arg"]
            self.phase_channels = ['O1-F-arg','O2-F-arg','O3-F-arg','O4-F-arg']
            self.amp_channels = ['O1-F-abs','O2-F-abs','O3-F-abs','O4-F-abs']
            self.real_channels = ['O1-F-real','O2-F-real','O3-F-real','O4-F-real',]
            self.height_channel = "MT-F-abs"
            self.preview_ampchannel = 'O2-F-abs'
            self.preview_phasechannel = 'O2-F-arg'
            self.height_indicator = 'MT'
            self.amp_indicator = 'abs'
            self.phase_indicator = 'arg'

    def _Create_Tag_Dict(self):
        # create tag_dict for each channel individually? if manipulated channels are loaded they might have different diffrent resolution
        # only center_pos, scan_area, pixel_area and rotation must be stored for each channel individually but rotation is not stored in the original .gsf files
        # but rotation could be added in the newly created .gsf files
        if self.parameters_type == File_Type.html:
            all_tables = pd.read_html("".join([self.directory_name,"/",self.filename,".html"]))
            tables = all_tables[0]
            self.measurement_tag_dict = {
                Tag_Type.scan_type: tables[2][0],
                Tag_Type.center_pos: [float(tables[2][4]), float(tables[3][4])],
                Tag_Type.rotation: int(tables[2][5]),
                Tag_Type.scan_area: [float(tables[2][6]), float(tables[3][6])],
                Tag_Type.pixel_area: [int(tables[2][7]), int(tables[3][7])],
                Tag_Type.integration_time: float(tables[2][9]),
                Tag_Type.tip_frequency: float(tables[2][13]),
                Tag_Type.tip_amplitude: float(tables[2][14]),
                Tag_Type.tapping_amplitude: float(tables[2][15])
            }
        elif self.parameters_type == File_Type.html_new:
            all_tables = pd.read_html("".join([self.directory_name,"/",self.filename,".html"]))
            tables = all_tables[0]
            self.measurement_tag_dict = {
                Tag_Type.scan_type: tables[2][0],
                Tag_Type.center_pos: [float(tables[2][4]), float(tables[3][4])],
                Tag_Type.rotation: int(tables[2][5]),
                Tag_Type.scan_area: [float(tables[2][6]), float(tables[3][6])],
                Tag_Type.pixel_area: [int(tables[2][7]), int(tables[3][7])],
                Tag_Type.integration_time: float(tables[2][9]),
                Tag_Type.tip_frequency: float(tables[2][14]),
                Tag_Type.tip_amplitude: float(tables[2][15]),
                Tag_Type.tapping_amplitude: float(tables[2][16])
            }
        elif self.parameters_type == File_Type.txt:

            parameters = self.directory_name + '/' + self.filename + '.parameters.txt'
            file = open(parameters, 'r')
            parameter_list = file.read()
            file.close()
            # print(parameter_list)
            parameter_list = parameter_list.split('\n')
            parameter_list = [element.split(': ') for element in parameter_list]
            center_pos = [float(parameter_list[7][1]), float(parameter_list[8][1])]
            rotation = float(parameter_list[9][1])
            scan_area = [float(parameter_list[0][1]), float(parameter_list[1][1])]
            pixel_area = [int(parameter_list[3][1]), int(parameter_list[4][1])]
            integration_time = float(parameter_list[6][1])
            tip_frequency = float(parameter_list[10][1])
            self.measurement_tag_dict = {
                Tag_Type.scan_type: None,
                Tag_Type.center_pos: center_pos,
                Tag_Type.rotation: rotation,
                Tag_Type.scan_area: scan_area,
                Tag_Type.pixel_area: pixel_area,
                Tag_Type.integration_time: integration_time,
                Tag_Type.tip_frequency: tip_frequency,
                Tag_Type.tip_amplitude: None,
                Tag_Type.tapping_amplitude: None
            }

    def _Create_Channels_Tag_Dict(self):
        if (self.file_type == File_Type.standard) or (self.file_type == File_Type.standard_new) or (self.file_type == File_Type.aachen_gsf):
            cod="latin1"
            # get the tag values from each .gsf file individually
            self.channel_tag_dict = []
            for channel in self.channels:
                if self.file_type == File_Type.standard_new and '_corrected' not in channel:
                    if channel == 'Z C':
                        filepath = self.directory_name + '/' + self.filename + ' ' + channel + '.gsf'
                    else:
                        filepath = self.directory_name + '/' + self.filename + ' ' + channel + ' raw.gsf'
                else:
                    filepath = self.directory_name + '/' + self.filename + ' ' + channel + '.gsf'
                file = open(filepath, 'r', encoding=cod)
                content = file.read()
                file.close()
                XRes = self._Get_Tagval(content, 'XRes')
                YRes = self._Get_Tagval(content, 'YRes')
                XReal = self._Get_Tagval(content, 'XReal')
                YReal = self._Get_Tagval(content, 'YReal')
                XOffset = self._Get_Tagval(content, 'XOffset')
                YOffset = self._Get_Tagval(content, 'YOffset')
                Rotation = 0
                try:
                    Rotation = self._Get_Tagval(content, 'Rotation')
                except:
                    print('Could not find the rotation tag value, proceeding anyways.')
                channel_dict = {
                    Tag_Type.center_pos: [float(XOffset), float(YOffset)],
                    Tag_Type.rotation: float(Rotation),
                    Tag_Type.pixel_area: [int(XRes), int(YRes)],
                    Tag_Type.scan_area: [float(XReal), float(YReal)] 
                }
                self.channel_tag_dict.append(channel_dict)
            pass
        else:
            if self.parameters_type == File_Type.txt:
                self.channel_tag_dict = []
                parameters = self.directory_name + '/' + self.filename + '.parameters.txt'
                file = open(parameters, 'r')
                parameter_list = file.read()
                file.close()
                # print(parameter_list)
                parameter_list = parameter_list.split('\n')
                parameter_list = [element.split(': ') for element in parameter_list]
                center_pos = [float(parameter_list[7][1]), float(parameter_list[8][1])]
                rotation = float(parameter_list[9][1])
                scan_area = [float(parameter_list[0][1]), float(parameter_list[1][1])]
                pixel_area = [int(parameter_list[3][1]), int(parameter_list[4][1])]
                channel_dict = {
                        Tag_Type.center_pos: center_pos,
                        Tag_Type.rotation: rotation,
                        Tag_Type.pixel_area: pixel_area,
                        Tag_Type.scan_area: scan_area
                    }
                # for this file type all channels must be of same size
                for channel in self.channels:
                    self.channel_tag_dict.append(channel_dict)
            else:
                print('channel tag dict for this filetype is not yet implemented')

    def _Get_Tagval(self, content, tag):
        """This function gets the value of the tag listed in the file header"""
        taglength=len(tag)
        tagdatastart=content.find(tag)+taglength+1   
        tagdatalength=content[tagdatastart:].find("\n") 
        tagdataend=tagdatastart+tagdatalength
        tagval=float(content[tagdatastart:tagdataend])
        return tagval

    def _Initialize_Data(self, channels) -> None:
        '''This function initializes the data in memory. If no channels are specified the already existing data is used,
        which is created automatically in the instance init method. If channels are specified, the instance data is overwritten.
        Channels must be specified as a list of channels.'''
        if channels == None:
            #none means the channels specified in the instance creation should be used
            pass
        else:
            self.all_data, self.all_data_dict = self._Load_Data(channels)
            self.channels = channels
            # update the channel tag dictionary, makes the program compatible with differrently sized datasets, like original data plus manipulated, eg. cut data
            self._Create_Channels_Tag_Dict()
            # reset all the instance variables dependent on the data, but nor the ones responsible for plotting
            self.scaling_factor = 1
            if self.autoscale == True:
                self.Quadratic_Pixels()
            self.mask_array = [] # not shure if it's best to reset the mask...
            self.upper_y_bound = None
            self.lower_y_bound = None
            self.align_points = None
            self.y_shifts = None
    
    def Initialize_Channels(self, channels:list) -> None:
        '''This function will load the data from the specified channels and replace the ones in memory.'''
        self._Initialize_Data(channels)

    def _Initialize_Logfile(self) -> str:
        logfile_path = self.directory_name + '/python_manipulation_log.txt'
        file = open(logfile_path, 'a') # the new logdata will be appended to the existing file
        now = datetime.now()
        current_datetime = now.strftime("%d/%m/%Y %H:%M:%S")
        file.write(current_datetime + '\n' + 'filename = ' + self.filename + '\n')
        file.close()
        return logfile_path

    def _Write_to_Logfile(self, parameter_name:str, parameter):
        file = open(self.logfile_path, 'a')
        file.write(f'{parameter_name} = {parameter}\n')
        file.close()

    def Print_test(self) -> None:
        '''Testfunction to print several instance values.'''
        # print('test:')
        # print('self.directory_name: ', self.directory_name)
        # print('self.filename: ', self.filename)
        # print('self.scaling_factor_y: ', self.scaling_factor_y)
        # print('self.scaling_factor_x: ', self.scaling_factor_x)
        # print('self.XRes: ', self.XRes)
        # print('self.YRes: ', self.YRes)
        print('self.channels: ', self.channels)
        print('self.all_subplots[-1]: ', [element[3] for element in self.all_subplots])

    def _Scale_Array(self, array, scaling) -> np.array:
        '''This function scales a given 2D Array, it thus creates 'scaling'**2 subpixels per pixel.
        The scaled array is returned.'''
        yres = len(array)
        xres = len(array[0])
        scaled_array = np.zeros((yres*scaling, xres*scaling))
        for i in range(len(array)):
            # line = np.zeros((xres))
            for j in range(len(array[0])):
                for k in range(scaling):
                    for l in range(scaling):
                        # scaled_dataset[h][i*scaling + k][j*scaling + l] = array[i][j]
                        scaled_array[i*scaling + k][j*scaling + l] = array[i][j]
        return scaled_array

    def Scale_Channels(self, channels:list=None, scaling:int=4) -> None:
        '''This function scales all the data in memory or the specified channels.
                
        Args:
            channels [list]: list of channels, will override the already existing channels
            scaling [int]: defines scaling factor per axis. Each pixel will be scaled to scaling**2 pixels.
        '''
        # ToDo: reimplement scaling dependent on axis, x and y independently
        self._Initialize_Data(channels)
        self._Write_to_Logfile('scaling', scaling)
        self.scaling_factor = scaling
        dataset = self.all_data
        yres = len(dataset[0])
        xres = len(dataset[0][0])
        # self.all_data = np.zeros((len(dataset), yres*scaling, xres*scaling))
        self.all_data = []
        for h in range(len(dataset)):
            array = dataset[h]
            scaled_array = self._Scale_Array(array, scaling)
            self.all_data.append(scaled_array)

    def _Gauss_Blurr_Data(self, array, sigma) -> np.array:
        '''Applies a gaussian blurr to the specified array, with a specified sigma. The blurred data is returned as a list.'''
        return gaussian_filter(array, sigma)

    def _Load_Data(self, channels) -> list:
        '''Loads all binary data of the specified channels and returns them in a list plus the dictionary with the channel information.
        Height data is automatically converted to nm. '''
        # datasize=int(self.XRes*self.YRes*4)
        #create a list containing all the lists of the individual channels
        all_binary_data = []
        #safe the information about which channel is which list in a dictionary
        data_dict = []
        # data_dict = {}
        # all_data = np.zeros((len(channels), self.YRes, self.XRes))
        all_data = []
        if self.file_type==File_Type.standard or self.file_type==File_Type.standard_new:
            for i in range(len(channels)):
                print(channels[i])
                if self.file_type==File_Type.standard_new and '_corrected' not in channels[i]:
                    if channels[i] == 'Z C':
                        f=open(f"{self.directory_name}/{self.filename} {channels[i]}.gsf","br")
                    else:
                        f=open(f"{self.directory_name}/{self.filename} {channels[i]} raw.gsf","br")
                else:
                    f=open(f"{self.directory_name}/{self.filename} {channels[i]}.gsf","br")
                binarydata=f.read()
                f.close()
                all_binary_data.append(binarydata)
                data_dict.append(channels[i])
            count = 0
            for channel in channels:
                # print(self.channel_tag_dict[self.channels.index(channel)])
                # print(self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area])
                # dont remember why i added the following but it leads to problems
                # if channel in self.channels:
                #     XRes, YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area]
                # else:
                XRes, YRes = self.measurement_tag_dict[Tag_Type.pixel_area]
                
                datasize=int(XRes*YRes*4)
                channel_data = np.zeros((YRes, XRes))
                reduced_binarydata=all_binary_data[count][-datasize:]
                phaseoffset = 0
                rounding_decimal = 2
                scaling = 1
                if 'A' in channel:
                    rounding_decimal = 5
                if 'Z' in channel:
                    print('channel: ', channel)
                    scaling = 1000000000#convert to nm
                if 'P' in channel:
                    # normal phase data ranges from -pi to pi and gets shifted by +pi
                    phaseoffset = np.pi
                    if '_corrected' in channel:
                        # if the data is from a corrected channel it is already shifted
                        phaseoffset = 0
                if 'R' in channel:
                    rounding_decimal = 4
                for y in range(0,YRes):
                    for x in range(0,XRes):
                        pixval=unpack("f",reduced_binarydata[4*(y*XRes+x):4*(y*XRes+x+1)])[0]
                        channel_data[y][x] = round(pixval*scaling + phaseoffset, rounding_decimal)
                all_data.append(channel_data)
                count+=1

        elif self.file_type == File_Type.aachen_ascii:
            for i in range(len(channels)):
                file = open(f"{self.directory_name}/{self.filename}_{channels[i]}.ascii", 'r')
                string_data = file.read()
                datalist = string_data.split('\n')
                datalist = [element.split(' ') for element in datalist]
                datalist = np.array(datalist[:-1], dtype=np.float)# convert list to np.array and strings to float
                channel = channels[i]
                phaseoffset = 0
                rounding_decimal = 2
                scaling = 1
                if ('abs' in channel) and ('MT' not in channel):
                    rounding_decimal = 5
                if 'arg' in channel:
                    phaseoffset = np.pi
                    flattened_data = datalist.flatten()# ToDo just for now, in future the voltages have to be converted
                    phaseoffset = min(flattened_data)
                    if '_corrected' in channel:
                        # if the data is from a corrected channel it is already shifted
                        phaseoffset = 0
                if 'MT' in channel:
                    scaling = pow(10, 9)
                XRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area][0]
                YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area][1]
                for y in range(0,YRes):
                    for x in range(0,XRes):
                        datalist[y][x] = round(datalist[y][x]*scaling + phaseoffset, rounding_decimal)
                all_data[i] = datalist
                data_dict.append(channels[i])

        elif self.file_type == File_Type.aachen_gsf:
            for i in range(len(channels)):
                f=open(f"{self.directory_name}/{self.filename}_{channels[i]}.gsf","br")
                binarydata=f.read()
                f.close()
                all_binary_data.append(binarydata)
                data_dict.append(channels[i])
            scaling = 1
            count = 0
            for channel in channels:
                XRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area][0]
                YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area][1]
                datasize=int(XRes*YRes*4)
                reduced_binarydata=all_binary_data[count][-datasize:]
                phaseoffset = 0
                rounding_decimal = 2
                if ('abs' in channel) and ('MT' not in channel):
                    rounding_decimal = 5
                    scaling = pow(10, 6)
                if 'arg' in channel:
                    phaseoffset = np.pi
                    
                    if '_corrected' in channel:
                        # if the data is from a corrected channel it is already shifted
                        phaseoffset = 0
                if 'MT' in channel:

                    # scaling_factor = 1000000000
                    scaling = pow(10, 9)
                    rounding_decimal = 5
                    # scaling_factor = 1
                    # print()
                for y in range(0,YRes):
                    for x in range(0,XRes):
                        pixval=unpack("f",reduced_binarydata[4*(y*XRes+x):4*(y*XRes+x+1)])[0]
                        all_data[count][y][x] = round(pixval*scaling + phaseoffset, rounding_decimal)
                count+=1
        return all_data, data_dict

    def _Load_Data_Binary(self, channels) -> list:
        '''Loads all binary data of the specified channels and returns them in a list plus the dictionary for access'''
        #create a list containing all the lists of the individual channels
        all_binary_data = []
        #safe the information about which channel is which list in a dictionary
        data_dict = []
        for i in range(len(channels)):
            f=open(f"{self.directory_name}/{self.filename} {channels[i]}.gsf","br")
            binarydata=f.read()
            f.close()
            all_binary_data.append(binarydata)
            data_dict.append(channels[i])
        return all_binary_data, data_dict

    def Set_Min_to_Zero(self, channles) -> None:
        '''This function sets the min value of the specified channels to zero.
                
        Args:
            channels [list]: list of channels, will override the already existing channels
        '''
        self._Write_to_Logfile('set_min_to_zero', True)
        for channel in channles:
            data = self.all_data[self.channels.index(channel)]
            flattened_data = data.flatten()
            data_min = min(flattened_data)
            self.all_data[self.channels.index(channel)] = data - data_min

    def _Add_Subplot(self, data, channel, scalebar=None) -> list:
        '''This function adds the specified data to the list of subplots. The list of subplots contains the data, the colormap,
        the colormap label and a title, which are generated from the channel information. The same array is also returned,
        so it can also be iterated by an other function to only plot the data of interest.'''
        if self.file_type == File_Type.standard or self.file_type == File_Type.standard_new:
            if 'A' in channel:
                cmap=SNOM_amplitude
                label = 'amplitude [a.u.]'
                title = f'Amplitude {channel}'
            elif 'P' in channel:
                if 'positive' in channel:
                    cmap = SNOM_phase
                    title = f'positively corrected phase O{channel[1]}P'
                elif 'negative' in channel:
                    cmap = SNOM_phase
                    title = f'negatively corrected phase O{channel[1]}P'
                else:
                    cmap=SNOM_phase
                    title = f'Phase {channel}'
                label = 'phase'
            elif 'Z' in channel:
                cmap='gray'
                label = 'height [nm]'
                title = f'Height {channel}'
            elif 'R' in channel:
                cmap=SNOM_realpart
                label = 'E [a.u.]'
                title = f'Realpart {channel}'
        elif self.file_type == File_Type.aachen_ascii or self.file_type == File_Type.aachen_gsf:
            if 'abs' in channel and not 'MT' in channel:
                cmap=SNOM_amplitude
                label = 'amplitude [a.u.]'
                title = f'Amplitude {channel}'
            elif 'arg' in channel:
                if 'positive' in channel:
                    cmap = SNOM_phase
                    title = f'positively corrected phase O{channel[1]}P' # ToDo
                elif 'negative' in channel:
                    cmap = SNOM_phase
                    title = f'negatively corrected phase O{channel[1]}P'
                else:
                    cmap=SNOM_phase
                    title = f'Phase {channel}'
                label = 'phase'
            elif 'MT' in channel:
                cmap='gray'
                label = 'height [nm]'
                title = f'Height {channel}'
        elif 'fft' in channel:
            cmap='viridis'
            label = 'intensity [a.u.]'
            title =  f'Fourier Transform {channel}'
        elif 'gauss' in channel:
            title = f'Gauss blurred {channel}'
        else:
            print('Unknown channel')
            exit()
        # subplots.append([data, cmap, label, title])
        if self.measurement_title != None:
            title = self.measurement_title + title
        if scalebar != None:
            self.all_subplots.append([np.copy(data), cmap, label, title, scalebar])
            return [data, cmap, label, title, scalebar]
        else:
            self.all_subplots.append([np.copy(data), cmap, label, title])
            return [data, cmap, label, title]
    
    def Remove_Subplots(self, index_array:list) -> None:
        '''This function removes the specified subplot from the memory.
        
        Args:
            index_array [list]: The indices of the subplots to remove from the plot list
        '''
        #sort the index array in descending order and delete the corresponding plots from the memory
        index_array.sort(reverse=True)
        for index in index_array:
            del self.all_subplots[index]

    def Remove_Last_Subplots(self, times:int=1) -> None:
        '''This function removes the last added subplots from the memory.
        Times specifies how often the last subplot should be removed.
        Times=1 means only the last, times=2 means the two last, ...
        
        Args:
            times [int]: how many subplots should be removed from the end of the list?
        '''
        for i in range(times):
            self.all_subplots.pop()

    def _Plot_Subplots(self, subplots) -> None:
        '''This function plots the subplots. The plots are created in a grid, by default the grid is optimized for 3 by 3. The layout changes dependent on the number of subplots
        of subplots and also the dimensions. Wider subplots are prefferably created vertically, otherwise they are plotted horizontally. Probably subject to future changes...'''
        number_of_axis = 9
        number_of_subplots = len(subplots)
        # print('Number of subplots: {}'.format(number_of_subplots))
        #specify the way the subplots are organized
        nrows = int((number_of_subplots-1)/np.sqrt(number_of_axis))+1

        if nrows >=2:
            ncols = int(np.sqrt(number_of_axis))
        elif nrows == 1:
            ncols = number_of_subplots
        else:
            print('Number of subplots must be lager than 0!')
            exit()
        changed_orientation = False
        if number_of_subplots == 4:
            ncols = 2
            nrows = 2
            changed_orientation = True
        data = subplots[0][0]
        # calculate the ratio (x/y) of the data, if the ratio is larger than 1 the images are wider than high,
        # and they will prefferably be positiond vertically instead of horizontally
        ratio = len(data[0])/len(data)
        if ((number_of_subplots == 2) or (number_of_subplots == 3)) and ratio >= 2:
            nrows = number_of_subplots
            ncols = 1
            changed_orientation = True
        #create the figure with subplots
        fig, ax = plt.subplots(nrows, ncols)    
        fig.set_figheight(self.figsizey)
        fig.set_figwidth(self.figsizex) 
        counter = 0

        
        for row in range(nrows):
            for col in range(ncols):
                if counter < number_of_subplots:
                    if nrows == 1:
                        if ncols == 1:
                            axis = ax
                        else:
                            axis = ax[col]
                    elif ncols == 1 and (nrows == 2 or nrows == 3) and changed_orientation == True:
                        axis = ax[row]
                    else:
                        axis = ax[row, col]
                    data = subplots[counter][0]
                    cmap = subplots[counter][1]
                    label = subplots[counter][2]
                    title = subplots[counter][3]
                    if len(subplots[counter]) == 5:
                        dx, units, dimension, scalebar_label, length_fraction, height_fraction, width_fraction, location, loc, pad, border_pad, sep, frameon, color, box_color, box_alpha, scale_loc, label_loc, font_properties, label_formatter, scale_formatter, fixed_value, fixed_units, animated, rotation = subplots[counter][4]
                        scalebar = ScaleBar(dx, units, dimension, scalebar_label, length_fraction, height_fraction, width_fraction,
                            location, loc, pad, border_pad, sep, frameon, color, box_color, box_alpha, scale_loc,
                            label_loc, font_properties, label_formatter, scale_formatter, fixed_value, fixed_units, animated, rotation) 
                        axis.add_artist(scalebar)
                        print('added a scalebar')
                    if 'R_corrected' in title:
                        flattened_data = data.flatten()
                        min_real = min(flattened_data)
                        max_real = max(flattened_data)
                        print('min: ', min_real)
                        print('max: ', max_real)
                        if abs(min_real) > abs(max_real):
                            limit = abs(min_real)
                        else: limit = abs(max_real)
                        print('limit: ', limit)
                        img = axis.pcolormesh(data, cmap=cmap, vmin=-limit, vmax=limit)
                    else:
                        img = axis.pcolormesh(data, cmap=cmap)
                    if ('Z' in title) and ('_masked' in title) and ('_reduced' not in title):
                        # create a white border around the masked area, but show the full unmasked height data
                        border_width = 1
                        yres = len(data)
                        xres = len(data[0])
                        white_pixels = np.zeros((yres, xres))
                        for y in range(border_width, yres - border_width):
                            for x in range(border_width, xres - border_width):
                                mean = self._Get_Mean_Value(self.mask_array, x, y, border_width)
                                if (self.mask_array[y][x] == 0) and (0 < mean) and (mean < 1):
                                    white_pixels[y, x] = 100
                        # The idea is to plot a second pcolormesh on the same axis as the height data
                        # Only the pixels with a nonzero value are displayed, all other are set to zero opacity (alpha value)
                        ncolors = 2
                        color_array = plt.get_cmap('Greys')(range(ncolors))

                        # change alpha values
                        color_array[:,-1] = np.linspace(0.0,1.0,ncolors)

                        # create a colormap object
                        map_object = LinearSegmentedColormap.from_list(name='rainbow_alpha',colors=color_array)

                        # register this new colormap with matplotlib
                        plt.register_cmap(cmap=map_object)
                        axis.pcolormesh(white_pixels, cmap='rainbow_alpha')
                    
                    
                    # elif '_shifted' in title:
                    #     XRes = len(data[0])
                    #     axis.plot(self.align_points, [element +int((self.upper_y_bound - self.lower_y_bound)/2) for element in self.y_shifts], color='red')
                    #     axis.hlines([self.upper_y_bound, self.lower_y_bound], xmin=0, xmax=XRes, color='white')

                    # axis = ax[col][row]
                    # invert y axis to fit to the scanning procedure which starts in the top left corner
                    axis.invert_yaxis()
                    # ratio = len(data[0])/len(data)
                    divider = make_axes_locatable(axis)
                    cax = divider.append_axes("right", size=f"{self.colorbar_width}%", pad=0.05) # size is the size of colorbar relative to original axis, 100% means same size, 10% means 10% of original
                    cbar = plt.colorbar(img, aspect=1, cax=cax)
                    cbar.ax.get_yaxis().labelpad = 15
                    cbar.ax.set_ylabel(label, rotation=270)
                    # print('label: ', label)
                    if self.hide_ticks == True:
                        # remove ticks on x and y axis, they only show pixelnumber anyways, better to add a scalebar
                        axis.set_xticks([])
                        axis.set_yticks([])
                    # adjust the colorbar range for realpart images, such that 0 is in the middle
                    # if 'R_corrected' in title:
                    #     flattened_data = data.flatten()
                    #     min_real = min(flattened_data)
                    #     max_real = max(flattened_data)
                    #     if abs(min_real) > abs(max_real):
                    #         limit = min_real
                    #     else: limit = max_real
                    #     cbar.set_clim(-limit, limit)
                    if self.show_titles == True:
                        axis.set_title(title)
                    axis.axis('scaled')
                    counter += 1
                    # add scalebar:
                    # ToDo
                    '''print('title: ', title)
                    print('channel: ', self.scalebar[counter][0])
                    if self.scalebar[counter][0] in title:
                        scalebar = self.scalebar[counter][1]
                        axis.add_artist(scalebar)'''

        #turn off all unneeded axis
        counter = 0
        for row in range(nrows):
            for col in range(int(np.sqrt(number_of_axis))):
                if  counter >= number_of_subplots > 3 and changed_orientation == False and number_of_subplots != 4: 
                    axis = ax[row, col]
                    axis.axis('off')
                counter += 1

        plt.subplots_adjust(hspace=self.hspace)
        if self.tight_layout == True:
            plt.tight_layout()
        plt.show()
    
    def Switch_Supplots(self, first_id:int=None, second_id:int=None) -> None:
        '''
        This function changes the position of the subplots.
        The first and second id corresponds to the positions of the two subplots which should be switched.
        This function will be repea you are satisfied.

        Args:
            first_id [int]: the first id of the two subplots which should be switched
            second_id [int]: the second id of the two subplots which should be switched
        '''
        if (first_id == None) or (second_id == None):
            first_id = int(input('Please enter the id of the first image: '))
            second_id = int(input('Please enter the id of the second image: '))
        first_subplot = self.all_subplots[first_id]
        self.all_subplots[first_id] = self.all_subplots[second_id]
        self.all_subplots[second_id] = first_subplot
        self.Display_All_Subplots()
        print('Are you happy with the new positioning?')
        user_input = self._User_Input_Bool()
        if user_input == False:
            print('Do you want to change the order again?')
            user_input = self._User_Input_Bool()
            if user_input == False:
                exit()
            else:
                self.Switch_Supplots()

    def _Display_Dataset(self, dataset, channels) -> None:
        '''Add all data contained in dataset as subplots to one figure.
        The data has to be shaped beforehand!
        channels should contain the information which channel is stored at which position in the dataset.
        '''
        subplots = []
        for i in range(len(dataset)):
            scalebar = None
            for j in range(len(self.scalebar)):
                if self.channels[i] == self.scalebar[j][0]:
                    scalebar = self.scalebar[j][1]
            subplots.append(self._Add_Subplot(dataset[i], channels[i], scalebar))
        self._Plot_Subplots(subplots)

    def Display_All_Subplots(self) -> None:
        '''
        This function displays all the subplots which have been created until this point.
        '''
        self._Plot_Subplots(self.all_subplots)
        
    def Display_Channels(self, channels:list=None) -> None:
        '''This function displays the channels in memory or the specified ones without any modifications
                
        Args:
            channels [list]: list of channels, will override the already existing channels

        '''
        if channels == None:
            dataset = self.all_data
            channels = self.all_data_dict
        else:
            dataset, dict = self._Load_Data(channels)
        self._Display_Dataset(dataset, channels)

    def Gauss_Filter_Channels(self, channels:list=None, sigma=2):
        self._Initialize_Data(channels)
        self._Write_to_Logfile('gaussian_filter_sigma', sigma)
        if self.scaling_factor == 1:
            print('The data is not yet scaled! Do you want to scale the data?')
            user_input = self._User_Input_Bool()
            if user_input == True:
                self.Scale_Channels()
        # start the blurring:
        for i in range(len(self.channels)):
            self.all_data[i] = self._Gauss_Blurr_Data(self.all_data[i], sigma)
            self.all_data_dict[i] += '_gauss'
    
    def _Create_Header(self, channel, data=None, filetype='gsf'):
        # data = self.all_data[self.channels.index(channel)]
        # load data instead, because sometimes the channel is not in memory
        if data is None:
            # channel is not in memory, so the standard values will be used
            data = self._Load_Data([channel])[0][0]
            XReal, YReal = self.measurement_tag_dict[Tag_Type.scan_area]
            rotation = self.measurement_tag_dict[Tag_Type.rotation]
            XOffset, YOffset = self.measurement_tag_dict[Tag_Type.center_pos]
        else: 
            # if channel is in memory it has to have a channel dict, where all necessary infos are stored
            XReal, YReal = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.scan_area]
            rotation = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.rotation]
            XOffset, YOffset = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.center_pos]
        XRes = len(data[0])
        YRes  = len(data)
        # print('xres: ', XRes)
        # print('yres: ', YRes)
        if filetype=='gsf':
            header = f'Gwyddion Simple Field 1.0\n'
        elif filetype=='txt':
            header = 'Simple Textfile, floats seperated by spaces\n'
        # round everything to nm
        header += f'XRes={round(XRes,9)}\nYRes={round(YRes,9)}\n'
        header += f'XReal={round(XReal,9)}\nYReal={round(YReal,9)}\n'
        header += f'XYUnits=m\n'
        header += f'XOffset={round(XOffset,9)}\nYOffset={round(YOffset,9)}\n'
        header += f'Rotation={round(rotation)}\n'
        if self.height_indicator in channel:
            header += 'ZUnits=nm\n'
        else:
            header += 'ZUnits=\n'
        header += f'Title={self.measurement_title}\n'
        # lenght = header.count('\n')
        length = len(header)
        number = 4 - ((length) % 4)
        NUL = b'\0'
        for i in range(number -1):
            NUL += b'\0' # add NUL terminator
        return header, NUL

    def Save_to_gsf(self, channels:list=None, appendix:str='_manipulated'):
        '''This function is ment to save all specified channels to external .gsf files.
        
        Args:
            channels [list]:    list of the channels to be saved, if not specified, all channels in memory are saved
                                Careful! The data will be saved as it is right now, so with all the manipulations.
                                Therefor the data will have an '_manipulated' appendix in the filename.
        '''
        self._Write_to_Logfile('save_to_gsf_appendix', appendix)
        if channels == None:
            channels = self.channels
        for channel in channels:
            filepath = self.directory_name + '/' + self.filename + ' ' + channel + appendix + '.gsf'
            data = self.all_data[self.channels.index(channel)]
            XRes = len(data[0])
            YRes  = len(data)
            header, NUL = self._Create_Header(channel, data)
            file = open(filepath, 'bw')
            file.write(header.encode('utf-8'))
            file.write(NUL) # the NUL marks the end of the header and konsists of 0 characters in the first dataline
            for y in range(YRes):
                for x in range(XRes):
                    file.write(pack('f', round(data[y][x], 2)))
            file.close()
            print(f'successfully saved channel {channel} to .gsf')
            print(filepath)

    def Save_to_txt(self, channels:list=None, appendix:str='_manipulated'):
        '''This function is ment to save all specified channels to external .gsf files.
        
        Args:
            channels [list]:    list of the channels to be saved, if not specified, all channels in memory are saved
                                Careful! The data will be saved as it is right now, so with all the manipulations.
                                Therefor the data will have an '_manipulated' appendix in the filename.
        '''
        self._Write_to_Logfile('save_to_txt_appendix', appendix)
        if channels == None:
            channels = self.channels
        for channel in channels:
            filepath = self.directory_name + '/' + self.filename + ' ' + channel + appendix + '.txt'
            data = self.all_data[self.channels.index(channel)]
            XRes = len(data[0])
            YRes  = len(data)
            header, NUL = self._Create_Header(channel, data, 'txt')
            file = open(filepath, 'w')
            file.write(header)
            # file.write(NUL) # the NUL marks the end of the header and konsists of 0 characters in the first dataline
            for y in range(YRes):
                for x in range(XRes):
                    file.write(f'{round(data[y][x], 2)} ')
            file.close()
            print(f'successfully saved channel {channel} to .txt')
            print(filepath)
    
    def _Get_Channel_Scaling(self, channel_id) -> int :
        '''This function checks if an instance channel is scaled and returns the scaling factor.'''
        channel_yres = len(self.all_data[channel_id])
        # channel_xres = len(self.all_data[channel_id][0])
        return int(channel_yres/self.YRes)

    def _Create_Height_Mask_Preview(self, mask_array) -> None:
        '''This function creates a preview of the height masking.
        The preview is based on all channels in the instance'''
        channels = self.channels
        dataset = self.all_data
        subplots = []
        for i in range(len(dataset)):
            masked_array = np.multiply(dataset[i], mask_array)
            subplots.append(self._Add_Subplot(np.copy(masked_array), channels[i]))
        self._Plot_Subplots(subplots)
        # remove the preview subplots from the memory
        self.Remove_Last_Subplots(3)
        
    def _User_Input_Bool(self) -> bool: 
        '''This function asks the user to input yes or no and returns a boolean value.'''
        user_input = input('Please type y for yes or n for no. \nInput: ')
        if user_input == 'y':
            user_bool = True
        elif user_input == 'n':
            user_bool = False
        return user_bool

    def _Create_Mask_Array(self, height_data, threshold) -> np.array:
        '''This function takes the height data and a threshold value to create a mask array containing 0 and 1 values.
        '''
        height_flattened = height_data.flatten()
        height_threshold = threshold*(max(height_flattened)-min(height_flattened))+min(height_flattened)

        # create an array containing 0 and 1 depending on wether the height value is below or above threshold
        mask_array = np.copy(height_data)
        yres = len(height_data)
        xres = len(height_data[0])
        for y in range(yres):
            for x in range(xres):
                value = 0
                if height_data[y][x] >= height_threshold:
                    value = 1
                mask_array[y][x] = value
        return mask_array

    def _Get_Height_Treshold(self, height_data, mask_array, threshold) -> float:
        '''This function returns the height threshold value dependent on the user input'''
        self._Create_Height_Mask_Preview(mask_array)
        print('Do you want to use these parameters to mask the data?')
        mask_data = self._User_Input_Bool()
        if mask_data == False:
            print('Do you want to change the treshold?')
            change_treshold = self._User_Input_Bool()
            if change_treshold == True:
                print(f'The old threshold was {threshold}')
                threshold = float(input('Please enter the new treshold value: '))
                mask_array = self._Create_Mask_Array(height_data, threshold)
                self._Get_Height_Treshold(height_data, mask_array, threshold)
            else:
                print('Do you want to abort the masking procedure?')
                abort = self._User_Input_Bool()
                if abort == True:
                    exit()
        return threshold

    def Heigth_Mask_Channels(self, channels:list=None, mask_channel=None, threshold=0.5, mask_data=False, export:bool=False) -> None:
        '''
        The treshold factor should be between 0 and 1. It sets the threshold for the height pixels.
        Every pixel below threshold will be set to 0. This also applies for all other channels. 
        You can either specify specific channels to mask or if you don't specify channels,
        all standard channels will be masked. If export is False only the channels in self.channels will be masked
        and nothing will be exported. 
        For this function to also work with scaled data the height channel has to be specified and scaled as well!
                
        Args:
            channels [list]: list of channels, will override the already existing channels
            threshold [float]: threshold value to create the height mask from
            mask_data [bool]: if you want to apply the mask directly with the specified threshold change to 'True',
                            otherwise you will be prompted with a preview and can then change the threshold iteratively
            export [bool]: if you want to apply the mask to all channels and export them change to 'True'
        '''
        if export == True:
            channels = self.all_channels
        self._Initialize_Data(channels)
        if (mask_channel == None) or (mask_channel not in self.channels):
            if self.height_channel in self.channels:
                height_data = self.all_data[self.channels.index(self.height_channel)]
                if 'leveled' not in self.all_data_dict[self.channels.index(self.height_channel)]:
                    leveled_height_data = self._Height_Levelling_3Point(height_data)
                else:
                    leveled_height_data = height_data # since height_data is already leveled
                
                # if the height channel is used in the instance for example with gaussian blurr, this data will be used and not the raw data
                # because the data might be scaled, careful to always use the corrected height channel 'Z C'
            else:
                height_data, trash = self._Load_Data([self.height_channel])
                leveled_height_data = self._Height_Levelling_3Point(height_data[0])
        else:
            leveled_height_data = self._Height_Levelling_3Point(self.all_data[self.channels.index(mask_channel)])

        mask_array = self._Create_Mask_Array(leveled_height_data, threshold)

        
        if mask_data == False:
            threshold = self._Get_Height_Treshold(leveled_height_data, mask_array, threshold)
            mask_data == True #unnecessary from now on...
        self._Write_to_Logfile('height_masking_threshold', threshold)
        mask_array = self._Create_Mask_Array(leveled_height_data, threshold)
        self.mask_array = mask_array
        if export == True:
            # open files for the masked data:
            for channel in channels:
                header, NUL= self._Create_Header(channel)
                data, trash = self._Load_Data([channel])
                datafile = open("".join([self.directory_name,"/",self.filename," ",channel,"_masked.gsf"]),"bw")
                datafile.write(header)
                datafile.write(NUL)
                masked_array = np.multiply(data[0], mask_array)
                flattened_data = masked_array.flatten()
                
                for pixeldata in flattened_data:
                    datafile.write(pack("f", pixeldata))
                datafile.close()
            print('All channels have been masked and exported!')
        elif export == False:
            dataset = self.all_data
            print('Channels in memory have been masked!')
            for i in range(len(dataset)):
                if self.height_channel not in self.all_data_dict[i]:
                    self.all_data[i] = np.multiply(dataset[i], mask_array)
                self.all_data_dict[i] = self.all_data_dict[i] + '_masked'

    def _Check_Pixel_Position(self, xres, yres, x, y) -> bool:
        '''This function checks if the pixel position is within the bounds'''
        if x < 0 or x > xres:
            return False
        elif y < 0 or y > yres:
            return False
        else: return True

    def _Get_Mean_Value(self, data, x_coord, y_coord, zone) -> float:
        '''This function returns the mean value of the pixel and its nearest neighbors.
        The zone specifies the number of neighbors. 1 means the pixel and the 8 nearest pixels.
        2 means zone 1 plus the next 16, so a total of 25 with the pixel in the middle. 
        '''
        xres = len(data[0])
        yres = len(data)
        size = 2*zone + 1
        mean = 0
        count = 0
        for y in range(size):
            for x in range(size):
                y_pixel = int(y_coord -(size-1)/2 + y)
                x_pixel = int(x_coord -(size-1)/2 + x)
                if self._Check_Pixel_Position(xres, yres, x_pixel, y_pixel) == True:
                    mean += data[y_pixel][x_pixel]
                    count += 1
        return mean/count

    def _Height_Levelling_3Point(self, height_data, zone=1) -> np.array:
        fig, ax = plt.subplots()
        ax.pcolormesh(height_data, cmap='gray')
        klicker = clicker(ax, ["event"], markers=["x"])
        ax.legend()
        ax.axis('scaled')
        plt.title('3 Point leveling: please click on three points\nto specify the underground plane.')
        plt.show()
        klicker_coords = klicker.get_positions()['event'] #klicker returns a dictionary for the events
        klick_coordinates = [[round(element[0]), round(element[1])] for element in klicker_coords]
        self._Write_to_Logfile('height_leveling_coordinates', klick_coordinates)
        if len(klick_coordinates) != 3:
            print('You need to specify 3 point coordinates! \nDo you want to try again?')
            user_input = self._User_Input_Bool()
            if user_input == True:
                self._Height_Levelling_3Point(zone)
            else:
                exit()
        # for the 3 point coordinates the height data is calculated over a small area around the clicked pixels to reduce deviations due to noise
        mean_values = [self._Get_Mean_Value(height_data, klick_coordinates[i][0], klick_coordinates[i][1], zone) for i in range(len(klick_coordinates))]
        matrix = [[klick_coordinates[i][0], klick_coordinates[i][1], mean_values[i]] for i in range(3)]
        A = matrix
        b = [100,100,100] # not sure why, 100 is a bit random, but 0 didn't work
        solution = np.linalg.solve(A, b)
        yres = len(height_data)
        xres = len(height_data[0])
        # create a plane with same dimensions as the height_data
        plane_data = np.zeros((yres, xres))
        for y in range(yres):
            for x in range(xres):
                plane_data[y][x] = -(solution[0]*x + solution[1]*y)/solution[2]
        leveled_height_data = np.zeros((yres, xres))
        # substract the plane_data from the height_data
        for y in range(yres):
            for x in range(xres):
                leveled_height_data[y][x] = height_data[y][x] - plane_data[y][x]
        
        return leveled_height_data

    def Level_Height_Channels(self, channels:list=None) -> None:
        '''This function levels all height channels which are either user specified or in the instance memory.
        The leveling will prompt the user with a preview to select 3 points for getting the coordinates of the leveling plane.
        
        Args:
            channels [list]: list of channels, will override the already existing channels
        '''
        self._Initialize_Data(channels)
        for channel in self.channels:
            if self.height_indicator in channel:
                self.all_data[self.channels.index(channel)] = self._Height_Levelling_3Point(self.all_data[self.channels.index(channel)])
                self.all_data_dict[self.channels.index(channel)] += '_leveled' 
                
    def _Fit_Horizontal_WG(self, data):
        YRes = len(data)
        XRes = len(data[0])
        #just calculate the shift for each pixel for now
        number_align_points = XRes #the number of intersections fitted with gaussian to find waveguide center along the x direction
        align_points = np.arange(0, XRes, int((XRes)/number_align_points), int)
        cutline_data_sets = []
        for element in align_points:
            cutline = []
            for i in range(YRes):
                cutline.append(data[i][element]) # *pow(10, 9) transform height data to nm
            cutline_data_sets.append(cutline)
        list_of_coefficients = []
        p0 = [100, (YRes)/2, 5, 0]
        bounds = ([0, -YRes, 0, -1000], [1000, YRes, YRes/2, 1000])
        for cutline in cutline_data_sets:
            coeff, var_matrix = curve_fit(Gauss_Function, range(0, YRes), cutline, p0=p0, bounds=bounds)
            list_of_coefficients.append(coeff)
            p0 = coeff #set the starting parameters for the next fit
        print("fit succsessful")
        return align_points, list_of_coefficients

    def _Shift_Data(self, data, y_shifts) -> np.array:
        YRes = len(data)
        XRes = len(data[0])
        min_shift = round(min(y_shifts))
        max_shift = round(max(y_shifts))
        new_YRes = YRes + int(abs(min_shift-max_shift))
        data_shifted = np.zeros((new_YRes, XRes))
        #create the realigned height
        for x in range(XRes):
            y_shift = int(-y_shifts[x] + abs(max_shift)) #the calculated shift has to be compensated by shifting the pixels
            for y in range(YRes):
                data_shifted[y + y_shift][x] = data[y][x]
        return data_shifted

    def Realign(self, channels:list=None):
        '''This function corrects the drift of the piezo motor. As of now it needs to be fitted to a region of the sample which is assumed to be straight.
        In the future this could be implemented with a general map containing the distortion created by the piezo motor, if it turns out to be constant...
        Anyways, you will be prompted with a preview of the height data, please select an area of the scan with only one 'straight' waveguide. 
        The bounds for the fitting routine are based on the lower and upper limit of this selection.

        Careful! Will not yet affect the scan size, so the pixelsize will be altered... ToDo
        
        Args:  
            channels [list]: list of channels, will override the already existing channels
        
        '''
        self._Initialize_Data(channels)
        # store the bounds in the instance so the plotting algorithm can access them
        # get the bounds from drawing a rectangle:
        if self.height_channel in self.channels:
            data = self.all_data[self.channels.index(self.height_channel)]
        else:
            data, trash = self._Load_Data([self.height_channel])
        coords = Select_Rectangle(data, self.height_channel)
        lower = coords[0][1]
        upper = coords[1][1]
        self.lower_y_bound = lower
        self.upper_y_bound = upper
        self._Write_to_Logfile('realign_bounds', [lower, upper])
        if self.height_channel in self.channels:
            height_data = self.all_data[self.channels.index(self.height_channel)]
        else:
            height_data_array, trash = self._Load_Data([self.height_channel])
            height_data = height_data_array[0]
            # if the channels have been scaled, the height has to be scaled as well
            scaling = self._Get_Channel_Scaling(0)
            if scaling != 1:
                height_data = self._Scale_Array(height_data, self.height_channel, scaling)
        YRes = len(height_data)
        XRes = len(height_data[0])
        reduced_height_data = np.zeros((upper-lower +1,XRes))
        for y in range(YRes):
            if (lower <= y) and (y <= upper):
                for x in range(XRes):
                    reduced_height_data[y-lower][x] = height_data[y][x]
        align_points, fit_coefficients = self._Fit_Horizontal_WG(reduced_height_data)
        y_shifts = [round(coeff[1],0) -int((upper - lower)/2) for coeff in fit_coefficients]
        # save the align points and y_shifts as instance variables so the plotting algorithm can access them
        self.align_points = align_points
        self.y_shifts = y_shifts

        # plot 
        fig, axis = plt.subplots()    
        fig.set_figheight(self.figsizey)
        fig.set_figwidth(self.figsizex) 
        cmap = 'gray'
        img = axis.pcolormesh(height_data, cmap=cmap)
        # axis.invert_yaxis()
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('height [nm]', rotation=270)
        axis.set_title('realinging')
        axis.axis('scaled')
        axis.plot(self.align_points, [element + lower + (upper-lower)/2 for element in self.y_shifts], color='red')
        axis.hlines([self.upper_y_bound, self.lower_y_bound], xmin=0, xmax=XRes, color='white')
        plt.show()

        # reinitialize the instance data to fit the new bigger arrays
        min_shift = round(min(y_shifts))
        max_shift = round(max(y_shifts))
        new_YRes = YRes + int(abs(min_shift-max_shift))
        all_data = self.all_data
        # self.all_data = np.zeros((len(all_data), new_YRes, XRes))
        self.all_data = []
        for i in range(len(self.channels)):
            shifted_data = self._Shift_Data(all_data[i], y_shifts)
            
            self.all_data.append(shifted_data)
            self.all_data_dict[i] += '_shifted'

    def Cut_Channels(self, channels:list=None, autocut:bool=False) -> None:
        '''This function cuts the specified channels to the specified region.

        Careful! Will not yet affect the scan size, so the pixelsize will be altered... ToDo

        Args:
            channels [list]: list of channels, will override the already existing channels 
            autocut [bool]: if set to 'True' the program will automatically try to remove zero lines and columns, which can result from masking
        '''
        self._Initialize_Data(channels)
        # check if height channel in channels and apply mask to it, until now it has not been masked in order to show the mask in the image
        if (self.height_channel in self.channels) and (len(self.mask_array) > 0):
            self.all_data[self.channels.index(self.height_channel)] = np.multiply(self.all_data[self.channels.index(self.height_channel)], self.mask_array)
        if autocut == True:
            self._Auto_Cut_Channels()
            self._Write_to_Logfile('auto_cut', True)
        else:
            if self.height_channel in self.channels:
                data = self.all_data[self.channels.index(self.height_channel)]
                channel = self.height_channel
            else:
                data = self.all_data[0]
                channel = self.channels[0]
            # get the coordinates of the selection rectangle
            coords = Select_Rectangle(data, channel)
            self._Write_to_Logfile('cut_coords', coords)
            # use the selection to create a mask and multiply to all channels, then apply auto_cut function
            yres = len(data)
            xres = len(data[0])
            mask_array = np.zeros((yres, xres))
            for y in range(yres):
                if y in range(coords[0][1], coords[1][1]):
                    for x in range(xres):
                        if x in range(coords[0][0], coords[1][0]):
                            mask_array[y][x] = 1
            for channel in self.channels:
                self.all_data[self.channels.index(channel)] = np.multiply(self.all_data[self.channels.index(channel)], mask_array)
            self._Auto_Cut_Channels()
        
    def _Auto_Cut_Channels(self) -> None:
        '''This function automatically cuts away all rows and lines which are only filled with zeros.
        This function applies to all channels in memory.
        '''
        # get the new size of the reduced channels
        reduced_data = self._Auto_Cut_Data(self.all_data[0])
        yres = len(reduced_data)
        xres = len(reduced_data[0])
        # copy old data to local variable
        all_data = self.all_data
        # reinitialize self.all_data, all channels must have the same size
        # self.all_data = np.zeros((len(all_data), yres, xres))
        self.all_data = []
        for i in range(len(self.channels)):
            reduced_data = self._Auto_Cut_Data(all_data[i])
            self.all_data.append(reduced_data)
            self.all_data_dict[i] += '_reduced'
 
    def _Auto_Cut_Data(self, data) -> np.array:
        '''This function cuts the data and removes zero values from the outside.'''
        xres = len(data[0])
        yres = len(data)
        # find empty columns and rows to delete:
        columns = []
        for x in range(xres):
            add_to_columns = True
            for y in range(yres):
                if data[y][x] != 0:
                    add_to_columns = False
            if add_to_columns == True:
                columns.append(x)
        rows = []
        for y in range(yres):
            add_to_rows = True
            for x in range(xres):
                if data[y][x] != 0:
                    add_to_rows = False
            if add_to_rows == True:
                rows.append(y)
        
        # create reduced data array
        x_reduced = xres - len(columns)
        y_reduced = yres - len(rows)
        data_reduced = np.zeros((y_reduced, x_reduced))
        # iterate through all pixels and check if they are in rows and columns, then add them to the reduced data array
        count_x = 0
        count_y = 0
        for y in range(yres):
            if y not in rows:
                for x in range(xres):
                    if x not in columns:
                        data_reduced[count_y][count_x] = data[y][x] 
                        count_x += 1
                count_x = 0
                count_y += 1
        return data_reduced

    def Scalebar(self, channels:list=[], units="m", dimension="si-length", label=None, length_fraction=None, height_fraction=None, width_fraction=None,
            location=None, loc=None, pad=None, border_pad=None, sep=None, frameon=None, color=None, box_color=None, box_alpha=None, scale_loc=None,
            label_loc=None, font_properties=None, label_formatter=None, scale_formatter=None, fixed_value=None, fixed_units=None, animated=False, rotation=None):
        '''channels contains the list of channels wich should have a scalebar
        Args:
            channels [list]: list of channels the scalebar should be added to, will not affect the channels in memory
            various definitions for the scalebar, please look up 'matplotlib_scalebar.scalebar' for more information
        '''
        
        # scalebar = ScaleBar(dx, units, dimension, label, length_fraction, height_fraction, width_fraction,
            # location, loc, pad, border_pad, sep, frameon, color, box_color, box_alpha, scale_loc,
            # label_loc, font_properties, label_formatter, scale_formatter, fixed_value, fixed_units, animated, rotation)
        
        
        count = 0
        for channel in self.channels:
            XRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area][0]
            XReal = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.scan_area][0]
            dx = XReal/(XRes*self.scaling_factor)
            scalebar_var = [dx, units, dimension, label, length_fraction, height_fraction, width_fraction,
                            location, loc, pad, border_pad, sep, frameon, color, box_color, box_alpha, scale_loc,
                            label_loc, font_properties, label_formatter, scale_formatter, fixed_value, fixed_units, animated, rotation]
            if (channel in channels) or (len(channels)==0):
                self.scalebar.append([channel, scalebar_var])                
            else:
                self.scalebar.append([channel, None])                
            count += 1

    def Rotate_90_deg(self, orientation:str = 'right'):
        '''This function will rotate all data in memory by 90 degrees.'''
        self._Write_to_Logfile('rotate_90_deg', orientation)
        if orientation == 'right':
            axes=(1,0)
            self._Write_to_Logfile('rotation', +90)
        elif orientation == 'left':
            axes=(0,1)
            self._Write_to_Logfile('rotation', -90)
        #rotate data:
        all_data = self.all_data
        # initialize data array
        # print(self.channels)
        # self.all_data = np.zeros((len(self.channels), self.XRes, self.YRes))
        self.all_data = []
        for channel in self.channels:
            # flip pixelarea and scanarea as well
            XReal, YReal = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.scan_area]
            self.channel_tag_dict[self.channels.index(channel)][Tag_Type.scan_area] = [YReal, XReal]
            XRes, YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area]
            self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area] = [YRes, XRes]
            # data = self.all_data[self.channels.index(channel)]
            self.all_data.append(np.rot90(all_data[self.channels.index(channel)], axes=axes))

    def _Get_Positions_from_Plot(self, data, channel) -> list:
        if self.phase_indicator in channel:
            cmap = SNOM_phase
        elif self.amp_indicator in channel:
            cmap = SNOM_amplitude
        elif self.height_indicator in channel:
            cmap = 'gray'

        fig, ax = plt.subplots()
        img = ax.pcolormesh(data, cmap=cmap)
        klicker = clicker(ax, ["event"], markers=["x"])
        ax.invert_yaxis()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(channel, rotation=270)
        ax.legend()
        ax.axis('scaled')
        plt.title('Please select one or more points to continue.')
        plt.tight_layout()
        plt.show()
        klicker_coords = klicker.get_positions()['event'] #klicker returns a dictionary for the events
        klick_coordinates = [[round(element[0]), round(element[1])] for element in klicker_coords]
        return klick_coordinates

    def _Get_Profile(self, data, coordinates, orientation, width) -> list:
        YRes = len(data)
        XRes = len(data[0])
        all_profiles = []
        for coord in coordinates:
            profile = []
            if orientation == Definitions.vertical:
                for y in range(YRes):
                    # count = 0
                    value = 0
                    for x in range(int(coord[0] - width/2), int(coord[0] + width/2)):
                        value += data[y][x]
                        # count += 1
                    value = value/width
                    profile.append(value)
                    # print('count: ', count)
                    # print('width: ', width)
            if orientation == Definitions.horizontal:
                for x in range(XRes):
                    value = 0
                    for y in range(int(coord[1] - width/2), int(coord[1] + width/2)):
                        value += data[y][x]
                    value = value/width
                    profile.append(value)
            all_profiles.append(profile)
        return all_profiles

    def Select_Profile(self, profile_channel:str, preview_channel:str=None, orientation:Definitions=Definitions.vertical, width=10, phase_orientation=1, coordinates=None):
        '''This function lets the user select a profile with given width in pixels and displays the data.'''
        if preview_channel is None:
            preview_channel = self.height_channel
        if coordinates == None:
            previewdata = self.all_data[self.channels.index(preview_channel)]
            coordinates = self._Get_Positions_from_Plot(previewdata, preview_channel)
            print('The coordinates you selected are:', coordinates)

        profiledata = self.all_data[self.channels.index(profile_channel)]

        cmap = 'gray'
        fig, ax = plt.subplots()
        img = ax.pcolormesh(profiledata, cmap=cmap)
        ax.invert_yaxis()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('height [nm]', rotation=270)
        ax.legend()
        ax.axis('scaled')
        xcoord = [coord[0] for coord in coordinates]
        ycoord = [coord[1] for coord in coordinates]
        if orientation == Definitions.vertical:
            ax.vlines(xcoord, ymin=0, ymax=len(profiledata))
        elif orientation == Definitions.horizontal:
            ax.hlines(ycoord, xmin=0, xmax=len(profiledata[0]))
        plt.title('You chose the following line profiles')
        plt.tight_layout()
        plt.show()
        # it would be nice to be able to add non pcolormesh plots to the subplotslist
        # self.all_subplots.append()
        #get the current values for the scan area and resolution
        XRes, YRes = self.channel_tag_dict[self.channels.index(profile_channel)][Tag_Type.pixel_area]
        XReal, YReal = self.channel_tag_dict[self.channels.index(profile_channel)][Tag_Type.scan_area]
        profiles = self._Get_Profile(profiledata, coordinates, orientation, width)
        for profile in profiles:
            if orientation == Definitions.horizontal:
                values = np.arange(len(profile))*XReal/XRes*pow(10,6)# convert m to µm
            elif orientation == Definitions.vertical:
                values = np.arange(len(profile))*YReal/YRes*pow(10,6)# convert m to µm
            plt.plot(values, profile)
        plt.title('Height profile')
        if orientation == Definitions.vertical:
            plt.xlabel('y-pos [µm]')
        elif orientation == Definitions.horizontal:
            plt.xlabel('x-pos [µm]')
        plt.ylabel('height [nm]')
        plt.tight_layout()
        plt.show()

    def _Scale_Data_XY(self, data, scale_x, scale_y) -> np.array:
        XRes = len(data[0])
        YRes = len(data)
        new_data = np.zeros((YRes*scale_y, XRes*scale_x))
        for y in range(YRes):
            for i in range(scale_y):
                for x in range(XRes):
                    for j in range(scale_x):
                        new_data[y*scale_y + i][x*scale_x + j]= data[y][x]
        return new_data

    def Quadratic_Pixels(self, channels:list=None):
        '''This function scales the data such that each pixel is quadratic, eg. the physical dimensions are equal.
        This is important because the pixels will be set to quadratic in the plotting function.
        
        Args:
            channels [list]: list of channels the scaling should be applied to. If not specified the scaling will be applied to all channels
        '''
        self._Write_to_Logfile('quadratic_pixels', True)
        if channels == None:
            channels = self.channels
        for channel in channels:
            if channel in self.channels:
                XRes, YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area]
                XReal, YReal = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.scan_area]
                pixel_size_x = round(XReal/XRes *1000000000) # pixel size in nm
                pixel_size_y = round(YReal/YRes *1000000000)
                # print('pixelsize x: ', pixel_size_x)
                # print('pixelsize y: ', pixel_size_y)
                scale_x = 1
                scale_y = 1
                if pixel_size_x < pixel_size_y:
                    # print('scale y: ', pixel_size_y/pixel_size_x)
                    scale_y = int(pixel_size_y/pixel_size_x)
                elif pixel_size_x > pixel_size_y:
                    # print('scale x: ', pixel_size_x/pixel_size_y)
                    scale_x = int(pixel_size_x/pixel_size_y)
                # print(pixel_size_x/scale_x, '!=', pixel_size_y/scale_y)
                if pixel_size_x/scale_x != pixel_size_y/scale_y:
                    print('The pixel size does not fit perfectly, you probably chose weired resolution values. You should probably not use this function then...')
                # print('scale x: ', scale_x)
                # print('scale y: ', scale_y)
                self.all_data[self.channels.index(channel)] = self._Scale_Data_XY(self.all_data[self.channels.index(channel)], scale_x, scale_y)
                self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area] = [XRes*scale_x, YRes*scale_y]
        
        

# needed for the Realign function
def Gauss_Function(x, A, mu, sigma, offset):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + offset



