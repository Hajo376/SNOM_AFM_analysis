##########################################################################
# This code was created by Hans-Joachim Schill, University of Bonn, 2022 #
##########################################################################
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
from struct import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_point_clicker import clicker# used for getting coordinates from images
from matplotlib_scalebar.scalebar import ScaleBar # used for creating scale bars
import numpy as np
import pandas as pd # used for getting data out of html files
from datetime import datetime
from enum import Enum, auto
from pathlib import Path, PurePath

# import own functionality
from SNOM_AFM_analysis.lib.snom_colormaps import *
from SNOM_AFM_analysis.lib.phase_slider import Get_Phase_Offset
from SNOM_AFM_analysis.lib.rectangle_selector import Select_Rectangle
from SNOM_AFM_analysis.lib.get_directionality import ChiralCoupler
from SNOM_AFM_analysis.lib import realign
from SNOM_AFM_analysis.lib import profile
from SNOM_AFM_analysis.lib import phase_analysis

class Definitions(Enum):
    vertical = auto()
    horizontal = auto()

class File_Type(Enum):
    """Different file types may vary accross different platforms."""
    standard = auto()
    standard_new = auto()
    aachen_ascii = auto()
    aachen_gsf = auto()
    neaspec_version_1_6_3359_1 = auto()
    # for the parameters:
    html = auto()
    txt = auto() # only for aachen style txt?!
    html_new = auto()# new software version creates slightly different html file
    html_neaspec_version_1_6_3359_1 = auto()
    comsol_gsf = auto()
    comsol_txt = auto()

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
    pixel_scaling = auto()

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
    # Define Plot font sizes
    font_size_default = 8
    font_size_axes_title = 12
    font_size_axes_label = 10
    font_size_tick_labels = 8
    font_size_legend = 8
    font_size_fig_title = 12
    #definitions for color bar ranges:
    height_cbar_range = True
    vmin_height = None
    vmax_height = None
    amp_cbar_range = True
    vmin_amp = None#1 # to make shure that the values will be initialized with the first plotting command
    vmax_amp = None#-1
    # phase_cbar_range = True
    full_phase_range = True # this will overwrite the cbar
    shared_phase_range = True # only used if full phase range is false
    vmin_phase = None
    vmax_phase = None
    real_cbar_range = True
    vlimit_real = None
    # vmin_real = None
    # vmax_real = None
    # show plot automatically? turn to false for gui programming
    show_plot = True
    

    

class Open_Measurement(File_Definitions, Plot_Definitions):
    """This is the main class. You need to specify the measurement folder containing the individual channels.
    You can also specify the channels you want to investigate in the initialisation or at a later point using the 'Initialize_Channels' method.
    You should safe the created measurement instance so you can then call functions on it.
    For example you can call the 'Display_Channels' method, this will simply display all channels currently in memory.
    You can use most methods without specifying channels, then the channels in memory will be used.
    If you specify new channels they will replace the old ones! 
    Only use methods which do not start with an underscore, those are only ment for usage within the class.

    All currently availabe methods:
        Initialize_Channels : Change current channels
        Scale_Channels : Increse pixel amount, important if you want to blurr your data or you want to do some shift.
            This way you can shift with higher precision.
        Set_Min_to_Zero : Sets the minimum hight in the channels containing hight to zero.
        Remove_Subplots : Delete specific (specified by index) subplot from memory.
        Remove_Last_Subplots : Delete the last subplot from memory.
        Switch_Supplots : Rearange subplots by switching positions.
        Display_All_Subplots : Display all suplots which are currently in memory.
        Display_Channels : Display all channels in memory.
        Gauss_Filter_Channels : Simple gaussian blurr, should only be used on amplitude, realpart and height channels.
        Gauss_Filter_Channels_complex : Works additionally for phase, but you need to specify matching amplitude and phase channels.
        Fourier_Filter_Channels : Simple fourier filter of channels, not yet fully implemented.
        Save_to_gsf : Save all channels to .gsf files, dont worry manipulated channels will have a filename extension so you dont overwrite the originals.
        Save_to_txt : Save all channels to .txt files, only use if necessary. gsf is more common and saves space.
        Synccorrection : Used to correct the phasedrift introduced by using the syncronized scan mode, only availabe in the transmission module. 
            You need to specify the used wavelength. All phase channels will be corrected and additionaly also the real parts will be created
        Heigth_Mask_Channels : Simple masking method, based on a height threshold. Works best for leveled data.
        Correct_Phase_Drift : Correct a linear phasedrift on the y-axis by selecting two points which should have the same phase.
        Level_Height_Channels : Simple height leveling based on a simple 3-point correction.
        Shift_Phase : Shift the current phase offset by slider.
        Realign : Used for realigning long measurements of straight wavegides. Relatively specific use case...
        Cut_Channels : Select a rectangle of your measurement to cut out.
            Can be combined with height masking, since if set to 'auto' all black on the outside will be removed automatically.
        Scalebar : Creates a scalebar in the specified channels. Only for plotting ...
        Rotate_90_deg : Rotate measurement by 90 deg.
        Select_Profile : Select a simple profile, for now only horizontal or vertical and over full measurement width or height
        Select_Profiles : Select multiple profiles ...
        Select_Profiles_SSH : specific use case ...
        Display_Profiles : Display all current profiles
        Display_Flattened_Profile : Display all current profiles and flatten them first, only for phase gradients
        Display_Phase_Difference : Display phase difference between phase profiles
        Quadratic_Pixels : This will scale the data automatically if possible, e.g. if you measured with unequal resotution in x and y.
            Only use if resolution difference is an integer multiple like 50 nm in x and 100 nm or 150 nm in y. 
        Overlay_Forward_and_Backward_Channels : You are tired of trowing away all the data in the backwards channels? Then use the overlay function!
            This will try to overlay forward and backward channels and will also try to shift them over each other.

    
    """
    # reorganize and put the following variables in the init function
    # ToDo
    all_subplots = []
    # mask_array = []
    # scalebar = [] # specifiy the channels to add the scalebar to, should contain the channel plus the scalbar object in a nested list
    # piezo motion distortion correction:
    # upper_y_bound = None
    # lower_y_bound = None
    # align_points = None
    # y_shifts = None
    # scaling_factor = 1 # change such that each channel has his own scaling factor

    def __init__(self, directory_name:str, channels:list=None, title:str=None, autoscale:bool=True) -> None:
        """Create a measurement object.

        Args:
            directory_name (str): the directory of the measurement.
            channels (list, optional): list of channels to load. If none are given, a set of default channels are chosen. Defaults to None.
            title (str, optional): title to display when plotting. Defaults to None.
            autoscale (bool, optional): tries to automatically scale the data to have quadratic pixels and undistorted dimensions.
                Only for data where the x and y resolutions are an integer multiple of each other, e.g. xres=50nm and yres=50nm or 100nm or 150nm. Defaults to True.
        """
        # self.directory_name = directory_name
        self.directory_name = Path(directory_name)
        # print('old directory: ', self.directory_name)
        # print('new directory: ', Path(self.directory_name))
        # self.filename = directory_name.split('/')[-1]
        self.filename = Path(PurePath(self.directory_name).parts[-1])
        # print('old filename: ', self.filename)
        # print('new filename: ', PurePath(Path(self.directory_name)).parts[-1])
        self.measurement_title = title # If a measurement_title is specified it will precede the automatically created title based on the channel dictionary
        self.autoscale = autoscale
        self.logfile_path = self._Initialize_Logfile()
        self._Initialize_File_Type()
        if channels == None: # the standard channels which will be used if no channels are specified
            if self.file_type == File_Type.comsol_gsf:
                channels = [self.preview_ampchannel, self.preview_phasechannel]
            else:
                channels = [self.preview_ampchannel, self.preview_phasechannel, self.height_channel]
        self.channels = channels
        self._Create_Measurement_Tag_Dict()
        # just initialize data here?
        # self._Create_Channels_Tag_Dict()
        self._Initialize_Data(self.channels)
        
        # Create a variable containing the data of the specified channels which can be varied later on with e.g. the gaussian_filter and subsequent methods
        # Aditionally a dictionary is created which contains the channel information
        # The dictionary and self.channels are not identical as self.channels only contains the raw channel information, whereas the dictionary is modified to contain
        # the information of which modifications like scaling, filtering ... have been applied.
        # self.all_data, self.channels_label = self._Load_Data(self.channels)
        # print(self.all_data)
        # if autoscale == True:
        #     self.Quadratic_Pixels()

    def _Find_Filetype(self) -> None:
        """This function aims at finding specific characteristics in the filename to idendify the filetype.
        For example the difference in File_Type.standard and File_Type.standard_new are an additional ' raw' at the end of the filename."""
        try:
            # f_1=open(f"{self.directory_name}/{self.filename} O1A.gsf","br")
            f_1=open(self.directory_name / Path(self.filename.name + ' O1A.gsf'),"br")
        except:
            # filetype is at least not standard
            try:
                # f_2=open(f"{self.directory_name}/{self.filename} O1A raw.gsf","br")
                f_2=open(self.directory_name / Path(self.filename.name + ' O1A raw.gsf'),"br")
            except:
                try:
                    # f_3=open(f"{self.directory_name}/{self.filename}_parameters.txt","r")
                    f_3=open(self.directory_name / Path(self.filename.name + '_parameters.txt'),"r")
                except:
                    try:
                        # f_4=open(f"{self.directory_name}/{self.filename}_O1-F-abs.ascii", 'r')
                        f_4=open(self.directory_name / Path(self.filename.name + '_O1-F-abs.ascii'), 'r')
                    except:
                        print("The correct filetype could not automatically be found. Please try again and specifiy the filetype.")
                    else:
                        f_4.close()
                        File_Definitions.file_type = File_Type.aachen_ascii
                        File_Definitions.parmeters_type = File_Type.txt
                else:
                    f_3.close()
                    File_Definitions.file_type = File_Type.comsol_gsf
                    File_Definitions.parmeters_type = File_Type.comsol_txt
            else:
                f_2.close()
                File_Definitions.file_type = File_Type.standard_new
                File_Definitions.parmeters_type = File_Type.html_new
        else:
            f_1.close()
            File_Definitions.file_type = File_Type.standard
            File_Definitions.parmeters_type = File_Type.html
        #alternative way: get the software version from the last entry in the .txt file
        try:
            # parameters = self.directory_name + '/' + self.filename + '.txt' #standard snom files, not for aachen files
            parameters = self.directory_name / Path(self.filename.name + '.txt') #standard snom files, not for aachen files
            file = open(parameters, 'r', encoding="utf-8")
        except:
            pass
        else:
            parameter_list = file.read()
            file.close()
            parameter_list = parameter_list.split('\n')
            # print(parameter_list)
            # print(parameter_list[-2])
            version_number = parameter_list[-2].split('\t')[-1]
            print('version_number:      ', version_number)
            if version_number != ['']:
                self.version_number = version_number
            if version_number == '1.6.3359.1':
                File_Definitions.file_type = File_Type.neaspec_version_1_6_3359_1
                File_Definitions.parmeters_type = File_Type.html_neaspec_version_1_6_3359_1

    def _Initialize_File_Type(self) -> None:
        self._Find_Filetype() # try to find the filetype automatically
        self.file_type = File_Definitions.file_type
        self.parameters_type = File_Definitions.parmeters_type
        if self.file_type == File_Type.standard or self.file_type == File_Type.standard_new or self.file_type == File_Type.neaspec_version_1_6_3359_1:
            self.phase_channels = ['O1P','O2P','O3P','O4P','O5P', 'R-O1P','R-O2P','R-O3P','R-O4P','R-O5P']
            self.amp_channels = ['O1A','O2A','O3A','O4A','O5A', 'R-O1A','R-O2A','R-O3A','R-O4A','R-O5A']
            self.real_channels = ['O1Re', 'O2Re', 'O3Re', 'O4Re', 'R-O5Re', 'R-O1Re', 'R-O2Re', 'R-O3Re', 'R-O4Re', 'R-O5Re']
            self.imag_channels = ['O1Im', 'O2Im', 'O3Im', 'O4Im', 'R-O5Im', 'R-O1Im', 'R-O2Im', 'R-O3Im', 'R-O4Im', 'R-O5Im']
            self.height_channel = 'Z C'
            self.height_channels = ['Z C', 'R-Z C']
            # self.all_channels = ['O1A','O1P','O2A','O2P','O3A','O3P','O4A','O4P','O5A','O5P','R-O1A','R-O1P','R-O2A','R-O2P','R-O3A','R-O3P','R-O4A','R-O4P','R-O5A','R-O5P']
            self.all_channels = self.phase_channels + self.amp_channels + self.height_channels
            self.preview_ampchannel = 'O2A'
            self.preview_phasechannel = 'O2P'
            self.height_indicator = 'Z'
            self.amp_indicator = 'A'
            self.phase_indicator = 'P'
            self.backwards_indicator = 'R-'
            self.real_indicator = 'Re'
            self.imag_indicator = 'Im'
        elif self.file_type == File_Type.aachen_ascii or self.file_type == File_Type.aachen_gsf:
            self.phase_channels = ['O1-F-arg','O2-F-arg','O3-F-arg','O4-F-arg', 'O1-B-arg','O2-B-arg','O3-B-arg','O4-B-arg']
            self.amp_channels = ['O1-F-abs','O2-F-abs','O3-F-abs','O4-F-abs', 'O1-B-abs','O2-B-abs','O3-B-abs','O4-B-abs']
            self.real_channels = ['O1-F-Re','O2-F-Re','O3-F-Re','O4-F-Re','O1-B-Re','O2-B-Re','O3-B-Re','O4-B-Re']
            self.imag_channels = ['O1-F-Im','O2-F-Im','O3-F-Im','O4-F-Im','O1-B-Im','O2-B-Im','O3-B-Im','O4-B-Im']
            self.height_channel = 'MT-F-abs'
            self.height_channels = ['MT-F-abs', 'MT-B-abs']
            # self.all_channels = ['O1-F-abs','O1-F-arg','O2-F-abs','O2-F-arg','O3-F-abs','O3-F-arg','O4-F-abs','O4-F-arg', 'O1-B-abs','O1-B-arg','O2-B-abs','O2-B-arg','O3-B-abs','O3-B-arg','O4-B-abs','O4-B-arg']
            self.all_channels = self.phase_channels + self.amp_channels + self.height_channels
            self.preview_ampchannel = 'O2-F-abs'
            self.preview_phasechannel = 'O2-F-arg'
            self.height_indicator = 'MT'
            self.amp_indicator = 'abs'
            self.phase_indicator = 'arg'
            self.real_indicator = 'Re'#not used
            self.imag_indicator = 'Im'#not used
            self.backwards_indicator = '-B-'
        elif self.file_type == File_Type.comsol_gsf:
            self.all_channels = ['abs', 'arg', 'real']
            self.phase_channels = ['arg']
            self.amp_channels = ['abs']
            self.real_channels = ['real']
            self.height_channel = 'none'
            self.preview_ampchannel = 'abs'
            self.preview_phasechannel = 'arg'
            self.height_indicator = 'none'
            self.amp_indicator = 'abs'
            self.phase_indicator = 'arg'
            self.real_indicator = 'real'
            self.imag_indicator = 'imag'
        
        # additional definitions independent of filetype:
        self.filter_gauss_indicator = 'gauss'
        self.filter_fourier_indicator = 'fft'


        #create also lists for the overlain channels
        self.overlain_phase_channels = [channel+'_overlain_manipulated' for channel in self.phase_channels]
        self.overlayed_amp_channels = [channel+'_overlain_manipulated' for channel in self.amp_channels]
        # self.corrected_phase_channels = ['O1P_corrected','O2P_corrected','O3P_corrected','O4P_corrected','O5P_corrected', 'R-O1P_corrected','R-O2P_corrected','R-O3P_corrected','R-O4P_corrected','R-O5P_corrected']
        self.corrected_phase_channels = [channel+'_corrected' for channel in self.phase_channels]
        self.corrected_overlain_phase_channels = [channel+'_overlain_manipulated' for channel in self.corrected_phase_channels]

    def _Create_Measurement_Tag_Dict(self):
        # create tag_dict for each channel individually? if manipulated channels are loaded they might have different diffrent resolution
        # only center_pos, scan_area, pixel_area and rotation must be stored for each channel individually but rotation is not stored in the original .gsf files
        # but rotation could be added in the newly created .gsf files
        print(f'self.parameters_type: {self.parameters_type}')
        print(f'self.file_type:       {self.file_type}')
        if self.parameters_type == File_Type.html:
            # all_tables = pd.read_html("".join([self.directory_name,"/",self.filename,".html"]))
            all_tables = pd.read_html(self.directory_name / Path(self.filename.name + ".html"))
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
            # all_tables = pd.read_html("".join([self.directory_name,"/",self.filename,".html"]))
            all_tables = pd.read_html(self.directory_name / Path(self.filename.name + ".html"))
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
        elif self.parameters_type == File_Type.html_neaspec_version_1_6_3359_1:
            # all_tables = pd.read_html("".join([self.directory_name,"/",self.filename,".html"]))
            all_tables = pd.read_html(self.directory_name / Path(self.filename.name + ".html"))
            tables = all_tables[0]
            self.measurement_tag_dict = {
                Tag_Type.center_pos: [float(tables[2][3]), float(tables[3][3])],
                Tag_Type.rotation: int(tables[2][4]),
                Tag_Type.scan_area: [float(tables[2][5]), float(tables[3][5])],
                Tag_Type.pixel_area: [int(tables[2][6]), int(tables[3][6])],
                Tag_Type.integration_time: float(tables[2][8]),
                Tag_Type.tip_frequency: float(tables[2][12]),
                Tag_Type.tip_amplitude: float(tables[2][13]),
                Tag_Type.tapping_amplitude: float(tables[2][14])
            }
        elif self.parameters_type == File_Type.txt:

            # parameters = self.directory_name + '/' + self.filename + '.parameters.txt'
            parameters = self.directory_name / Path(self.filename.name + '.parameters.txt')
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
        elif self.parameters_type == File_Type.comsol_txt:
            # parameters = self.directory_name + '/' + self.filename + '_parameters.txt'
            parameters = self.directory_name / Path(self.filename.name + '_parameters.txt')
            file = open(parameters, 'r')
            parameter_list = file.read()
            file.close()
            # print(parameter_list)
            parameter_list = parameter_list.split('\n')
            parameter_list = [element.split('=') for element in parameter_list]
            # center_pos = [float(parameter_list[7][1]), float(parameter_list[8][1])]
            # rotation = float(parameter_list[9][1])
            scan_area = [float(parameter_list[2][1]), float(parameter_list[3][1])]
            pixel_area = [int(parameter_list[0][1]), int(parameter_list[1][1])]
            # integration_time = float(parameter_list[6][1])
            # tip_frequency = float(parameter_list[10][1])
            self.measurement_tag_dict = {
                Tag_Type.scan_type: None,
                Tag_Type.center_pos: None,
                Tag_Type.rotation: None,
                Tag_Type.scan_area: scan_area,
                Tag_Type.pixel_area: pixel_area,
                Tag_Type.integration_time: None,
                Tag_Type.tip_frequency: None,
                Tag_Type.tip_amplitude: None,
                Tag_Type.tapping_amplitude: None
            }
        # only used by synccorrection, every other function should use the channels tag dict version, as pixel resolution could vary
        self.XRes, self.YRes = self.measurement_tag_dict[Tag_Type.pixel_area]
        self.XReal, self.YReal = self.measurement_tag_dict[Tag_Type.scan_area] # in µm

    def _Create_Channels_Tag_Dict(self, channels:list=None):
        # ToDo optimize everything so new filetypes dont need so much extra copies
        if channels == None:
            channels = self.channels
        if (self.file_type == File_Type.standard) or (self.file_type == File_Type.standard_new) or (self.file_type == File_Type.aachen_gsf) or (self.file_type == File_Type.comsol_gsf) or (self.file_type == File_Type.neaspec_version_1_6_3359_1):
            cod="latin1"
            # get the tag values from each .gsf file individually
            if channels == self.channels:
                self.channel_tag_dict = []
            for channel in channels:
                if (self.file_type == File_Type.standard_new or self.file_type==File_Type.neaspec_version_1_6_3359_1) and '_corrected' not in channel:
                    # if ' C' in channel or '_manipulated' in channel: #channel == 'Z C' or channel == 'R-Z C':
                    #     filepath = self.directory_name / Path(self.filename.name + ' ' + channel + '.gsf')
                    # else:
                    #     filepath = self.directory_name / Path(self.filename.name + ' ' + channel + ' raw.gsf')
                    if channel in self.all_channels:
                        filepath = self.directory_name / Path(self.filename.name + ' ' + channel + ' raw.gsf')
                    else:
                        filepath = self.directory_name / Path(self.filename.name + ' ' + channel + '.gsf')

                # elif self.file_type == File_Type.aachen_ascii:
                #     if '_corrected' in channel or '_manipulated' in channel: 
                #         filepath = self.directory_name + '/' + self.filename + ' ' + channel + '.gsf'
                #     else:
                #         filepath = self.directory_name + '/' + self.filename + '_' + channel + '.ascii'
                #         cod = None
                elif self.file_type == File_Type.comsol_gsf:
                    filepath = self.directory_name / Path(self.filename.name + '_' + channel + '.gsf')

                else:
                    # filepath = self.directory_name + '/' + self.filename + ' ' + channel + '.gsf'
                    filepath = self.directory_name / Path(self.filename.name + ' ' + channel + '.gsf')
                print(self.file_type, filepath, self.channels)
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
                    Tag_Type.scan_area: [float(XReal), float(YReal)],
                    Tag_Type.pixel_scaling: 1
                }
                self.channel_tag_dict.append(channel_dict)
            pass
        else:
            if self.parameters_type == File_Type.txt: #only for aachen files
                if channels == self.channels:
                    self.channel_tag_dict = []
                # parameters = self.directory_name + '/' + self.filename + '.parameters.txt'
                parameters = self.directory_name / Path(self.filename.name + '.parameters.txt')
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
                        Tag_Type.scan_area: scan_area,
                        Tag_Type.pixel_scaling: 1
                    }
                # for this file type all channels must be of same size
                for channel in channels:
                    self.channel_tag_dict.append(channel_dict)
            else:
                print('channel tag dict for this filetype is not yet implemented')
    '''
    # old version picky on the ' ' characters
    def _Get_Tagval(self, content, tag):
        """This function gets the value of the tag listed in the file header"""
        taglength=len(tag)
        tagdatastart=content.find(tag)+taglength+1   
        tagdatalength=content[tagdatastart:].find("\n") 
        tagdataend=tagdatastart+tagdatalength
        tagval=float(content[tagdatastart:tagdataend])
        return tagval
    '''
    def _Get_Tagval(self, content, tag):
        """This function gets the value of the tag listed in the file header"""
        # print('trying to split the content')
        # print(content)
        content_array = content.split('\n')
        # print(content_array[0:5])
        tag_array = []
        tagval = 0# if no tag val can be found return 0
        for element in content_array:
            if len(element) > 50: # its probably not part of the header anymore...
                break
            elif '=' not in element:
                pass
            else:
                tag_pair = element.split('=')
                # print(f'tag_pair = {tag_pair}')
                tag_name = tag_pair[0].replace(' ', '')# remove possible ' ' characters
                tag_val = tag_pair[1].replace(' ', '')# remove possible ' ' characters
                tag_array.append([tag_name, tag_val])
        for i in range(len(tag_array)):
            if tag_array[i][0] == tag:
                tagval = float(tag_array[i][1])
        return tagval
    
    def _Initialize_Data(self, channels=None) -> None:
        """This function initializes the data in memory. If no channels are specified the already existing data is used,
        which is created automatically in the instance init method. If channels are specified, the instance data is overwritten.
        Channels must be specified as a list of channels."""
        # print(f'initialising channels: {channels}')
        if channels == None:
            #none means the channels specified in the instance creation should be used
            pass
        else:
            self.channels = channels
            # update the channel tag dictionary, makes the program compatible with differrently sized datasets, like original data plus manipulated, eg. cut data
            self._Create_Channels_Tag_Dict()
            self.all_data, self.channels_label = self._Load_Data(channels)
            xres = len(self.all_data[0][0])
            yres = len(self.all_data[0])
            # reset all the instance variables dependent on the data, but nor the ones responsible for plotting
            # self.scaling_factor = 1
            if self.autoscale == True:
                self.Quadratic_Pixels()
            # initialize instance variables:
            self.mask_array = [] # not shure if it's best to reset the mask...
            self.upper_y_bound = None
            self.lower_y_bound = None
            self.align_points = None
            self.y_shifts = None
            self.scalebar = []    

    def Initialize_Channels(self, channels:list) -> None:
        """This function will load the data from the specified channels and replace the ones in memory.
        
        Args:
            channels [list]: a list containing the channels you want to initialize
        """
        self._Initialize_Data(channels)

    def Add_Channels(self, channels:list) -> None:
        """This function will add the specified channels to memory without changing the already existing ones.

        Args:
            channels (list): Channels to add to memory.
        """
        self.channels += channels
        # update the channel tag dictionary, makes the program compatible with differrently sized datasets, like original data plus manipulated, eg. cut data
        self._Create_Channels_Tag_Dict(channels)
        all_data, channels_label = self._Load_Data(channels)
        for i in range(len(channels)):
            self.all_data.append(all_data[i])
            self.channels_label.append(channels_label[i])
        # reset all the instance variables dependent on the data, but nor the ones responsible for plotting
        # self.scaling_factor = 1
        if self.autoscale == True:
            self.Quadratic_Pixels(channels)

    def _Initialize_Logfile(self) -> str:
        # logfile_path = self.directory_name + '/python_manipulation_log.txt'
        logfile_path = self.directory_name / Path('python_manipulation_log.txt')
        file = open(logfile_path, 'a') # the new logdata will be appended to the existing file
        now = datetime.now()
        current_datetime = now.strftime("%d/%m/%Y %H:%M:%S")
        file.write(current_datetime + '\n' + 'filename = ' + self.filename.name + '\n')
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
        # print('self.XRes: ', self.XRes)
        # print('self.YRes: ', self.YRes)
        # print('self.channels: ', self.channels)
        # print('self.all_subplots[-1]: ', [element[3] for element in self.all_subplots])

    def _Scale_Array(self, array, scaling) -> np.array:
        """This function scales a given 2D Array, it thus creates 'scaling'**2 subpixels per pixel.
        The scaled array is returned."""
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
        """This function scales all the data in memory or the specified channels.
                
        Args:
            channels (list, optional): List of channels to scale. If not specified all channels in memory will be scaled. Defaults to None.
            scaling (int, optional): Defines scaling factor. Each pixel will be scaled to scaling**2 subpixels. Defaults to 4.
        """
        # ToDo: reimplement scaling dependent on axis, x and y independently
        # self._Initialize_Data(channels)
        if channels is None:
            channels = self.channels
        self._Write_to_Logfile('scaling', scaling)
        # self.scaling_factor = scaling
        # dataset = self.all_data
        # channel_tag_dict = self.channel_tag_dict
        # yres = len(dataset[0])
        # xres = len(dataset[0][0])
        # self.all_data = np.zeros((len(dataset), yres*scaling, xres*scaling))
        # re initialize data storage and channel_tag_dict, since resolution is changed
        # self.all_data = []
        # self.channel_tag_dict = []
        for channel in channels:
            if channel in self.channels:
                self.all_data[self.channels.index(channel)] = self._Scale_Array(self.all_data[self.channels.index(channel)], scaling)
                XReal, YReal = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area]
                self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area] = [XReal*scaling, YReal*scaling]
                self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_scaling] = scaling
            else:
                print(f'Channel {channel} is not in memory! Please initiate the channels you want to use first!')

    def _Load_Data(self, channels:list) -> list:
        """Loads all binary data of the specified channels and returns them in a list plus the dictionary with the channel information.
        Height data is automatically converted to nm. """
        if self.file_type == File_Type.comsol_gsf:
            return self._Load_Data_comsol(channels)
        # print('load data for self.channels: ', self.channels)
        # datasize=int(self.XRes*self.YRes*4)
        #create a list containing all the lists of the individual channels
        all_binary_data = []
        #safe the information about which channel is which list in a dictionary
        data_dict = []
        # data_dict = {}
        # all_data = np.zeros((len(channels), self.YRes, self.XRes))
        all_data = []
        if self.file_type==File_Type.standard or self.file_type==File_Type.standard_new or self.file_type ==File_Type.neaspec_version_1_6_3359_1:
            for i in range(len(channels)):
                # print(channels[i])
                if (self.file_type==File_Type.standard_new or self.file_type==File_Type.neaspec_version_1_6_3359_1) and '_corrected' not in channels[i]:
                    # if ' C' in channels[i] or '_manipulated' in channels[i]: #channels[i] == 'Z C' or channels[i] == 'R-Z C':
                    #     f=open(self.directory_name / Path(self.filename.name + f' {channels[i]}.gsf'),"br")
                    # else:
                    #     f=open(self.directory_name / Path(self.filename.name + f' {channels[i]} raw.gsf'),"br")
                    if channels[i] in self.all_channels:
                        f=open(self.directory_name / Path(self.filename.name + f' {channels[i]} raw.gsf'),"br")
                    else:
                        f=open(self.directory_name / Path(self.filename.name + f' {channels[i]}.gsf'),"br")

                else:
                    f=open(self.directory_name / Path(self.filename.name + f' {channels[i]}.gsf'),"br")
                binarydata=f.read()
                f.close()
                all_binary_data.append(binarydata)
                data_dict.append(channels[i])
            count = 0
            for channel in channels:
                # if _manipulated in channel name use channel dict, because resolution etc could be different to original data
                XRes, YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area]
                '''if '_manipulated' in channel:
                    XRes, YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area]
                # dont remember why i added the following but it leads to problems
                # if channel in self.channels:
                #     XRes, YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area]
                else:
                    XRes, YRes = self.measurement_tag_dict[Tag_Type.pixel_area]'''
                
                datasize=int(XRes*YRes*4)
                channel_data = np.zeros((YRes, XRes))
                reduced_binarydata=all_binary_data[count][-datasize:]
                phaseoffset = 0
                rounding_decimal = 2
                scaling = 1
                if self.amp_indicator in channel:
                    rounding_decimal = 5
                if self.height_indicator in channel:
                    # print('channel: ', channel)
                    scaling = 1000000000#convert to nm
                if self.phase_indicator in channel:
                    # normal phase data ranges from -pi to pi and gets shifted by +pi
                    phaseoffset = np.pi
                    if '_corrected' in channel:
                        # if the data is from a corrected channel it is already shifted
                        phaseoffset = 0
                if self.real_indicator in channel or self.imag_indicator in channel:
                    rounding_decimal = 4
                for y in range(0,YRes):
                    for x in range(0,XRes):
                        pixval=unpack("f",reduced_binarydata[4*(y*XRes+x):4*(y*XRes+x+1)])[0]
                        channel_data[y][x] = round(pixval*scaling + phaseoffset, rounding_decimal)
                all_data.append(channel_data)
                count+=1

        elif self.file_type == File_Type.aachen_ascii:
            for i in range(len(channels)):
                # file = open(f"{self.directory_name}/{self.filename}_{channels[i]}.ascii", 'r')
                file = open(self.directory_name / Path(self.filename.name + f'_{channels[i]}.ascii'), 'r')
                string_data = file.read()
                datalist = string_data.split('\n')
                datalist = [element.split(' ') for element in datalist]
                datalist = np.array(datalist[:-1], dtype=float)#, dtype=np.float convert list to np.array and strings to float
                channel = channels[i]
                phaseoffset = 0
                rounding_decimal = 2
                scaling = 1
                if (self.amp_indicator in channel) and (self.height_indicator not in channel):
                    rounding_decimal = 5
                if self.phase_indicator in channel:
                    phaseoffset = np.pi
                    flattened_data = datalist.flatten()# ToDo just for now, in future the voltages have to be converted
                    phaseoffset = min(flattened_data)
                    if '_corrected' in channel:
                        # if the data is from a corrected channel it is already shifted
                        phaseoffset = 0
                if self.height_indicator in channel:
                    scaling = pow(10, 9)
                XRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area][0]
                YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area][1]
                for y in range(0,YRes):
                    for x in range(0,XRes):
                        datalist[y][x] = round(datalist[y][x]*scaling + phaseoffset, rounding_decimal)
                # all_data[i] = datalist #old
                all_data.append(datalist)
                data_dict.append(channels[i])

        elif self.file_type == File_Type.aachen_gsf:
            for i in range(len(channels)):
                # f=open(f"{self.directory_name}/{self.filename}_{channels[i]}.gsf","br")
                f=open(self.directory_name / Path(self.filename.name + f'_{channels[i]}.gsf'),"br")
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
                if (self.amp_indicator in channel) and (self.height_indicator not in channel):
                    rounding_decimal = 5
                    scaling = pow(10, 6)
                if self.phase_indicator in channel:
                    phaseoffset = np.pi
                    
                    if '_corrected' in channel:
                        # if the data is from a corrected channel it is already shifted
                        phaseoffset = 0
                if self.height_indicator in channel:

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

    def _Load_Data_comsol(self, channels):
        # datasize=int(self.XRes*self.YRes*4)
        #create a list containing all the lists of the individual channels
        all_binary_data = []
        #safe the information about which channel is which list in a dictionary
        data_dict = []
        # data_dict = {}
        # all_data = np.zeros((len(channels), self.YRes, self.XRes))
        all_data = []
        for i in range(len(channels)):
            # print(channels[i])
            # f=open(f"{self.directory_name}/{self.filename} {channels[i]}.gsf","br")
            f=open(self.directory_name / Path(self.filename.name + f'_{channels[i]}.gsf'),"br")
            binarydata=f.read()
            f.close()
            all_binary_data.append(binarydata)
            data_dict.append(channels[i])
        count = 0
        for channel in channels:
            XRes, YRes = self.measurement_tag_dict[Tag_Type.pixel_area]
            
            datasize=int(XRes*YRes*4)
            channel_data = np.zeros((YRes, XRes))
            reduced_binarydata=all_binary_data[count][-datasize:]
            phaseoffset = 0
            rounding_decimal = 2
            scaling = 1
            if self.amp_indicator in channel:
                rounding_decimal = 5
            if self.phase_indicator in channel:
                # normal phase data ranges from -pi to pi and gets shifted by +pi
                phaseoffset = np.pi
            if self.real_indicator in channel or self.imag_indicator in channel:
                rounding_decimal = 4
            for y in range(0,YRes):
                for x in range(0,XRes):
                    pixval=unpack("f",reduced_binarydata[4*(y*XRes+x):4*(y*XRes+x+1)])[0]
                    channel_data[y][x] = round(pixval*scaling + phaseoffset, rounding_decimal)
            all_data.append(channel_data)
            count+=1
        return all_data, data_dict

    def _Load_Data_Binary(self, channels) -> list:
        """Loads all binary data of the specified channels and returns them in a list plus the dictionary for access"""
        #create a list containing all the lists of the individual channels
        all_binary_data = []
        #safe the information about which channel is which list in a dictionary
        data_dict = []
        for i in range(len(channels)):
            # f=open(f"{self.directory_name}/{self.filename} {channels[i]}.gsf","br")
            f=open(self.directory_name / Path(self.filename.name + f' {channels[i]}.gsf'),"br")
            binarydata=f.read()
            f.close()
            all_binary_data.append(binarydata)
            data_dict.append(channels[i])
        return all_binary_data, data_dict

    def Set_Min_to_Zero(self, channels:list=None) -> None:
        """This function sets the min value of the specified channels to zero.
                
        Args:
            channels (list, optional): List of channels to set min value to zero. If not specified this will apply to all height channels in memory. Defaults to None.
        """
        if channels is None:
            channels = []
            for channel in self.channels:
                if self.height_indicator in channel:
                    channels.append(channel)

        self._Write_to_Logfile('set_min_to_zero', True)
        for channel in channels:
            if channel in self.channels:
                data = self.all_data[self.channels.index(channel)]
                flattened_data = data.flatten()
                data_min = min(flattened_data)
                self.all_data[self.channels.index(channel)] = data - data_min
            else:
                print('At least one of the specified channels is not in memory! You probably should initialize the channels first.')

    def _Add_Subplot(self, data, channel, scalebar=None) -> list:
        """This function adds the specified data to the list of subplots. The list of subplots contains the data, the colormap,
        the colormap label and a title, which are generated from the channel information. The same array is also returned,
        so it can also be iterated by an other function to only plot the data of interest."""
        if self.amp_indicator in channel and self.height_indicator not in channel:
            cmap=SNOM_amplitude
            label = 'Amplitude [a.u.]'
            title = f'Amplitude {channel}'
        elif self.phase_indicator in channel:
            if 'positive' in channel:
                cmap = SNOM_phase
                title = f'Positively corrected phase O{channel[1]}P'
            elif 'negative' in channel:
                cmap = SNOM_phase
                title = f'Negatively corrected phase O{channel[1]}P'
            else:
                cmap=SNOM_phase
                title = f'Phase {channel}'
            label = 'Phase'
        elif self.height_indicator in channel:
            cmap=SNOM_height
            label = 'Height [nm]'
            title = f'Height {channel}'
        elif self.real_indicator in channel or self.imag_indicator in channel:
            cmap=SNOM_realpart
            label = 'E [a.u.]'
            if self.real_indicator in channel:
                title = f'Real part {channel}'
            else:
                title = f'Imaginary part {channel}'
        
        
        
        elif self.filter_fourier_indicator in channel:
            cmap='viridis'
            label = 'Intensity [a.u.]'
            title =  f'Fourier Transform {channel}'
        elif self.filter_gauss_indicator in channel:
            title = f'Gauss blurred {channel}'
        else:
            print('In _Add_Subplot(), encountered unknown channel')
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
        
        #old:
        '''
        if self.file_type == File_Type.standard or self.file_type == File_Type.standard_new or self.file_type==File_Type.neaspec_version_1_6_3359_1:
            if self.amp_indicator in channel:
                cmap=SNOM_amplitude
                label = 'Amplitude [a.u.]'
                title = f'Amplitude {channel}'
            elif self.phase_indicator in channel:
                if 'positive' in channel:
                    cmap = SNOM_phase
                    title = f'Positively corrected phase O{channel[1]}P'
                elif 'negative' in channel:
                    cmap = SNOM_phase
                    title = f'Negatively corrected phase O{channel[1]}P'
                else:
                    cmap=SNOM_phase
                    title = f'Phase {channel}'
                label = 'Phase'
            elif self.height_indicator in channel:
                cmap=SNOM_height
                label = 'Height [nm]'
                title = f'Height {channel}'
            elif self.real_indicator in channel or self.imag_indicator in channel:
                cmap=SNOM_realpart
                label = 'E [a.u.]'
                if self.real_indicator in channel:
                    title = f'Real part {channel}'
                else:
                    title = f'Imaginary part {channel}'

        elif self.file_type == File_Type.aachen_ascii or self.file_type == File_Type.aachen_gsf:
            if 'abs' in channel and not 'MT' in channel:
                cmap=SNOM_amplitude
                label = 'Amplitude [a.u.]'
                title = f'Amplitude {channel}'
            elif 'arg' in channel:
                if 'positive' in channel:
                    cmap = SNOM_phase
                    title = f'Positively corrected phase O{channel[1]}P' # ToDo
                elif 'negative' in channel:
                    cmap = SNOM_phase
                    title = f'Negatively corrected phase O{channel[1]}P'
                else:
                    cmap=SNOM_phase
                    title = f'Phase {channel}'
                label = 'phase'
            elif 'MT' in channel:
                cmap=SNOM_height
                label = 'Height [nm]'
                title = f'Height {channel}'
        elif self.file_type == File_Type.comsol_gsf:
            if 'abs' in channel:
                cmap=SNOM_amplitude
                label = 'Amplitude [a.u.]'
                title = f'Amplitude {channel}'
            elif 'arg' in channel:
                cmap=SNOM_phase
                title = f'Phase {channel}'
                label = 'phase'
            elif 'real' in channel:
                cmap=SNOM_realpart
                label = 'E [a.u.]'
                title = f'Realpart {channel}'
        '''
    
    def Remove_Subplots(self, index_array:list) -> None:
        """This function removes the specified subplot from the memory.
        
        Args:
            index_array [list]: The indices of the subplots to remove from the plot list
        """
        #sort the index array in descending order and delete the corresponding plots from the memory
        index_array.sort(reverse=True)
        for index in index_array:
            del self.all_subplots[index]

    def Remove_Last_Subplots(self, times:int=1) -> None:
        """This function removes the last added subplots from the memory.
        Times specifies how often the last subplot should be removed.
        Times=1 means only the last, times=2 means the two last, ...
        
        Args:
            times [int]: how many subplots should be removed from the end of the list?
        """
        for i in range(times):
            self.all_subplots.pop()

    def _Plot_Subplots(self, subplots) -> None:
        """This function plots the subplots. The plots are created in a grid, by default the grid is optimized for 3 by 3. The layout changes dependent on the number of subplots
        of subplots and also the dimensions. Wider subplots are prefferably created vertically, otherwise they are plotted horizontally. Probably subject to future changes..."""
        number_of_axis = 9
        number_of_subplots = len(subplots)
        # print('Number of subplots: {}'.format(number_of_subplots))
        #specify the way the subplots are organized
        nrows = int((number_of_subplots-1)/np.sqrt(number_of_axis))+1
        # set the plotting font sizes:
        plt.rc('font', size=self.font_size_default)          # controls default text sizes
        plt.rc('axes', titlesize=self.font_size_axes_title)     # fontsize of the axes title
        plt.rc('axes', labelsize=self.font_size_axes_label)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=self.font_size_tick_labels)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=self.font_size_tick_labels)    # fontsize of the tick labels
        plt.rc('legend', fontsize=self.font_size_legend)    # legend fontsize
        plt.rc('figure', titlesize=self.font_size_fig_title)  # fontsize of the figure title

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

        #start the plotting process
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
                        # print('added a scalebar')
                    #center the colorscale for real data around 0
                    # get minima and maxima from data:
                    flattened_data = data.flatten()
                    min_data = min(flattened_data)
                    max_data = max(flattened_data)
                    # print('min: ', min_data)
                    # print('max: ', max_data)

                    if self.real_indicator in title or self.imag_indicator in title: # for real part or imaginary part data
                        if self.file_type == File_Type.comsol_gsf:
                            data = Set_nan_to_zero(data) #comsol data can contain nan values which are problematic for min and max
                        data_limit = Get_Largest_Abs(min_data, max_data)
                        if Plot_Definitions.vlimit_real is None: Plot_Definitions.vlimit_real = data_limit
                        
                        

                        
                        if Plot_Definitions.real_cbar_range is True:
                            if Plot_Definitions.vlimit_real < data_limit: Plot_Definitions.vlimit_real = data_limit
                            img = axis.pcolormesh(data, cmap=cmap, vmin=-Plot_Definitions.vlimit_real, vmax=Plot_Definitions.vlimit_real)
                        else:
                            img = axis.pcolormesh(data, cmap=cmap, vmin=-data_limit, vmax=data_limit)
                    else:
                        if cmap == SNOM_phase and Plot_Definitions.full_phase_range is True: # for phase data
                            # print('plotting full range phase')
                            vmin = 0
                            vmax = 2*np.pi
                            img = axis.pcolormesh(data, cmap=cmap, vmin=vmin, vmax=vmax)
                        elif cmap == SNOM_phase and Plot_Definitions.full_phase_range is False:
                            if Plot_Definitions.vmin_phase is None: Plot_Definitions.vmin_phase = min_data
                            if Plot_Definitions.vmax_phase is None: Plot_Definitions.vmax_phase = max_data
                            if Plot_Definitions.shared_phase_range is True:
                                if Plot_Definitions.vmin_phase > min_data: Plot_Definitions.vmin_phase = min_data
                                if Plot_Definitions.vmax_phase < max_data: Plot_Definitions.vmax_phase = max_data
                            else:
                                Plot_Definitions.vmin_phase = min_data
                                Plot_Definitions.vmax_phase = max_data
                                # vmin = min_data
                                # vmax = max_data
                            img = axis.pcolormesh(data, cmap=cmap, vmin=Plot_Definitions.vmin_phase, vmax=Plot_Definitions.vmax_phase)
                            
                        elif cmap == SNOM_amplitude and Plot_Definitions.amp_cbar_range is True:
                            if Plot_Definitions.vmin_amp is None: Plot_Definitions.vmin_amp = min_data
                            if Plot_Definitions.vmax_amp is None: Plot_Definitions.vmax_amp = max_data
                            if min_data < Plot_Definitions.vmin_amp: Plot_Definitions.vmin_amp = min_data # update the min and max values in Plot_Definitions if new values are outside of range
                            if max_data > Plot_Definitions.vmax_amp: Plot_Definitions.vmax_amp = max_data
                            vmin = Plot_Definitions.vmin_amp
                            vmax = Plot_Definitions.vmax_amp
                            img = axis.pcolormesh(data, cmap=cmap, vmin=vmin, vmax=vmax)
                        elif cmap == SNOM_height and Plot_Definitions.height_cbar_range is True:
                            if Plot_Definitions.vmin_height is None: Plot_Definitions.vmin_height = min_data # initialize for the first time
                            if Plot_Definitions.vmax_height is None: Plot_Definitions.vmax_height = max_data
                            if min_data < Plot_Definitions.vmin_height: Plot_Definitions.vmin_height = min_data # update the min and max values in Plot_Definitions if new values are outside of range
                            if max_data > Plot_Definitions.vmax_height: Plot_Definitions.vmax_height = max_data
                            vmin = Plot_Definitions.vmin_height
                            vmax = Plot_Definitions.vmax_height
                            img = axis.pcolormesh(data, cmap=cmap, vmin=vmin, vmax=vmax)
                        else:
                            # print('not plotting full range phase')
                            img = axis.pcolormesh(data, cmap=cmap)
                    
                    # legacy method to draw white pixels around masked areas, currently out of service because 
                    # the mask is not stored in the plot variable but for the whole measurement.
                    # repeated calls of the measurement instance would lead to problems
                    '''
                    if (cmap == SNOM_height) and ('_masked' in title) and ('_reduced' not in title):
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
                        try:
                            plt.register_cmap(cmap=map_object)
                        except: pass
                        axis.pcolormesh(white_pixels, cmap='rainbow_alpha')
                    '''
                    
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

        #turn off all unneeded axes
        counter = 0
        for row in range(nrows):
            for col in range(int(np.sqrt(number_of_axis))):
                if  counter >= number_of_subplots > 3 and changed_orientation == False and number_of_subplots != 4: 
                    axis = ax[row, col]
                    axis.axis('off')
                counter += 1

        plt.subplots_adjust(hspace=self.hspace)
        if self.tight_layout is True:
            plt.tight_layout()
        if Plot_Definitions.show_plot is True:
            plt.show()
    
    def Switch_Supplots(self, first_id:int=None, second_id:int=None) -> None:
        """
        This function changes the position of the subplots.
        The first and second id corresponds to the positions of the two subplots which should be switched.
        This function will be repea you are satisfied.

        Args:
            first_id [int]: the first id of the two subplots which should be switched
            second_id [int]: the second id of the two subplots which should be switched
        """
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
        """Add all data contained in dataset as subplots to one figure.
        The data has to be shaped beforehand!
        channels should contain the information which channel is stored at which position in the dataset.
        """
        subplots = []
        for i in range(len(dataset)):
            scalebar = None
            for j in range(len(self.scalebar)):
                if self.channels[i] == self.scalebar[j][0]:
                    scalebar = self.scalebar[j][1]
            subplots.append(self._Add_Subplot(dataset[i], channels[i], scalebar))
        self._Plot_Subplots(subplots)

    def Display_All_Subplots(self) -> None:
        """
        This function displays all the subplots which have been created until this point.
        """
        self._Plot_Subplots(self.all_subplots)

    def Display_Channels(self, channels:list=None) -> None: #, show_plot:bool=True
        """This function displays the channels in memory or the specified ones.
                
        Args:
            channels (list, optional): List of channels to display. If not specified all channels from memory will be plotted. Defaults to None.

        """
        # self.show_plot = show_plot
        if channels == None:
            dataset = self.all_data
            plot_channels_dict = self.channels_label
            # plot_channels = self.channels
        else:
            dataset = []
            plot_channels_dict = []
            # plot_channels = []
            for channel in channels:
                if channel in self.channels:
                    dataset.append(self.all_data[self.channels.index(channel)])
                    plot_channels_dict.append(self.channels_label[self.channels.index(channel)])
                    # plot_channels.append(channel)
                else: 
                    print(f'Channel {channel} is not in memory! Please initiate the channels you want to display first!')
                    print(self.channels)

            # dataset, dict = self._Load_Data(channels)
        self._Display_Dataset(dataset, plot_channels_dict)
        # self._Display_Dataset(dataset, plot_channels)

    def _Gauss_Blurr_Data(self, array, sigma) -> np.array:
        """Applies a gaussian blurr to the specified array, with a specified sigma. The blurred data is returned as a list."""
        return gaussian_filter(array, sigma)

    def Gauss_Filter_Channels(self, channels:list=None, sigma=2):
        """This function will gauss filter the specified channels. If no channels are specified, the ones in memory will be used.

        Args:
            channels (list, optional): List of channels to blurr, if not specified all channels will be blurred. Should not be used for phase. Defaults to None.
            sigma (int, optional): The 'width' of the gauss blurr in pixels, you should scale the data before blurring. Defaults to 2.
        """
        # self._Initialize_Data(channels) # remove initialization and only filter specified channels
        if channels is None:
            channels = self.channels
        self._Write_to_Logfile('gaussian_filter_sigma', sigma)
        
        # start the blurring:
        for channel in channels:
            if channel in self.channels:
                channel_index = self.channels.index(channel)
                # check pixel scaling from channel tag dict for each channel
                pixel_scaling = self.channel_tag_dict[channel_index][Tag_Type.pixel_scaling]
                if pixel_scaling == 1:
                    if Plot_Definitions.show_plot:
                        print(f'The data in channel {channel} is not yet scaled! Do you want to scale the data?')
                        user_input = self._User_Input_Bool()
                        if user_input == True:
                            self.Scale_Channels([channel])
                self.all_data[channel_index] = self._Gauss_Blurr_Data(self.all_data[channel_index], sigma)
                self.channels_label[channel_index] += '_' + self.filter_gauss_indicator
                # self.channels[channel_index] = channel + '_' + self.filter_gauss_indicator
            else: 
                print(f'Channel {channel} is not in memory! Please initiate the channels you want to use first!')

    def Gauss_Filter_Channels_complex(self, channels:list=None, sigma=2) -> None:
        """This function gauss filters the instance channels. For optical channels, amplitude and phase have to be specified!
        Please make shure you scale your data prior to calling this function rather improve the visibility than loosing to much information
                
        Args:
            channels [list]: list of channels to blurr, must contain amplitude and phase of same channels.
            sigma [int]: the sigma used for blurring the data, bigger sigma means bigger blurr radius

        """
        # self._Initialize_Data(channels) # remove initialization and only filter specified channels
        self._Write_to_Logfile('gaussian_filter_complex_sigma', sigma)
        #check if data is scaled, this should be done prior to blurring
        # check pixel scaling from channel tag dict for each channel
            
        if channels is None:
            channels = self.channels
        for channel in channels:
            if channel not in self.channels:
                print(f'Channel {channel} is not in memory! Please initiate the channels you want to use first!')
        channels_to_filter = []
        # if optical channels should be blurred, the according amp and phase data are used to get the complex values and blurr those
        # before backconversion to amp and phase, the realpart could also be returned in future... ToDo
        # print(f'self.channels: {self.channels}')
        # print(f'self.overlayed_amp_channels: {self.overlayed_amp_channels}')
        # print(f'self.overlain_phase_channels: {self.overlain_phase_channels}')
        # print(f'self.corrected_overlain_phase_channels: {self.corrected_overlain_phase_channels}')
        for i in range(len(self.phase_channels)):
            if (self.amp_channels[i] in channels):
                if (self.phase_channels[i] in channels):
                    channels_to_filter.append(self.channels.index(self.amp_channels[i]))
                    channels_to_filter.append(self.channels.index(self.phase_channels[i]))
                elif (self.corrected_phase_channels[i] in channels):
                    channels_to_filter.append(self.channels.index(self.amp_channels[i]))
                    channels_to_filter.append(self.channels.index(self.corrected_phase_channels[i]))
            elif self.overlayed_amp_channels[i] in channels:
                if self.overlain_phase_channels[i] in channels:
                    channels_to_filter.append(self.channels.index(self.overlayed_amp_channels[i]))
                    channels_to_filter.append(self.channels.index(self.overlain_phase_channels[i]))
                elif self.corrected_overlain_phase_channels[i] in channels:
                    channels_to_filter.append(self.channels.index(self.overlayed_amp_channels[i]))
                    channels_to_filter.append(self.channels.index(self.corrected_overlain_phase_channels[i]))
        
        # should not be necessary anymore since backwards channesl are now included in standart channle lists
        # also for backwards direction:
        for i in range(len(self.phase_channels)):
            if (self.backwards_indicator + self.amp_channels[i] in channels):
                if (self.backwards_indicator + self.phase_channels[i] in channels):
                    channels_to_filter.append(self.channels.index(self.backwards_indicator + self.amp_channels[i]))
                    channels_to_filter.append(self.channels.index(self.backwards_indicator + self.phase_channels[i]))
                elif (self.backwards_indicator + self.corrected_phase_channels[i] in channels):
                    channels_to_filter.append(self.channels.index(self.backwards_indicator + self.amp_channels[i]))
                    channels_to_filter.append(self.channels.index(self.backwards_indicator + self.corrected_phase_channels[i]))
        # print(f'channels to filter: {channels_to_filter}')
        # for i in range(len(channels_to_filter)):
        #     print(self.channels[channels_to_filter[i]])
        # if not at least one pair is found:
        if len(channels_to_filter) == 0:
            print('In order to apply the gaussian_filter amplitude and phase of the same channel number must be in the channels list!')
            print('Otherwise only height cannels will be filtered!')
        # add all channels which are not optical to the 'to filter' list
        # print('channels_to_filter:', channels_to_filter)
        for channel in channels:
            # if (channel not in self.amp_channels) and ((self.phase_channels[i] not in self.channels) or (self.corrected_phase_channels[i] not in self.channels)):
            if self.height_indicator in channel:
                channels_to_filter.append(self.channels.index(channel))
                channels_to_filter.append(self.channels.index(channel))# just add twice for now, change later! ToDo
            elif self.real_indicator in channel or self.imag_indicator in channel:
                channels_to_filter.append(self.channels.index(channel))
                channels_to_filter.append(self.channels.index(channel))# just add twice for now, change later! ToDo
            elif self.channels.index(channel) not in channels_to_filter:
                print(f'You wanted to blurr {channel}, but that is not implemented! 1')
        # print('channels_to_filter:', channels_to_filter)
        
        for i in range(int(len(channels_to_filter)/2)):
            # print(f'channel {self.channels[channels_to_filter[2*i]]} is blurred!')
            if (self.channels[channels_to_filter[2*i]] in self.amp_channels) or (self.channels[channels_to_filter[2*i]] in [self.backwards_indicator + element for element in self.amp_channels]) or (self.channels[channels_to_filter[2*i]] in self.overlayed_amp_channels):
                # print(self.channel_tag_dict) # ToDo remove
                # print(self.file_type) # ToDo remove
                # print(self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_scaling] = scaling)
                pixel_scaling_amp = self.channel_tag_dict[channels_to_filter[2*i]][Tag_Type.pixel_scaling]
                pixel_scaling_phase = self.channel_tag_dict[channels_to_filter[2*i+1]][Tag_Type.pixel_scaling]
                if pixel_scaling_amp == 1 and pixel_scaling_phase == 1:
                    if Plot_Definitions.show_plot:# is only false if no plots should be shown or if user inputs are unwanted, eg. for gui
                        print('The data is not yet scaled! Do you want to scale the data?')
                        user_input = self._User_Input_Bool()
                        if user_input == True:
                            self.Scale_Channels([self.channels[channels_to_filter[2*i]], self.channels[channels_to_filter[2*i+1]]])
                amp = self.all_data[channels_to_filter[2*i]]
                phase = self.all_data[channels_to_filter[2*i+1]]
                # compl = np.add(amp*np.cos(phase), 1J*amp*np.sin(phase))
                real = amp*np.cos(phase)
                imag = amp*np.sin(phase)

                # compl_blurred = self._Gauss_Blurr_Data(compl, sigma)
                real_blurred = self._Gauss_Blurr_Data(real, sigma)
                imag_blurred = self._Gauss_Blurr_Data(imag, sigma)
                compl_blurred = np.add(real_blurred, 1J*imag_blurred)
                amp_blurred = np.abs(compl_blurred)
                phase_blurred = self._Get_Compl_Angle(compl_blurred)

                self.all_data[channels_to_filter[2*i]] = amp_blurred
                self.channels_label[channels_to_filter[2*i]] = self.channels_label[channels_to_filter[2*i]] + '_' + self.filter_gauss_indicator
                # self.channels[channels_to_filter[2*i]] += '_' + self.filter_gauss_indicator
                self.all_data[channels_to_filter[2*i+1]] = phase_blurred
                self.channels_label[channels_to_filter[2*i+1]] = self.channels_label[channels_to_filter[2*i+1]] + '_' + self.filter_gauss_indicator
                # self.channels[channels_to_filter[2*i+1]] += '_' + self.filter_gauss_indicator

            elif self.height_indicator in self.channels[channels_to_filter[2*i]]:
                pixel_scaling = self.channel_tag_dict[channels_to_filter[2*i]][Tag_Type.pixel_scaling]
                if pixel_scaling == 1:
                    if Plot_Definitions.show_plot:
                        print('The data is not yet scaled! Do you want to scale the data?')
                        user_input = self._User_Input_Bool()
                        if user_input == True:
                            self.Scale_Channels([self.channels[channels_to_filter[2*i]]])
                height = self.all_data[channels_to_filter[2*i]]
                height_blurred = self._Gauss_Blurr_Data(height, sigma)
                self.all_data[channels_to_filter[2*i]] = height_blurred
                self.channels_label[channels_to_filter[2*i]] = self.channels_label[channels_to_filter[2*i]] + '_' + self.filter_gauss_indicator
                # self.channels[channels_to_filter[2*i]] += '_' + self.filter_gauss_indicator
            elif self.real_indicator in self.channels[channels_to_filter[2*i]] or self.imag_indicator in self.channels[channels_to_filter[2*i]]:
                pixel_scaling = self.channel_tag_dict[channels_to_filter[2*i]][Tag_Type.pixel_scaling]
                if pixel_scaling == 1:
                    if Plot_Definitions.show_plot:
                        print('The data is not yet scaled! Do you want to scale the data?')
                        user_input = self._User_Input_Bool()
                        if user_input == True:
                            self.Scale_Channels([self.channels[channels_to_filter[2*i]]])
                data = self.all_data[channels_to_filter[2*i]]
                data_blurred = self._Gauss_Blurr_Data(data, sigma)
                self.all_data[channels_to_filter[2*i]] = data_blurred
                self.channels_label[channels_to_filter[2*i]] = self.channels_label[channels_to_filter[2*i]] + '_' + self.filter_gauss_indicator
                # self.channels[channels_to_filter[2*i]] += '_' + self.filter_gauss_indicator
            
            else:
                print(f'You wanted to blurr {self.channels[channels_to_filter[2*i]]}, but that is not implemented! 2')
                
    def _Get_Compl_Angle(self, compl_number_array) -> np.array:
        """This function returns the angles of a clomplex number array."""
        YRes = len(compl_number_array)
        XRes = len(compl_number_array[0])
        realpart = compl_number_array.real
        imagpart = compl_number_array.imag
        r = np.sqrt(pow(imagpart, 2) + pow(realpart, 2))
        phase = np.arctan2(r*imagpart, r*realpart) #returns values between -pi and pi, add pi for the negative values
        for i in range(YRes):
            for j in range(XRes):
                if phase[i][j] < 0:
                    phase[i][j]+=2*np.pi
        return phase

    def _Fourier_Filter_Array(self, complex_array) -> np.array:
        '''
        Takes a complex array and returns the fourier transformed complex array
        '''
        FS_compl = np.fft.fftn(complex_array)
        return FS_compl
    
    def Fourier_Filter_Channels(self, channels:list=None) -> None:
        """This function applies the Fourier filter to all data in memory or specified channels
                
        Args:
            channels [list]: list of channels, will override the already existing channels
        """
        self._Initialize_Data(channels)
        self._Write_to_Logfile('fourier_filter', True)
        channels_to_filter = []
        for i in range(len(self.amp_channels)):
            if (self.amp_channels[i] in self.channels) and (self.phase_channels[i] in self.channels):
                channels_to_filter.append(self.channels.index(self.amp_channels[i]))
                channels_to_filter.append(self.channels.index(self.phase_channels[i]))
            else:
                print('In order to apply the fourier_filter amplitude and phase of the same channel number must be in the channels list!')
        for i in range(int(len(channels_to_filter)/2)):
            amp = self.all_data[channels_to_filter[i]]
            phase = self.all_data[channels_to_filter[i+1]]
            compl = np.add(amp*np.cos(phase), 1J*amp*np.sin(phase))
            FS_compl = self._Fourier_Filter_Array(compl)
            FS_compl_abs = np.absolute(FS_compl)
            FS_compl_angle = self._Get_Compl_Angle(FS_compl)
            self.all_data[channels_to_filter[i]] = np.log(np.abs(np.fft.fftshift(FS_compl_abs))**2)
            self.channels_label[channels_to_filter[i]] = self.channels_label[channels_to_filter[i]] + '_fft'
            self.all_data[channels_to_filter[i+1]] = FS_compl_angle
            self.channels_label[channels_to_filter[i+1]] = self.channels_label[channels_to_filter[i+1]] + '_fft'

    def _Create_Header(self, channel, data=None, filetype='gsf'):
        # data = self.all_data[self.channels.index(channel)]
        # load data instead, because sometimes the channel is not in memory
        if data is None:
            # channel is not in memory, so the standard values will be used
            data = self._Load_Data([channel])[0][0]
            # XReal, YReal = self.measurement_tag_dict[Tag_Type.scan_area]# change to self.channel_dat_dict?
            XReal, YReal = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.scan_area]# change to self.channel_dat_dict?
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
        header += f'XRes={int(XRes)}\nYRes={int(YRes)}\n'
        header += f'XReal={round(XReal,9)}\nYReal={round(YReal,9)}\n'
        header += f'XYUnits=m\n'
        header += f'XOffset={round(XOffset*pow(10, -6),9)}\nYOffset={round(YOffset*pow(10, -6),9)}\n'
        header += f'Rotation={round(rotation)}\n'
        if self.height_indicator in channel:
            header += 'ZUnits=m\n'
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
        """This function is ment to save all specified channels to external .gsf files.
        
        Args:
            channels [list]:    list of the channels to be saved, if not specified, all channels in memory are saved
                                Careful! The data will be saved as it is right now, so with all the manipulations.
                                Therefor the data will have an '_manipulated' appendix in the filename.
        """
        self._Write_to_Logfile('save_to_gsf_appendix', appendix)
        if channels == None:
            channels = self.channels
        for channel in channels:
            # filepath = self.directory_name + '/' + self.filename + ' ' + channel + appendix + '.gsf'
            filepath = self.directory_name / Path(self.filename.name + f' {channel}{appendix}.gsf')
            data = self.all_data[self.channels.index(channel)]
            XRes = len(data[0])
            YRes  = len(data)
            header, NUL = self._Create_Header(channel, data)
            file = open(filepath, 'bw')
            file.write(header.encode('utf-8'))
            file.write(NUL) # the NUL marks the end of the header and konsists of 0 characters in the first dataline
            if self.height_indicator in channel:
                for y in range(YRes):
                    for x in range(XRes):
                        file.write(pack('f', round(data[y][x],5)*pow(10,-9)))
            else:
                for y in range(YRes):
                    for x in range(XRes):
                        file.write(pack('f', round(data[y][x], 5)))
            file.close()
            print(f'successfully saved channel {channel} to .gsf')
            print(filepath)

    def Save_to_txt(self, channels:list=None, appendix:str='_manipulated'):
        """This function is ment to save all specified channels to external .txt files.
        
        Args:
            channels [list]:    list of the channels to be saved, if not specified, all channels in memory are saved
                                Careful! The data will be saved as it is right now, so with all the manipulations.
                                Therefor the data will have an '_manipulated' appendix in the filename.
        """
        self._Write_to_Logfile('save_to_txt_appendix', appendix)
        if channels == None:
            channels = self.channels
        for channel in channels:
            # filepath = self.directory_name + '/' + self.filename + ' ' + channel + appendix + '.txt'
            filepath = self.directory_name / Path(self.filename.name + f' {channel}{appendix}.txt')
            data = self.all_data[self.channels.index(channel)]
            XRes = len(data[0])
            YRes  = len(data)
            header, NUL = self._Create_Header(channel, data, 'txt')
            file = open(filepath, 'w')
            file.write(header)
            # file.write(NUL) # the NUL marks the end of the header and konsists of 0 characters in the first dataline
            for y in range(YRes):
                for x in range(XRes):
                    file.write(f'{round(data[y][x], 5)} ')
            file.close()
            print(f'successfully saved channel {channel} to .txt')
            print(filepath)
    
    def _Create_Synccorr_Preview(self, channel, wavelength, nouserinput=False) -> None:
        """
        This function is part of the Synccorrection and creates a preview of the corrected data.
        channel specifies which channel should be used for the preview.
        Wavelength must be given in µm.
        Scanangle is the rotation angle of the scan in radians.
        """
        scanangle = self.measurement_tag_dict[Tag_Type.rotation]*np.pi/180
        phasedir_positive = 1
        phasedir_negative = -1
        # phase_data = self._Load_Data([channel])[0][0]
        phase_data = self.all_data[self.channels.index(channel)]
        YRes = len(phase_data)
        XRes = len(phase_data[0])
        phase_positive = np.zeros((YRes, XRes))
        phase_negative = np.zeros((YRes, XRes))
        phase_no_correction = np.zeros((YRes, XRes))
        for y in range(0,YRes):
            for x in range(0,XRes):
                xreal=x*self.XReal/XRes
                yreal=y*self.YReal/YRes
                #phase accumulated by movement of parabolic mirror only depends on 'x' direction
                phase_no_correction[y][x] = phase_data[y][x]# + np.pi
                phase_positive[y][x] = np.mod(phase_data[y][x] - phasedir_positive*(np.cos(-scanangle)*xreal + np.sin(-scanangle)*yreal)/wavelength*2*np.pi, 2*np.pi)
                # phase_positive[y][x] = np.mod(phase_data[y][x] + np.pi - phasedir_positive*(np.cos(-scanangle)*xreal + np.sin(-scanangle)*yreal)/wavelength*2*np.pi, 2*np.pi)
                phase_negative[y][x] = np.mod(phase_data[y][x] - phasedir_negative*(np.cos(-scanangle)*xreal + np.sin(-scanangle)*yreal)/wavelength*2*np.pi, 2*np.pi)
                # phase_negative[y][x] = np.mod(phase_data[y][x] + np.pi - phasedir_negative*(np.cos(-scanangle)*xreal + np.sin(-scanangle)*yreal)/wavelength*2*np.pi, 2*np.pi)
        #create plots of the uncorrected and corrected images
        subplots = []
        subplots.append(self._Add_Subplot(phase_no_correction, channel))
        subplots.append(self._Add_Subplot(phase_positive, channel + '_positive'))
        subplots.append(self._Add_Subplot(phase_negative, channel + '_negative'))
        self._Plot_Subplots(subplots)
        # remove the preview subplots from the subplot memory after plotting
        self.Remove_Last_Subplots(3)
        #ask the user to chose a correction direction
        if nouserinput is False:
            phasedir = self._Gen_From_Input_Phasedir()
            return phasedir
        #start the correction
        # self.Synccorrection(wavelength, phasedir)

    def Synccorrection(self, wavelength:float, phasedir:int=None) -> None:
        """This function corrects all the phase channels for the linear phase gradient
        which stems from the synchronized measurement mode.
        The wavelength must be given in µm.
        The phasedir is either 1 or -1. If you are unshure about the direction just leave the parameter out.
        You will be shown a preview for both directions and then you must choose the correct one.
                
        Args:
            wavelenght (float): please enter the wavelength in µm.
            phasedir (int): the phase direction, leave out if not known and you will be prompted with a preview and can select the appropriate direction.

        """
        if self.autoscale == True:
            print('careful! The synccorretion does not work when autoscale is enabled.')
            exit()
        all_channels = self.phase_channels + self.amp_channels
        # print(all_channels)
        self._Initialize_Data(all_channels)
        # print(self.channels)
        scanangle = self.measurement_tag_dict[Tag_Type.rotation]*np.pi/180
        if phasedir == None:
            phasedir = self._Create_Synccorr_Preview(self.preview_phasechannel, wavelength)
        self._Write_to_Logfile('synccorrection_wavelength', wavelength)
        self._Write_to_Logfile('synccorrection_phasedir', phasedir)
        # all_channels = self.all_channels
        # dataset, dict = self._Load_Data(all_channels)
        header, NUL = self._Create_Header(self.preview_phasechannel) # channel for header just important to distinguish z axis unit either m or nothing
        for channel in self.phase_channels:
            i = self.phase_channels.index(channel)
            # phasef = open(self.directory_name + '/' + self.filename + ' ' + channel + '_corrected.gsf', 'bw')
            phasef = open(self.directory_name / Path(self.filename.name + f' {channel}_corrected.gsf'), 'bw')
            # realf = open(self.directory_name + '/' + self.filename + ' ' + self.real_channels[i] + '_corrected.gsf', 'bw')
            realf = open(self.directory_name / Path(self.filename.name + f' {self.real_channels[i]}_corrected.gsf'), 'bw')
            phasef.write(header.encode('utf-8'))
            realf.write(header.encode('utf-8'))
            phasef.write(NUL) # add NUL terminator
            realf.write(NUL)
            # XRes = len(dataset[0][0])
            # XRes = len(dataset[0])
            for y in range(0,self.YRes):
                for x in range(0,self.XRes):
                    #convert pixel number to realspace coordinates in µm
                    xreal=x*self.XReal/self.XRes
                    yreal=y*self.YReal/self.YRes
                    #open the phase, add pi to change the range from 0 to 2 pi and then substract the linear phase gradient, which depends on the scanangle!
                    # amppixval = dataset[2*i][y][x]
                    amppixval = self.all_data[self.channels.index(self.amp_channels[i])][y][x]
                    phasepixval = self.all_data[self.channels.index(self.phase_channels[i])][y][x]
                    phasepixval_corr = np.mod(phasepixval + np.pi - phasedir*(np.cos(-scanangle)*xreal + np.sin(-scanangle)*yreal)/wavelength*2*np.pi, 2*np.pi)
                    realpixval = amppixval*np.cos(phasepixval_corr)
                    phasef.write(pack("f",phasepixval_corr))
                    realf.write(pack("f",realpixval))
            phasef.close()
            realf.close()

    def _Gen_From_Input_Phasedir(self) -> int:
        """
        This function asks the user to input a phase direction, input must be either n or p, for negative or positive respectively.
        """
        phasedir = input('Did you prefer the negative or positive phase direction? Please enter either \'n\' or \'p\'\n')
        if phasedir == 'n':
            return -1
        elif phasedir == 'p':
            return 1
        else:
            print('Wrong letter! Please try again.')
            self._Gen_From_Input_Phasedir()
    
    def _Get_Channel_Scaling(self, channel_id) -> int :
        """This function checks if an instance channel is scaled and returns the scaling factor."""
        channel_yres = len(self.all_data[channel_id])
        # channel_xres = len(self.all_data[channel_id][0])
        return int(channel_yres/self.YRes)

    def _Create_Height_Mask_Preview(self, mask_array) -> None:
        """This function creates a preview of the height masking.
        The preview is based on all channels in the instance"""
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
        """This function asks the user to input yes or no and returns a boolean value."""
        user_input = input('Please type y for yes or n for no. \nInput: ')
        if user_input == 'y':
            user_bool = True
        elif user_input == 'n':
            user_bool = False
        return user_bool

    def _User_Input(self, message:str):
        """This function confronts the user with the specified message and returns the user input

        Args:
            message (str): the message to display
        """
        return input(message)

    def _Create_Mask_Array(self, height_data, threshold) -> np.array:
        """This function takes the height data and a threshold value to create a mask array containing 0 and 1 values.
        """
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
        """This function returns the height threshold value dependent on the user input"""
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
        """
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
        """
        if export == True:
            channels = self.all_channels
        self._Initialize_Data(channels)
        if (mask_channel == None) or (mask_channel not in self.channels):
            if self.height_channel in self.channels:
                height_data = self.all_data[self.channels.index(self.height_channel)]
                if 'leveled' not in self.channels_label[self.channels.index(self.height_channel)]:
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
        self.mask_array = mask_array # todo, mask array must be saved as part of the image, otherwise multiple measurement creations will use the same mask
        if export == True:
            # open files for the masked data:
            for channel in channels:
                header, NUL= self._Create_Header(channel)
                data, trash = self._Load_Data([channel])
                # datafile = open("".join([self.directory_name,"/",self.filename," ",channel,"_masked.gsf"]),"bw")
                datafile = open(self.directory_name / Path(self.filename.name + f' {channel}_masked.gsf'),"bw")
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
                if self.height_channel not in self.channels_label[i]:
                    self.all_data[i] = np.multiply(dataset[i], mask_array)
                self.channels_label[i] = self.channels_label[i] + '_masked'

    def _Check_Pixel_Position(self, xres, yres, x, y) -> bool:
        """This function checks if the pixel position is within the bounds"""
        if x < 0 or x > xres:
            return False
        elif y < 0 or y > yres:
            return False
        else: return True

    def _Get_Mean_Value(self, data, x_coord, y_coord, zone) -> float:
        """This function returns the mean value of the pixel and its nearest neighbors.
        The zone specifies the number of neighbors. 1 means the pixel and the 8 nearest pixels.
        2 means zone 1 plus the next 16, so a total of 25 with the pixel in the middle. 
        """
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
        ax.pcolormesh(height_data, cmap=SNOM_height)
        klicker = clicker(ax, ["event"], markers=["x"])
        ax.legend()
        ax.axis('scaled')
        plt.title('3 Point leveling: please click on three points\nto specify the underground plane.')
        if Plot_Definitions.show_plot:
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
    
    def _level_height_data(self, height_data, klick_coordinates, zone):
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

    def _get_klicker_coordinates(data, cmap):
        fig, ax = plt.subplots()
        ax.pcolormesh(data, cmap=cmap)
        klicker = clicker(ax, ["event"], markers=["x"])
        ax.legend()
        ax.axis('scaled')
        plt.title('3 Point leveling: please click on three points\nto specify the underground plane.')
        # if Plot_Definitions.show_plot:
        plt.show()
        klicker_coords = klicker.get_positions()['event'] #klicker returns a dictionary for the events
        klick_coordinates = [[round(element[0]), round(element[1])] for element in klicker_coords]
        return klick_coordinates

    def _Height_Levelling_3Point_forGui(self, height_data, zone=1) -> np.array:
        klick_coordinates = self._get_klicker_coordinates(height_data, SNOM_height)
        if len(klick_coordinates) != 3:
            print('You need to specify 3 point coordinates! Data was not leveled!')
            return height_data
        #     klick_coordinates = get_coordinates()
        #     user_input = self._User_Input_Bool()
        #     if user_input == True:
        #         self._Height_Levelling_3Point(zone)
        #     else:
        #         exit()
        # for the 3 point coordinates the height data is calculated over a small area around the clicked pixels to reduce deviations due to noise
        self._Write_to_Logfile('height_leveling_coordinates', klick_coordinates)
        return self._level_height_data(klick_coordinates, zone)


    def _Level_Phase_Slope(self, data, slope) -> np.array:
        """This function substracts a linear phase gradient in y direction from the specified phase data.
        """
        yres = len(data)
        xres = len(data[0])
        for y in range(yres):
            for x in range(xres):
                data[y][x] -= y*slope
        return self._Shift_Phase_Data(data, 0)

    def Correct_Phase_Drift(self, channels:list=None, export:bool=False, phase_slope=None, zone:int=1) -> None:
        """This function asks the user to click on two points which should have the same phase value.
        Only the slow drift in y-direction will be compensated. Could in future be extended to include a percentual drift compensation along the x-direction.
        But should usually not be necessary.
                
        Args:
            channels [list]: list of channels, will override the already existing channels
            export [bool]: do you want to aply the correction to all phase channels and export them?
            phase_slope [float]: if you already now the phase slope you can enter it, otherwise leave it out
                                and it will prompt you with a preview to select two points to calculate the slope from
            zone [int]: defines the area which is used to calculate the mean around the click position in the preview,
                        0 means only the click position, 1 means the nearest 9 ...
        """
        # zone = int(zone*self.scaling_factor/4) #automatically enlargen the zone if the data has been scaled by more than a factor of 4
        self._Initialize_Data(channels)
        phase_data = None
        if self.preview_phasechannel in self.channels:
            phase_data = np.copy(self.all_data[self.channels.index(self.preview_phasechannel)])
            phase_channel = self.preview_phasechannel
        else:
            phase_data = self._Load_Data([self.preview_phasechannel])[0][0]
            phase_channel = self.preview_phasechannel
        # for i in range(len(self.channels)):
            # if '3P' in self.channels[i]:
            #     phase_data = np.copy(self.all_data[i])
            #     phase_channel = self.channels[i]
            # elif ('2P' in self.channels[i]) and ('3P' not in self.channels[i]):
            #     phase_data = np.copy(self.all_data[i])
            #     phase_channel = self.channels[i]
            # elif ('4P' in self.channels[i]) and ('3P' not in self.channels[i])  and ('2P' not in self.channels[i]):
            #     phase_data = np.copy(self.all_data[i])
            #     phase_channel = self.channels[i]
        if export == True:
            # ToDo
            # do something with the phase slope...
            print('You want to export a phase slope correction, but nothing happens!')
            pass
        else:
            if phase_slope != None:
                #level all phase channels in memory...
                self._Write_to_Logfile('phase_driftcomp_slope', phase_slope)
                for i in range(len(self.channels)):
                    if 'P' in self.channels[i]:
                        self.all_data[i] = self._Level_Phase_Slope(self.all_data[i], phase_slope)
                        self.channels_label[i] += '_driftcomp'
            else:
                fig, ax = plt.subplots()
                img = ax.pcolormesh(phase_data, cmap=SNOM_phase)
                klicker = clicker(ax, ["event"], markers=["x"])
                ax.invert_yaxis()
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = plt.colorbar(img, cax=cax)
                cbar.ax.get_yaxis().labelpad = 15
                cbar.ax.set_ylabel('phase', rotation=270)
                ax.legend()
                ax.axis('scaled')
                plt.title('Phase leveling: please click on two points\nto specify the phase drift.')
                plt.show()
                klicker_coords = klicker.get_positions()['event'] #klicker returns a dictionary for the events
                klick_coordinates = [[round(element[0]), round(element[1])] for element in klicker_coords]
                if len(klick_coordinates) != 2:
                    print('You must specify two points which should have the same phase, along the y-direction')
                    print('Do you want to try again?')
                    user_input = self._User_Input_Bool()
                    if user_input == True:
                        self.Correct_Phase_Drift(channels, export, None)
                    else: 
                        exit()
                mean_values = [self._Get_Mean_Value(phase_data, klick_coordinates[i][0], klick_coordinates[i][1], zone) for i in range(len(klick_coordinates))]
                #order points from top to bottom
                if klick_coordinates[0][1] > klick_coordinates[1][1]:
                    second_corrd = klick_coordinates[0]
                    second_mean = mean_values[0]
                    klick_coordinates[0] = klick_coordinates[1]
                    klick_coordinates[1] = second_corrd
                    mean_values[0] = mean_values[1]
                    mean_values[1] = second_mean
                phase_slope = (mean_values[1] - mean_values[0])/(klick_coordinates[1][1] - klick_coordinates[0][1])
                leveled_phase_data = self._Level_Phase_Slope(phase_data, phase_slope)
                fig, ax = plt.subplots()
                ax.pcolormesh(leveled_phase_data, cmap=SNOM_phase)
                ax.invert_yaxis()
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = plt.colorbar(img, cax=cax)
                cbar.ax.get_yaxis().labelpad = 15
                cbar.ax.set_ylabel('phase', rotation=270)
                ax.legend()
                ax.axis('scaled')
                plt.title('Leveled Pase: ' + phase_channel)
                plt.show()
                print('Are you satisfied with the phase leveling?')
                user_input = self._User_Input_Bool()
                if user_input == True:
                    #use the phase slope to level all phase channels in memory
                    self.Correct_Phase_Drift(None, False, phase_slope)
                else:
                    print('Do you want to repeat the leveling?')
                    user_input = self._User_Input_Bool()
                    if user_input == True:
                        #start the leveling process again
                        self.Correct_Phase_Drift()
                    else:
                        exit()


    def Level_Height_Channels(self, channels:list=None) -> None:
        """This function levels all height channels which are either user specified or in the instance memory.
        The leveling will prompt the user with a preview to select 3 points for getting the coordinates of the leveling plane.
        
        Args:
            channels (list, optional): List of channels to level. If not specified all channels in memory will be used. Defaults to None.
        """
        # self._Initialize_Data(channels)
        if channels is None:
            channels = self.channels
        for channel in channels:
            if self.height_indicator in channel:
                self.all_data[self.channels.index(channel)] = self._Height_Levelling_3Point(self.all_data[self.channels.index(channel)])
                self.channels_label[self.channels.index(channel)] += '_leveled' 

    def Level_Height_Channels_forGui(self, channels:list=None):# todo not used?
        """This function levels all height channels which are either user specified or in the instance memory.
        The leveling will prompt the user with a preview to select 3 points for getting the coordinates of the leveling plane.
        This function is specifically for use with GUI.
        
        Args:
            channels (list, optional): List of channels to level. If not specified all channels in memory will be used. Defaults to None.
        """
        # self._Initialize_Data(channels)
        if channels is None:
            channels = self.channels
        for channel in channels:
            if self.height_indicator in channel:
                self.all_data[self.channels.index(channel)] = self._Height_Levelling_3Point_forGui(self.all_data[self.channels.index(channel)])
                self.channels_label[self.channels.index(channel)] += '_leveled' 

    def _Shift_Phase_Data(self, data, shift) -> np.array:
        """This function adds a phaseshift to the specified phase data. The phase data is automatically kept in the 0 to 2 pi range.
        Could in future be extended to show a live view of the phase data while it can be modified by a slider...
        e.g. by shifting the colorscale in the preview rather than the actual data..."""
        yres = len(data)
        xres = len(data[0])
        for y in range(yres):
            for x in range(xres):
                data[y][x] = (data[y][x] + shift) % (2*np.pi)
        return data

    def Shift_Phase(self, shift:float=None, channels:list=None) -> None:
        """This function will prompt the user with a preview of the first phase channel in memory.
        Under the preview is a slider, by changing the slider value the phase preview will shift accordingly.
        If you are satisfied with the shift, hit the 'accept' button. The preview will close and the shift will
        be applied to all phase channels in memory.

        Args:
            shift (float, optional): If you know the shift value already, you can enter values between 0 and 2*Pi
            channels (list, optional): List of channels to apply the shift to, only phase channels will be shifted though.
                If not specified all channels in memory will be used. Defaults to None.
        """
        if channels is None:
            channels = self.channels
        # self._Initialize_Data(channels)
        if shift == None:
            shift_known = False
        else:
            shift_known = True
        if shift_known == False:
            if self.preview_phasechannel in channels:
                    phase_data = np.copy(self.all_data[self.channels.index(self.preview_phasechannel)])
            else:
                # check if corrected phase channel is present
                # just take the first phase channel in memory
                for channel in channels:
                    if self.phase_indicator in channel:
                        phase_data = np.copy(self.all_data[self.channels.index(channel)])
                        # print(len(phase_data))
                        # print(len(phase_data[0]))
                        break
            shift = Get_Phase_Offset(phase_data)
            print('The phase shift you chose is:', shift)
            shift_known = True

        # export shift value to logfile
        self._Write_to_Logfile('phase_shift', shift)
        # shift all phase channels in memory
        # could also be implemented to shift each channel individually...
        for channel in channels:
            if self.phase_indicator in channel:
                self.all_data[self.channels.index(channel)] = self._Shift_Phase_Data(self.all_data[self.channels.index(channel)], shift)
                
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
        # print("fit succsessful")
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
        """This function corrects the drift of the piezo motor. As of now it needs to be fitted to a region of the sample which is assumed to be straight.
        In the future this could be implemented with a general map containing the distortion created by the piezo motor, if it turns out to be constant...
        Anyways, you will be prompted with a preview of the height data, please select an area of the scan with only one 'straight' waveguide. 
        The bounds for the fitting routine are based on the lower and upper limit of this selection.

        Careful! Will not yet affect the scan size, so the pixelsize will be altered... ToDo
        
        Args:  
            channels [list]: list of channels, will override the already existing channels
        
        """
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
        cmap = SNOM_height
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
            self.channels_label[i] += '_shifted'

    def Cut_Channels(self, channels:list=None, preview_channel:str=None, autocut:bool=False, coords:list=None, reset_mask:bool=False) -> None:
        """This function cuts the specified channels to the specified region. If no coordinates are specified you will be prompted with a window to select an area.
        If you created a mask previously for this instance the old mask will be reused! Otherwise you should manually change the reset_mask parameter to True.

        Args:
            channels (list, optional): List of channels you want to cut. If not specified all channels in memory will be cut. Defaults to None.
            preview_channel (str, optional): The channel to display for the area selection. If not specified the height channel will be used if it is in memory,
                otherwise the first of the specified channels will be used. Defaults to None
            autocut (bool, optional): If set to 'True' the program will automatically try to remove zero lines and columns, which can result from masking.
            coords (list, optional): If you already now the coordinates ([[x1,y1], [x2,y2], [x3,y3], [x4,y4]]) to which you want to cut your data. Defaults to None.
            reset_mask (bool, optional): If you dont want to reuse an old mask set to True. Defaults to False.
        """
        # self._Initialize_Data(channels)
        if channels is None:
            channels = self.channels # if nothing is specified, the cut will be applied to all channels in memory!
        # check if height channel in channels and apply mask to it, until now it has not been masked in order to show the mask in the image
        if preview_channel is None:
            if (self.height_channel in channels):
                preview_channel = self.height_channel
            else:
                preview_channel = channels[0]

        # apply the already existing mask if possible.  
        if reset_mask == False:  
            if (len(self.mask_array) > 0):
                for channel in channels:
                    index = self.channels.index(channel)
                    self.all_data[index] = np.multiply(self.all_data[index], self.mask_array)
                    # self.channels[index] += '_reduced'
            else:
                print('There does not seem to be an old mask...')
        else:
            if autocut == True:
                self._Auto_Cut_Channels(channels)
                self._Write_to_Logfile('auto_cut', True)
            else:
                # if self.height_channel in self.channels:
                #     data = self.all_data[self.channels.index(self.height_channel)]
                #     channel = self.height_channel
                # else:
                #     data = self.all_data[0]
                #     channel = self.channels[0]
                data = self.all_data[self.channels.index(preview_channel)]
                # get the coordinates of the selection rectangle
                if coords is None:
                    coords = Select_Rectangle(data, preview_channel)
                self._Write_to_Logfile('cut_coords', coords)
                # use the selection to create a mask and multiply to all channels, then apply auto_cut function
                yres = len(data)
                xres = len(data[0])
                self.mask_array = np.zeros((yres, xres))
                for y in range(yres):
                    if y in range(coords[0][1], coords[1][1]):
                        for x in range(xres):
                            if x in range(coords[0][0], coords[1][0]):
                                self.mask_array[y][x] = 1
                for channel in channels:
                    index = self.channels.index(channel)
                    # set all values outside of the mask to zero and then cut all zero away from the outside with _Auto_Cut_Channels(channels)
                    self.all_data[index] = np.multiply(self.all_data[index], self.mask_array)
                    # self.channels[index] += '_reduced'
                self._Auto_Cut_Channels(channels)
        
    def _Auto_Cut_Channels(self, channels:list=None) -> None:
        """This function automatically cuts away all rows and lines which are only filled with zeros.
        This function applies to all channels in memory.
        """
        if channels is None:
            channels = self.channels
        
        # get the new size of the reduced channels
        reduced_data = self._Auto_Cut_Data(self.all_data[0])
        yres = len(reduced_data)
        xres = len(reduced_data[0])
        # copy old data to local variable
        all_data = self.all_data
        # reinitialize self.all_data, all channels must have the same size
        # self.all_data = np.zeros((len(all_data), yres, xres))
        # self.all_data = []
        # for i in range(len(self.channels)):
        #     reduced_data = self._Auto_Cut_Data(all_data[i])
        #     self.all_data.append(reduced_data)
        #     self.channels_label[i] += '_reduced'
        for channel in channels:
            index = self.channels.index(channel)
            # get the old size of the data
            xres, yres = self.channel_tag_dict[index][Tag_Type.pixel_area]
            xreal, yreal = self.channel_tag_dict[index][Tag_Type.scan_area]
            self.all_data[index] = self._Auto_Cut_Data(self.all_data[index])
            xres_new = len(self.all_data[index][0])
            yres_new = len(self.all_data[index])
            xreal_new = xreal*xres_new/xres
            yreal_new = yreal*yres_new/yres
            # save new resolution and scan area in channel tag dict:
            self.channel_tag_dict[index][Tag_Type.pixel_area] = [xres_new, yres_new]
            self.channel_tag_dict[index][Tag_Type.scan_area] = [xreal_new, yreal_new]
            # add new appendix to channel
            self.channels_label[index] += '_reduced'
        self._Write_to_Logfile('cut', 'autocut')

    def _Auto_Cut_Data(self, data) -> np.array:
        """This function cuts the data and removes zero values from the outside."""
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
        """Adds a scalebar to all specified channels.
        Args:
            channels (list): List of channels the scalebar should be added to.
            various definitions for the scalebar, please look up 'matplotlib_scalebar.scalebar' for more information
        """
        
        # scalebar = ScaleBar(dx, units, dimension, label, length_fraction, height_fraction, width_fraction,
            # location, loc, pad, border_pad, sep, frameon, color, box_color, box_alpha, scale_loc,
            # label_loc, font_properties, label_formatter, scale_formatter, fixed_value, fixed_units, animated, rotation)
        
        
        count = 0
        for channel in self.channels:
            XRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area][0]
            XReal = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.scan_area][0]
            pixel_scaling = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_scaling]
            # dx = XReal/(XRes*pixel_scaling)
            dx = XReal/(XRes)
            scalebar_var = [dx, units, dimension, label, length_fraction, height_fraction, width_fraction,
                            location, loc, pad, border_pad, sep, frameon, color, box_color, box_alpha, scale_loc,
                            label_loc, font_properties, label_formatter, scale_formatter, fixed_value, fixed_units, animated, rotation]
            if (channel in channels) or (len(channels)==0):
                self.scalebar.append([channel, scalebar_var])                
            else:
                self.scalebar.append([channel, None])                
            count += 1

    def Rotate_90_deg(self, orientation:str = 'right'):
        """This function will rotate all data in memory by 90 degrees.

        Args:
            orientation (str, optional): rotate clockwise ('right') or counter clockwise ('left'). Defaults to 'right'.
        """
        # self._Write_to_Logfile('rotate_90_deg', orientation)
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

    def _Get_Positions_from_Plot(self, channel, data, coordinates:list=None, orientation=None) -> list:
        if self.phase_indicator in channel:
            cmap = SNOM_phase
        elif self.amp_indicator in channel:
            cmap = SNOM_amplitude
        elif self.height_indicator in channel:
            cmap = SNOM_height

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
        if coordinates != None and orientation != None:
            self._Plot_Profile_Lines(data, ax, coordinates, orientation)
        plt.title('Please select one or more points to continue.')
        plt.tight_layout()
        plt.show()
        klicker_coords = klicker.get_positions()['event'] #klicker returns a dictionary for the events
        klick_coordinates = [[round(element[0]), round(element[1])] for element in klicker_coords]
        return klick_coordinates

    def _Get_Profile(self, data, coordinates:list, orientation:Definitions, width:int) -> list:
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

    def Select_Profile(self, profile_channel:str, preview_channel:str=None, orientation:Definitions=Definitions.vertical, width:int=10, phase_orientation:int=1, coordinates:list=None):
        """This function lets the user select a profile with given width in pixels and displays the data.

        Args:
            profile_channel (str): channel to use for profile data extraction
            preview_channel (str, optional): channel to preview the profile positions. If not specified the height channel will be used for that. Defaults to None.
            orientation (Definitions, optional): profiles can be horizontal or vertical. Defaults to Definitions.vertical.
            width (int, optional): width of the profile in pixels, will calculate the mean. Defaults to 10.
            phase_orientation (int, optional): only relevant for phase profiles. Necessary for the flattening to work properly. Defaults to 1.
            coordinates (list, optional): if you already now the position of your profile you can also specify the coordinates and skip the selection. Defaults to None.
        """
        if preview_channel is None:
            preview_channel = self.height_channel
        if coordinates == None:
            previewdata = self.all_data[self.channels.index(preview_channel)]
            coordinates = self._Get_Positions_from_Plot(preview_channel, previewdata)
            # print('The coordinates you selected are:', coordinates)

        profiledata = self.all_data[self.channels.index(profile_channel)]

        cmap = SNOM_phase
        fig, ax = plt.subplots()
        img = ax.pcolormesh(profiledata, cmap=cmap)
        ax.invert_yaxis()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('phase', rotation=270)
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

        profiles = self._Get_Profile(profiledata, coordinates, orientation, width)
        for profile in profiles:
            xvalues = np.linspace(0, 10, len(profile))
            plt.plot(xvalues, profile, 'x')
        plt.title('Phase profiles')
        plt.tight_layout()
        plt.show()

        flattened_profiles = [phase_analysis.Flatten_Phase_Profile(profile, phase_orientation) for profile in profiles]
        for profile in flattened_profiles:
            xvalues = np.linspace(0, 10, len(profile))
            plt.plot(xvalues, profile)
        plt.title('Flattened phase profiles')
        plt.tight_layout()
        plt.show()

        difference_profile = phase_analysis.Get_Profile_Difference(profiles[0], profiles[1])
        # difference_profile = Get_Profile_Difference(flattened_profiles[0], flattened_profiles[1])
        xres, yres = self.channel_tag_dict[self.channels.index(profile_channel)][Tag_Type.pixel_area]
        xreal, yreal = self.channel_tag_dict[self.channels.index(profile_channel)][Tag_Type.scan_area]
        pixel_scaling = self.channel_tag_dict[self.channels.index(profile_channel)][Tag_Type.pixel_scaling]
        xvalues = [i*yreal/yres/pixel_scaling for i in range(yres*pixel_scaling)]
        # xvalues = np.linspace(0, 10, len(difference_profile))
        plt.plot(xvalues, difference_profile)
        plt.xlabel('Y [µm]')
        plt.ylabel('Phase difference')
        plt.ylim(ymin=0, ymax=2*np.pi)
        plt.title('Phase difference')
        plt.tight_layout()
        plt.show()

    def _Plot_Data_and_Profile_pos(self, channel, data, coordinates, orientation):
        if self.phase_indicator in channel:
            cmap = SNOM_phase
        elif self.amp_indicator in channel:
            cmap = SNOM_amplitude
        elif self.height_indicator in channel:
            cmap = SNOM_height
        fig, ax = plt.subplots()
        img = ax.pcolormesh(data, cmap=cmap)
        ax.invert_yaxis()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('phase', rotation=270)
        ax.legend()
        ax.axis('scaled')
        self._Plot_Profile_Lines(data, ax, coordinates, orientation)
        plt.title('You chose the following line profiles')
        plt.tight_layout()
        plt.show()

    def _Plot_Profile_Lines(self, data, ax, coordinates, orientation):
        xcoord = [coord[0] for coord in coordinates]
        ycoord = [coord[1] for coord in coordinates]
        if orientation == Definitions.vertical:
            ax.vlines(xcoord, ymin=0, ymax=len(data))
        elif orientation == Definitions.horizontal:
            ax.hlines(ycoord, xmin=0, xmax=len(data[0]))

    def _Get_Profiles_Coordinates(self, profile_channel, profiledata, preview_channel, previewdata, orientation, redo:bool=False, coordinates=None, redo_coordinates=None):
        if redo == False:
            coordinates = self._Get_Positions_from_Plot(preview_channel, previewdata)
        else:
            display_coordinates = [coordinates[i] for i in range(len(coordinates)) if i not in redo_coordinates]# remove coordinates to redo and plot the other ones while selecton is active
            redone_coordinates = self._Get_Positions_from_Plot(preview_channel, previewdata, display_coordinates, orientation)
            count = 0
            for index in redo_coordinates:
                coordinates[index] = redone_coordinates[count]
                count += 1

        self._Plot_Data_and_Profile_pos(profile_channel, profiledata, coordinates, orientation)
        print('Are you satisfied with the profile positions? Or would you like to change one ore more profile positions?')
        user_input_bool = self._User_Input_Bool() 
        if user_input_bool == False:
            user_input = self._User_Input('Please enter the indices of the profiles you like to redo, separated by a space character e.g. (0 1 3 11 ...)\nYour indices: ') 
            redo_coordinates = user_input.split(' ')
            redo_coordinates = [int(coord) for coord in redo_coordinates]
            print('coordinates to redo: ', redo_coordinates)
            print('Please select the new positons only for the indices you selected and in the same ordering, those were: ', redo_coordinates)
            coordinates = self._Get_Profiles_Coordinates(profile_channel, profiledata, preview_channel, previewdata, orientation, redo=True, coordinates=coordinates, redo_coordinates=redo_coordinates)
        
        return coordinates

    def Select_Profiles(self, profile_channel:str, preview_channel:str=None, orientation:Definitions=Definitions.vertical, width:int=10, coordinates:list=None):
        """This function lets the user select a profile with given width in pixels and displays the data.

        Args:
            profile_channel (str): channel to use for profile data extraction
            preview_channel (str, optional): channel to preview the profile positions. If not specified the height channel will be used for that. Defaults to None.
            orientation (Definitions, optional): profiles can be horizontal or vertical. Defaults to Definitions.vertical.
            width (int, optional): width of the profile in pixels, will calculate the mean. Defaults to 10.
            coordinates (list, optional): if you already now the position of your profile you can also specify the coordinates and skip the selection. Defaults to None.

        """
        if preview_channel is None:
            preview_channel = self.height_channel
        if preview_channel not in self.channels and profile_channel not in self.channels:
            print('The channels for preview and the profiles were not found in the memory, they will be loaded automatically.\nBe aware that all prior modifications will get deleted.')  
            self._Initialize_Data([profile_channel, preview_channel])#this will negate any modifications done prior like blurr...
        profiledata = self.all_data[self.channels.index(profile_channel)]
        previewdata = self.all_data[self.channels.index(preview_channel)]

        if coordinates == None:
            coordinates = self._Get_Profiles_Coordinates(profile_channel, profiledata, preview_channel, previewdata, orientation)
        
        print('The final profiles are shown in this plot.')
        self._Plot_Data_and_Profile_pos(profile_channel, profiledata, coordinates, orientation)
        # get the profile data and save to class variables
        # additional infos are also stored and can be used by plotting and analysis functions
        self.profiles = self._Get_Profile(profiledata, coordinates, orientation, width)
        self.profile_channel = profile_channel
        self.profile_orientation = orientation
        return self.profiles

    def Select_Profiles_SSH(self, profile_channel_amp:str, profile_channel_phase:str, preview_channel:str=None, orientation:Definitions=Definitions.vertical, width_amp:int=10, width_phase:int=1, coordinates:list=None):
        """This function lets the user select a profile with given width in pixels and displays the data.
        Specific function for ssh model measurements. This will create a plot of field per waveguide index for the topological array.
        The field is calculated from the amplitude profiles times the cosine of the phasedifference to the central waveguide. 

        Args:
            profile_channel_amp (str): amplitude channel for profile data
            profile_channel_phase (str): phase channel for profile data
            preview_channel (str, optional): channel to preview the profile positions. If not specified the height channel will be used for that. Defaults to None.
            orientation (Definitions, optional): profiles can be horizontal or vertical. Defaults to Definitions.vertical.
            width_amp (int, optional): width of the amplitude profile in pixels. Defaults to 10.
            width_phase (int, optional): width of the phase profile in pixels. Defaults to 1.
            coordinates (list, optional): if you already now the position of your profile you can also specify the coordinates and skip the selection. Defaults to None.
        """
        if preview_channel is None:
            preview_channel = self.height_channel
        if preview_channel not in self.channels or profile_channel_amp not in self.channels or profile_channel_phase not in self.channels:
            print('The channels for preview and the profiles were not found in the memory, they will be loaded automatically.\nBe aware that all prior modifications will get deleted.')  
            self._Initialize_Data([profile_channel_amp, profile_channel_phase, preview_channel])#this will negate any modifications done prior like blurr...
        profiledata_amp = self.all_data[self.channels.index(profile_channel_amp)]
        profiledata_phase = self.all_data[self.channels.index(profile_channel_phase)]
        previewdata = self.all_data[self.channels.index(preview_channel)]
        # get the profile coordinates
        if coordinates == None:
            coordinates = self._Get_Profiles_Coordinates(profile_channel_phase, profiledata_phase, preview_channel, previewdata, orientation)
        print(f'You selected the following coordinates: ', coordinates)
        print('The final profiles are shown in this plot.')
        self._Plot_Data_and_Profile_pos(profile_channel_phase, profiledata_phase, coordinates, orientation)
        self._Plot_Data_and_Profile_pos(profile_channel_amp, profiledata_amp, coordinates, orientation)
        self.profile_channel = profile_channel_phase
        self.profile_orientation = orientation

        # get the profile data for amp and phase
        phase_profiles = self._Get_Profile(profiledata_phase, coordinates, orientation, width_phase)
        # test:
        self._Display_Profile([phase_profiles[6], phase_profiles[16]])

        amp_profiles = self._Get_Profile(profiledata_amp, coordinates, orientation, width_amp)
        mean_amp = [np.mean(amp) for amp in amp_profiles]
        reference_index = int((len(phase_profiles)-1)/2)
        # phase_difference_profiles = [Phase_Analysis.Get_Profile_Difference(phase_profiles[reference_index], phase_profiles[i]) for i in range(len(phase_profiles))]
        flattened_profiles = [phase_analysis.Flatten_Phase_Profile(profile, +1) for profile in phase_profiles]
        self._Display_Profile(flattened_profiles, linestyle='-', title='Flattened phase profiles') # display the flattened profiles
        # phase_difference_profiles = [Phase_Analysis.Get_Profile_Difference_2(phase_profiles[reference_index], phase_profiles[i]) for i in range(len(phase_profiles))]
        phase_difference_profiles = [phase_analysis.Get_Profile_Difference_2(flattened_profiles[reference_index], flattened_profiles[i]) for i in range(len(flattened_profiles))]
        self._Display_Profile(phase_difference_profiles, linestyle='-', title='Phase difference to center wg') # display the phase difference profiles, no jumps close to 2 pi should occure or the average will lead to false values!
        # mean_phase_differences = [np.mean(diff) for diff in phase_difference_profiles]# todo this does not work!
        mean_phase_differences = [np.mean(diff) if np.mean(diff)>0 else np.mean(diff) + np.pi*2 for diff in phase_difference_profiles]# todo this does not work!
        real_per_wg_index = [mean_amp[i]*np.cos(mean_phase_differences[i]) for i in range(len(phase_profiles))]
        intensity_per_wg_index = [val**2 for val in real_per_wg_index]
        wg_indices = np.arange(-reference_index, reference_index+1)
        # print(wg_indices)
        fig = plt.figure(figsize=[4,2])
        plt.plot(wg_indices, real_per_wg_index, '-o', label='Real per wg index')
        plt.hlines(0, xmin=-10, xmax=10, linestyles='--')
        plt.ylabel(r'E$_z$ [arb.u]')
        plt.xlabel('Waveguide index')
        # plt.ylim([-0.04,0.04])
        
        plt.xticks(range(-reference_index, reference_index, 2))
        plt.legend()
        plt.tight_layout()
        plt.show()

        #same for intensity: hm not thought throu...
        # plt.plot(wg_indices, real_per_wg_index, '-o', label='Intensity per wg index')
        # plt.hlines(0, xmin=-10, xmax=10, linestyles='--')
        # plt.ylabel(r'I$_z$ [arb.u]')
        # plt.xlabel('Waveguide index')
        # # plt.ylim([-0.04,0.04])
        
        # plt.xticks(range(-reference_index, reference_index, 2))
        # plt.legend()
        # plt.tight_layout()
        # plt.show()
        
    def _Display_Profile(self, profiles, ylabel=None, labels=None, linestyle='x', title=None):
        if self.profile_orientation == Definitions.horizontal:
            xrange = self.channel_tag_dict[self.channels.index(self.profile_channel)][Tag_Type.scan_area][0]
            x_center_pos = self.channel_tag_dict[self.channels.index(self.profile_channel)][Tag_Type.center_pos][0]
            xres = self.channel_tag_dict[self.channels.index(self.profile_channel)][Tag_Type.pixel_area][0]# for now only profiles with lenght equal to scan dimensions are allowed
            xvalues = [x_center_pos - xrange/2 + x*(xrange/xres) for x in range(xres)]
            xlabel = 'X [µm]'
            if title == None:
                title = 'Horizontal profiles of channel ' + self.profile_channel
        elif self.profile_orientation == Definitions.vertical:
            yrange = self.channel_tag_dict[self.channels.index(self.profile_channel)][Tag_Type.scan_area][1]
            y_center_pos = self.channel_tag_dict[self.channels.index(self.profile_channel)][Tag_Type.center_pos][1]
            yres = self.channel_tag_dict[self.channels.index(self.profile_channel)][Tag_Type.pixel_area][1]# for now only profiles with lenght equal to scan dimensions are allowed
            xvalues = [y_center_pos - yrange/2 + y*(yrange/yres) for y in range(yres)]
            xlabel = 'Y [µm]'
            if title == None:
                title = 'Vertical profiles of channel ' + self.profile_channel
        # find out y label:
        if ylabel == None:
            if self.phase_indicator in self.profile_channel:
                ylabel = 'Phase'
            elif self.amp_indicator in self.profile_channel:
                ylabel = 'Amplitude [arb.u.]'
            elif self.height_indicator in self.profile_channel:
                ylabel = 'Height [nm]'
        for profile in profiles:
            index = profiles.index(profile)
            if labels == None:
                plt.plot(xvalues, profile, linestyle, label=f'Profile index: {index}')
            else:
                plt.plot(xvalues, profile, linestyle, label=labels[profiles.index(profile)])
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        # if labels != None:
        plt.legend()
        plt.tight_layout()
        plt.show()

    def Display_Profiles(self, ylabel:str=None, labels:list=None):
        """This function will display all current profiles from memory.

        Args:
            ylabel (str, optional): label of the y axis. The x axis label is in µm per default. Defaults to None.
            labels (list, optional): the description of the profiles. Will be displayed in the legend. Defaults to None.
        """
        self._Display_Profile(self.profiles)

    def Display_Flattened_Profile(self, phase_orientation:int):
        """This function will flatten all profiles in memory and display them. Only useful for phase profiles!

        Args:
            phase_orientation (int): direction of the phase, must be '1' or '-1'
        """
        flattened_profiles = [phase_analysis.Flatten_Phase_Profile(profile, phase_orientation) for profile in self.profiles]
        self._Display_Profile(flattened_profiles)

    def Display_Phase_Difference(self, reference_index:int):
        """This function will calculate the phase difference of all profiles relative to the profile specified by the reference index.

        Args:
            reference_index (int): index of the reference profile. Basically the nth-1 selected profile.
        """
        difference_profiles = [phase_analysis.Get_Profile_Difference(self.profiles[reference_index], self.profiles[i]) for i in range(len(self.profiles)) if i != reference_index]
        labels = ['Wg index ' + str(i) for i in range(len(difference_profiles))]
        self._Display_Profile(difference_profiles, 'Phase difference', labels)

    def _Get_Mean_Phase_Difference(self, profiles, reference_index:int):
        difference_profiles = [phase_analysis.Get_Profile_Difference(profiles[reference_index], profiles[i]) for i in range(len(profiles)) if i != reference_index]
        mean_differences = [np.mean(diff) for diff in difference_profiles]
        return mean_differences

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
        """This function scales the data such that each pixel is quadratic, eg. the physical dimensions are equal.
        This is important because the pixels will be set to quadratic in the plotting function.
        
        Args:
            channels [list]: list of channels the scaling should be applied to. If not specified the scaling will be applied to all channels
        """
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

    def Overlay_Forward_and_Backward_Channels(self, height_channel_forward:str, height_channel_backward:str, channels:list=None):
        """This function is ment to overlay the backwards and forwards version of the specified channels.
        You should only specify the forward version of the channels you want to overlay. The function will create a mean version
        which can then be displayed and saved. Note that the new version will be larger then the previous ones.

        Args:
            height_channel_forward (str): usual corrected height channel
            height_channel_backward (str): backwards height channel
            channels (list, optional): a list of all channels to be overlayed. Defaults to None.
        """
        all_channels = []
        for channel in channels:
            # print('extension: ', channel, self.backwards_indicator + channel)
            all_channels.extend([channel, self.backwards_indicator + channel])
        all_channels.extend([height_channel_forward, height_channel_backward])
        # print('specified channels to overlay: ', channels)
        # print('all channels:', all_channels)
        self._Initialize_Data(all_channels)
        # print('after initialisation' , self.channels)

        self.Set_Min_to_Zero([height_channel_forward, height_channel_backward])
        
        #scale and blurr channels for better overlap
        self.Scale_Channels()
        # self.Gauss_Filter_Channels()
        # self.Gauss_Filter_Channels_complex()

        height_data_forward = self.all_data[self.channels.index(height_channel_forward)]
        height_data_backward = self.all_data[self.channels.index(height_channel_backward)]
        
        #gauss blurr the data used for the alignment, so it might be a litte more precise
        height_channel_forward_blurr = self._Gauss_Blurr_Data(height_data_forward, 2)
        height_channel_backward_blurr = self._Gauss_Blurr_Data(height_data_backward, 2)

        # array_1 = height_data_forward[0]
        # array_2 = height_data_backward[0]

        '''
        mean_deviation_array = Realign.Calculate_Squared_Deviation(array_1, array_2)
        mean_deviation = np.mean(mean_deviation_array)
        x = range(len(array_1))
        plt.plot(x, array_1, label='array_2')
        plt.plot(x, array_2, label='array_1')
        plt.plot(x, mean_deviation_array, label="Mean deviation_array")
        plt.hlines(mean_deviation, label="Mean deviation", xmin=0, xmax=len(array_1))
        plt.legend()
        plt.show()
        '''

        # try to optimize by shifting second array and minimizing mean deviation
        pixel_scaling = self.channel_tag_dict[0][Tag_Type.pixel_scaling] # does not matter which channel to get the scaling from since all have been scaled
        N = 5*pixel_scaling #maximum iterations, scaled if pixelnumber was increased

        # Realign.Minimize_Deviation_1D(array_1, array_2, n_tries=N)
        # Realign.Minimize_Deviation_2D(height_data_forward, height_data_backward, n_tries=N)

        # get the index which minimized the deviation of the height channels
        # index = Realign.Minimize_Deviation_2D(height_data_forward, height_data_backward, N, False)
        index = realign.Minimize_Deviation_2D(height_channel_forward_blurr, height_channel_backward_blurr, N, False)
        # self.all_data[self.channels.index(height_channel_forward)], self.all_data[self.channels.index(height_channel_backward)] = Realign.Shift_Array_2D_by_Index(height_data_forward, height_data_backward, index)


        # print('trying to create mean data')
        # print(self.channels)
        for channel in channels:
            if self.backwards_indicator not in channel:
                #test:
                if self.height_indicator in channel:
                    # get current res and size and add the additional res and size due to addition of zeros while shifting
                    XRes, YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area]
                    XReal, YReal = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.scan_area]
                    XRes_new = XRes + abs(index)# absolute value? index can be negative, but resolution can only increase, same for real dimensions
                    XReal_new = XReal + XReal/XRes*abs(index)
                    
                    # create channel_dict for new mean data 
                    self.channel_tag_dict.append(self.channel_tag_dict[self.channels.index(channel)])
                    self.channel_tag_dict[-1][Tag_Type.pixel_area] = [XRes_new, YRes]
                    self.channel_tag_dict[-1][Tag_Type.scan_area] = [XReal_new, YReal]

                    # also create data dict entry
                    self.channels_label.append(self.channels_label[self.channels.index(channel)] + '_overlain')

                    # add new channel to channels
                    self.channels.append(channel + '_overlain')

                    #test realign (per scan) based on minimization of differences 
                    #not usable right now, drift compensation might lead to differently sized data
                    # self.all_data[self.channels.index(height_channel_forward)] = Realign.Minimize_Drift(self.all_data[self.channels.index(height_channel_forward)], display=False)
                    # self.all_data[self.channels.index(height_channel_backward)] = Realign.Minimize_Drift(self.all_data[self.channels.index(height_channel_backward)])

                    # shift the data of the forward and backwards channel to match
                    self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)] = realign.Shift_Array_2D_by_Index(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)], index)
        

                    # create mean data and append to all_data
                    self.all_data.append(realign.Create_Mean_Array(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)]))
                else:
                    # get current res and size and add the additional res and size due to addition of zeros while shifting
                    XRes, YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area]
                    XReal, YReal = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.scan_area]
                    XRes_new = XRes + abs(index)# absolute value? index can be negative, but resolution can only increase, same for real dimensions
                    XReal_new = XReal + XReal/XRes*abs(index)
                    
                    # create channel_dict for new mean data 
                    # print(self.channels)
                    # print(self.channel_tag_dict)
                    # print(self.channel_tag_dict[self.channels.index(channel)])
                    self.channel_tag_dict.append(self.channel_tag_dict[self.channels.index(channel)])
                    # print('old data dict: ', self.channel_tag_dict[-2])
                    # print('n#ew data dict: ', self.channel_tag_dict[-1])
                    # print('new data dict pixel area: ', self.channel_tag_dict[-1][Tag_Type.pixel_area])
                    self.channel_tag_dict[-1][Tag_Type.pixel_area] = [XRes_new, YRes]
                    self.channel_tag_dict[-1][Tag_Type.scan_area] = [XReal_new, YReal]

                    # also create data dict entry
                    self.channels_label.append(self.channels_label[self.channels.index(channel)] + '_overlain')

                    # add new channel to channels
                    self.channels.append(channel + '_overlain')
                    
                    #test realign (per scan) based on minimization of differences 
                    # self.all_data[self.channels.index(channel)] = Realign.Minimize_Drift(self.all_data[self.channels.index(channel)], display=False)
                    # self.all_data[self.channels.index(self.backwards_indicator+ channel)] = Realign.Minimize_Drift(self.all_data[self.channels.index(self.backwards_indicator+ channel)])

                    # shift the data of the forward and backwards channel to match
                    self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)] = realign.Shift_Array_2D_by_Index(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)], index)

                    # create mean data and append to all_data
                    self.all_data.append(realign.Create_Mean_Array(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)]))

                    # XRes, YRes = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.pixel_area]
                    # XReal, YReal = self.channel_tag_dict[self.channels.index(channel)][Tag_Type.scan_area]
            
    def Overlay_Forward_and_Backward_Channels_V2(self, height_channel_forward:str, height_channel_backward:str, channels:list=None):
        """
        Caution! This variant is ment to keep the scan size identical!

        This function is ment to overlay the backwards and forwards version of the specified channels.
        You should only specify the forward version of the channels you want to overlay. The function will create a mean version
        which can then be displayed and saved.

        Args:
            height_channel_forward (str): Usual corrected height channel
            height_channel_backward (str): Backwards height channel
            channels (list, optional): List of all channels to be overlayed. Only specify the forward direction. Defaults to None.
            If not specified only the amp channels and the height channel will be overlaid.
        """
        if channels is None:
            channels = [channel for channel in self.amp_channels if self.backwards_indicator not in channel]
            channels.append(self.height_channel)
        all_channels = []
        for channel in channels:
            all_channels.extend([channel, self.backwards_indicator + channel])
        if height_channel_forward not in channels:
            all_channels.extend([height_channel_forward, height_channel_backward])
        self._Initialize_Data(all_channels)
        self.Set_Min_to_Zero([height_channel_forward, height_channel_backward])
        
        #scale channels for more precise overlap
        self.Scale_Channels()
        height_data_forward = self.all_data[self.channels.index(height_channel_forward)]
        height_data_backward = self.all_data[self.channels.index(height_channel_backward)]
        
        #gauss blurr the data used for the alignment, so it might be a litte more precise
        height_channel_forward_blurr = self._Gauss_Blurr_Data(height_data_forward, 2)
        height_channel_backward_blurr = self._Gauss_Blurr_Data(height_data_backward, 2)

        # try to optimize by shifting second array and minimizing mean deviation
        pixel_scaling = self.channel_tag_dict[0][Tag_Type.pixel_scaling] # does not matter which channel to get the scaling from since all have been scaled
        N = 5*pixel_scaling #maximum iterations, scaled if pixelnumber was increased

        # get the index which minimized the deviation of the height channels
        index = realign.Minimize_Deviation_2D(height_channel_forward_blurr, height_channel_backward_blurr, N, False)

        for channel in channels:
            if self.backwards_indicator not in channel:
                if self.height_indicator in channel:
                    # create channel_dict for new mean data 
                    self.channel_tag_dict.append(self.channel_tag_dict[self.channels.index(channel)])

                    # also create data dict entry
                    self.channels_label.append(self.channels_label[self.channels.index(channel)] + '_overlain')

                    # add new channel to channels
                    self.channels.append(channel + '_overlain')
        
                    # create mean data and append to all_data
                    self.all_data.append(realign.Create_Mean_Array_V2(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)], index))
                else:
                    # create channel_dict for new mean data 
                    self.channel_tag_dict.append(self.channel_tag_dict[self.channels.index(channel)])

                    # also create data dict entry
                    self.channels_label.append(self.channels_label[self.channels.index(channel)] + '_overlain')

                    # add new channel to channels
                    self.channels.append(channel + '_overlain')
                    
                    # create mean data and append to all_data
                    self.all_data.append(realign.Create_Mean_Array_V2(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)], index))

    def Manually_Create_Complex_Channel(self, amp_channel:str, phase_channel:str, complex_type:str=None) -> None:
        """This function will manually create a realpart channel depending on the amp and phase channel you give.
        The channels don't have to be in memory. If they are not they will be loaded but not added to memory, only the realpart will be added.
        Carful, only for expert users!

        Args:
            amp_channel (str): Amplitude channel.
            phase_channel (str): Phase channel.
            complex_type (str, optional): Type of the data you want to create. 'real' creates the realpart, 'imag' the imaginary part.
                If not specified both will be created. Defaults to None.

        Returns:
            None
        """
        # check if channels match, check for data type (amp, phase) and demodulation order
        if self.amp_indicator not in amp_channel or self.phase_indicator not in phase_channel:
            print('The specified channels are not specified as needed!')
            exit()
        demodulation = amp_channel[1:2]
        if demodulation not in phase_channel:
            print('The channels you specified are not from the same demodulation order!\nProceeding anyways...')
        # check if channels are in memory, if not load the data
        if amp_channel not in self.channels:
            amp_data, amp_dict = self._Load_Data(amp_channel)
        else:
            amp_data = self.all_data[self.channels.index(amp_channel)]
            amp_dict = self.channel_tag_dict[self.channels.index(amp_channel)]
        if phase_channel not in self.channels:
            phase_data, phase_dict = self._Load_Data(phase_channel)
        else:
            phase_data = self.all_data[self.channels.index(phase_channel)]
            phase_dict = self.channel_tag_dict[self.channels.index(phase_channel)]
        # check if size is identical:
        xres_amp, yres_amp = amp_dict[Tag_Type.pixel_area]
        xres_phase, yres_phase = phase_dict[Tag_Type.pixel_area]
        if xres_amp != xres_phase or yres_amp != yres_phase:
            print('The data of the specified channels has different resolution!')
            exit()
        
        # create complex data:
        real_data = np.zeros((yres_amp, xres_amp))
        imag_data = np.zeros((yres_amp, xres_amp))
        for y in range(yres_amp):
            for x in range(xres_amp):
                real_data[y][x] = amp_data[y][x]*np.cos(phase_data[y][x])
                imag_data[y][x] = amp_data[y][x]*np.sin(phase_data[y][x])
        # create realpart and imaginary part channel and dict and add to memory
        real_channel = f'O{demodulation}' + self.real_indicator# + '_manipulated' # make shure not to overwrite the realpart created by the Synccorrection
        imag_channel = f'O{demodulation}' + self.imag_indicator# + '_manipulated' # make shure not to overwrite the imagpart created by the Synccorrection
        real_channel_dict = amp_dict
        imag_channel_dict = amp_dict

        if complex_type == 'real':
            self.channels.append(real_channel)
            self.all_data.append(real_data)
            self.channel_tag_dict.append(real_channel_dict)
            self.channels_label.append(real_channel)
        elif complex_type == 'imag':
            self.channels.append(imag_channel)
            self.all_data.append(imag_data)
            self.channel_tag_dict.append(imag_channel_dict)
            self.channels_label.append(imag_channel)
        elif complex_type is None:
            # just save both
            self.channels.append(real_channel)
            self.all_data.append(real_data)
            self.channel_tag_dict.append(real_channel_dict)
            self.channels_label.append(real_channel)

            self.channels.append(imag_channel)
            self.all_data.append(imag_data)
            self.channel_tag_dict.append(imag_channel_dict)
            self.channels_label.append(imag_channel)



# could be exported to external file

def Set_nan_to_zero(data) -> np.array:
    xres = len(data[0])
    yres = len(data)
    for y in range(yres):
        for x in range(xres):
            if str(data[y][x]) == 'nan':
                data[y][x] = 0
    return data       

# needed for the Realign function
def Gauss_Function(x, A, mu, sigma, offset):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + offset

def Get_Largest_Abs(val1, val2):
    if abs(val1) > abs(val2): return abs(val1)
    else: return abs(val2)
# not in use, delete?
'''
    def Center_Realpart(self) -> None:
        if self.file_type== File_Type.standard or self.file_type == File_Type.standard_new:
            for channel in self.channels:
                if ('O' in channel) and ('R' in channel):
                    data = self.all_data[self.channels.index(channel)]
                    
        else:
            print('Not yet implemented')
            pass

    def Create_Horizontal_Profile(self, channel:str, filepath:str=None):
        # self.Display_Channels()
        # print(f'1self.channels: {self.channels}')
        self._Initialize_Data([channel])
        # print(f'2self.channels: {self.channels}')
        # self.Display_Channels()
        # print(self.channel_tag_dict)

        if filepath == None:
            filepath = self.directory_name + '/' + self.filename + ' ' + channel + '_profile.txt'
        #create profile
        channel_index = self.channels.index(channel)
        # print(f'channel_index: {channel_index}')
        # print(f'trying to create a profile for channel: {channel}')
        # print(f'self.all_data[channel_index]: {self.all_data[channel_index]}')
        # print(f'self.all_data: {self.all_data}')
        # xres = len(self.all_data[channel_index][0])
        # yres = len(self.all_data[channel_index])
        # print(f'xres: {xres}')
        # print(f'yres: {yres}')
        profile = Profile.Horizontal_Profile(self.all_data[channel_index])
        
        xres, yres = self.channel_tag_dict[channel_index][Tag_Type.pixel_area]
        xreal, yreal = self.channel_tag_dict[channel_index][Tag_Type.scan_area]
        xcenter, ycenter = self.channel_tag_dict[channel_index][Tag_Type.center_pos]
        x_array = [xcenter-xreal/2 + xreal/xres*i for i in range(len(profile))]
        # print(f'xres , yres : {xres, yres}')
        # print(f'xres , yres : {xreal, yreal}')
        
        #create file
        file = open(filepath, 'w')
        if 'Z' in channel:
            file.write('#X[m]\tY[nm]\n')
            scaling = pow(10,-9)
        else: 
            file.write('#X[m]\tY[arb.u.]\n')
            scaling = 1
        for i in range(len(profile)):
            file.write(f'{round(x_array[i],9)}\t{round(profile[i]*scaling,5)}\n')
        file.close()

            
    def Test_Directionality(self, filename, channels:list=None):
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
        
        CC_instance = ChiralCoupler(height_data, self.all_data, self.channels)
        # the height data is used to fit, the result can be used for all following functions
        CC_instance.Integrate_Amplitude(amp_channel)
        # CC_instance.Display_Phase_Profile(phase_channel)
        # CC_instance.Display_Phase_Difference(phase_channel)
        # CC_instance.Export_Data(amp_channel, filename)
        CC_instance.Export_Mean_Data(filename)
        

        # Fit_Decaying_Amplitude()
        # ChiralCoupler(height_data, amp_data[0], amp_channels[0]).Integrate_Amplitude()
        # ChiralCoupler(height_data, phase_data[0], phase_channels[0]).Display_Phase_Profile()
        # ChiralCoupler(height_data, phase_data[0], phase_channels[0]).Display_Phase_Difference()        
            '''