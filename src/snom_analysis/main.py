'''This module contains the basic classes and functions for the snom analysis.'''

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

from scipy.ndimage import gaussian_filter # one could implement a bunch more filters
from scipy.optimize import curve_fit
from struct import unpack, pack
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_point_clicker import clicker# used for getting coordinates from images
from matplotlib_scalebar.scalebar import ScaleBar # used for creating scale bars
from matplotlib import patches # used for creating rectangles 
import numpy as np
from datetime import datetime
from pathlib import Path, PurePath
import os
import pickle as pkl # for saving and loading pickle files, the plot memory is saved in a pickle file
import gc # garbage collector to free memory
import json # for saving and loading json files like the plotting parameters, easy to view and edit by the user
import ast # for string to list, dict ... conversion
# for gif creation
import imageio # for creating/viewing gifs
from matplotlib.animation import FuncAnimation
# for old version
from PIL import Image
# for config file
from configparser import ConfigParser

# import own functionality
from .lib.snom_colormaps import SNOM_height, SNOM_amplitude, SNOM_phase, SNOM_realpart, all_colormaps
from .lib.phase_slider import get_phase_offset
from .lib.rectangle_selector import select_rectangle
from .lib.data_range_selector import select_data_range
from .lib import realign
from .lib import profile
from .lib import phase_analysis
from .lib.file_handling import get_parameter_values, find_index, convert_header_to_dict
from .lib.profile_selector import select_profile
# import additional functions
from .lib.additional_functions import set_nan_to_zero, gauss_function, get_largest_abs, calculate_colorbar_size, mean_index_array
# import definitions such as measurement and channel tags
from .lib.definitions import Definitions, MeasurementTags, ChannelTags, PlotDefinitions, MeasurementTypes
from .lib.height_masking import get_height_treshold
 
# new version is based on filehandler to do basic stuff and then a class for each different measurement type like snom/afm, approach curves, spectra etc.
class FileHandler(PlotDefinitions):
    """This class handles the measurement filetype and all toplevel functionality.
    This class will be called by each measurement type class to handle the filetype, measurement and channel dictionaries and the config file.
    When creating a new instance of this class the config file will be loaded and the filetype will be determined.
    Also the measurement tag dictionary will be created.
    
    Args:
        directory_name (str): The path of the directory where the measurement files are stored.
        title (str, optional): The title of the measurement. Defaults to None.
    """
    def __init__(self, directory_name:str, title:str=None) -> None:
        self.measurement_type = MeasurementTypes.NONE
        self.directory_name = Path(directory_name)
        self.filename = Path(PurePath(self.directory_name).parts[-1])
        self._generate_savefolder()
        self.measurement_title = title # If a measurement_title is specified it will precede the automatically created title based on the channel dictionary
        self.logfile_path = self._initialize_logfile()
        # testing the new config file:
        
        if self.config_path.exists():
            self._load_config() # load the config file
        else:
            print('Config file not found, creating a new one.')
            self._create_default_config() # create a default config file if not existing
        
        self._initialize_file_type()
        
    def _generate_savefolder(self):
        """Generate savefolder if not already existing. Careful, has to be the same one as for the snom plotter gui app.
        """
        # create parent folder in the user directory, both snom analysis and plotting cofig files will be saved there
        parent_folder = Path(os.path.expanduser('~')) / Path('SNOM_Config')
        if not Path.exists(parent_folder):
            os.makedirs(parent_folder)
        # create a save folder for the snom analysis config files
        self.save_folder = Path(os.path.expanduser('~')) / parent_folder / Path('SNOM_Analysis')
        if not Path.exists(self.save_folder):
            os.makedirs(self.save_folder)
        # define the paths for the different files
        self.all_subplots_path = self.save_folder / Path('all_subplots.p')
        self.plotting_parameters_path = self.save_folder / Path('plotting_parameters.json') # probably not a good idea to use the same folder as the snom plotter app
        self.config_path = self.save_folder / Path('config.ini')
        self.mpl_style_path = self.save_folder / Path('snom_analysis.mplstyle')
     
    def _initialize_file_type(self) -> None:
        # try to find the filetype automatically
        self._find_filetype() 

    def _create_default_config(self):
        """This function creates a default config file in case the script is run for the first time or the old config file is missing.
        This can also be called to reset the config file to default settings. But all manual changes will be lost.
        """
        config = ConfigParser()
        # the order is important, as the script will try to find the filetype in the order specified here
        # these are just the filetypes i have encountered so far, more can be added
        config['FILETYPES'] = {
            'filetype1': '<FILETYPE1>', # 1.10.9592.0 standard_new
            'filetype2': '<FILETYPE2>', # 1.8.5017.0 standard
            'filetype3': '<FILETYPE3>', # aachen ascii
            'filetype4': '<FILETYPE4>', # aachen gsf # not supported yet
            'filetype5': '<FILETYPE5>', # version 1.6.3359.1
            'filetype6': '<FILETYPE6>', # comsol
        }
        # old not needed anymore
        config['PARAMETERTYPES'] = {
            'PARAMETERTYPE1': 'html',
            'PARAMETERTYPE2': 'txt',
            'PARAMETERTYPE3': 'html_new',
            'PARAMETERTYPE4': 'html_neaspec_version_1_6_3359_1',
            'PARAMETERTYPE5': 'comsol_txt',
            'PARAMETERTYPE6': 'new_parameters_txt',
        }
        config['FILETYPE1'] = {
            'filetype': '<standard_new>',
            'parametertype': '<new_parameters_txt>',
            'phase_channels': ['O1P','O2P','O3P','O4P','O5P', 'R-O1P','R-O2P','R-O3P','R-O4P','R-O5P'],
            'amp_channels': ['O1A','O2A','O3A','O4A','O5A', 'R-O1A','R-O2A','R-O3A','R-O4A','R-O5A'],
            'real_channels': ['O1Re', 'O2Re', 'O3Re', 'O4Re', 'R-O5Re', 'R-O1Re', 'R-O2Re', 'R-O3Re', 'R-O4Re', 'R-O5Re'],
            'imag_channels': ['O1Im', 'O2Im', 'O3Im', 'O4Im', 'R-O5Im', 'R-O1Im', 'R-O2Im', 'R-O3Im', 'R-O4Im', 'R-O5Im'],
            'height_channel': '<Z C>',
            'height_channels': ['Z C', 'R-Z C'],
            'mechanical_channels': ['M0A', 'M0P', 'M1A', 'M1P', 'M2A', 'M2P', 'M3A', 'M3P', 'M4A', 'M4P', 'M5A', 'M5P', 'R-M0A', 'R-M0P', 'R-M1A', 'R-M1P', 'R-M2A', 'R-M2P', 'R-M3A', 'R-M3P', 'R-M4A', 'R-M4P', 'R-M5A', 'R-M5P'],
            'preview_ampchannel': '<O2A>',
            'preview_phasechannel': '<O2P>',
            'preview_channels': ['O2A', 'O2P', 'Z C'],
            'height_indicator': '<Z>',
            'amp_indicator': '<A>',
            'phase_indicator': '<P>',
            'backwards_indicator': '<R->',
            'real_indicator': '<Re>',
            'imag_indicator': '<Im>',
            'optical_indicator': '<O>',
            'mechanical_indicator': '<M>',
            'channel_prefix_default': '< >',
            'channel_prefix_custom': '< >',
            'channel_suffix_default': '< raw>',
            'channel_suffix_custom': '<>',
            'channel_suffix_synccorrected_phase': '<_synccorrected>',
            'channel_suffix_manipulated': '<_manipulated>',
            'channel_suffix_overlain': '<_overlain>',
            'parameters_name': '<.txt>', # measurement_directory + parameters_name
            'parameters_header_indicator': '<# >',
            'parameters_separator': '<:>',
            'file_ending': '<.gsf>',
            'phase_offset_default': np.pi, # shift raw data to the interval [0, 2pi]
            'phase_offset_custom': 0, # assume custom data is already in the interval [0, 2pi]
            'rounding_decimal_amp_default': 5,
            'rounding_decimal_amp_custom': 5,
            'rounding_decimal_phase_default': 5,
            'rounding_decimal_phase_custom': 5,
            'rounding_decimal_complex_default': 5,
            'rounding_decimal_complex_custom': 5,
            'rounding_decimal_height_default': 2, # when in nm
            'rounding_decimal_height_custom': 2, # when in nm
            'height_scaling_default': 10**9, # data is in m convert to nm
            'height_scaling_custom': 10**9, # data is in m convert to nm
            'measurement_tags': {
                # carful the keys will be used to reference enums, so they should be unique and uppercase, they also must be identical for all filetypes
                # the values are the tags in the parameter file so they should match the file format
                'SCAN': 'Scan', # scan type, afm, snom, approach curve, 2d/3d, PsHet...
                'PROJECT': 'Project',
                'DESCRIPTION': 'Description',
                'DATE': 'Date',
                'SCANNERCENTERPOSITION': 'Scanner Center Position (X, Y)',
                'ROTATION': 'Rotation',
                'SCANAREA': 'Scan Area (X, Y, Z)',
                'PIXELAREA': 'Pixel Area (X, Y, Z)',
                'AVERAGING': 'Averaging',
                'INTEGRATIONTIME': 'Integration time',
                'LASERSOURCE': 'Laser Source',
                'DETECTOR': 'Detector',
                'TARGETWAVELENGTH': 'Target Wavelength',
                'DEMODULATIONMODE': 'Demodulation Mode',
                'TIPFREQUENCY': 'Tip Frequency',
                'TIPAMPLITUTDE': 'Tip Amplitude',
                'TAPPINGAMPLITUDE': 'Tapping Amplitude',
                'MODULATIONFREQUENCY': 'Modulation Frequency',
                'MODULATIONAMPLITUDE': 'Modulation Amplitude',
                'MODULATIONOFFSET': 'Modulation Offset',
                'SETPOINT': 'Setpoint',
                'REGULATOR': 'Regulator (P, I, D)',
                'TIPPOTENTIAL': 'Tip Potential',
                'M1ASCALING': 'M1A Scaling',
                'Q-FACTOR': 'Q-Factor',
                'VERSION': 'Version',
            },
            'channel_tags': {
                'PIXELAREA': ['XRes', 'YRes'],
                'YINCOMPLETE': 'YResIncomplete',
                'ROTATION': 'Neaspec_Angle',
                'SCANAREA': ['XReal', 'YReal'],
                'SCANNERCENTERPOSITION': ['XOffset', 'YOffset'],
                'XYUNIT': 'XYUnits',
                'ZUNIT': 'ZUnits',
                'WAVENUMBERSCALING': 'Neaspec_WavenumberScaling',
            },
        }
        config['FILETYPE2'] = {
            'filetype': '<standard>',
            'parametertype': '<new_parameters_txt>',
            'phase_channels': ['O1P','O2P','O3P','O4P','O5P', 'R-O1P','R-O2P','R-O3P','R-O4P','R-O5P'],
            'amp_channels': ['O1A','O2A','O3A','O4A','O5A', 'R-O1A','R-O2A','R-O3A','R-O4A','R-O5A'],
            'real_channels': ['O1Re', 'O2Re', 'O3Re', 'O4Re', 'R-O5Re', 'R-O1Re', 'R-O2Re', 'R-O3Re', 'R-O4Re', 'R-O5Re'],
            'imag_channels': ['O1Im', 'O2Im', 'O3Im', 'O4Im', 'R-O5Im', 'R-O1Im', 'R-O2Im', 'R-O3Im', 'R-O4Im', 'R-O5Im'],
            'height_channel': '<Z C>',
            'height_channels': ['Z C', 'R-Z C'],
            'mechanical_channels': ['M0A', 'M0P', 'M1A', 'M1P', 'M2A', 'M2P', 'M3A', 'M3P', 'M4A', 'M4P', 'M5A', 'M5P', 'R-M0A', 'R-M0P', 'R-M1A', 'R-M1P', 'R-M2A', 'R-M2P', 'R-M3A', 'R-M3P', 'R-M4A', 'R-M4P', 'R-M5A', 'R-M5P'],
            'preview_ampchannel': '<O2A>',
            'preview_phasechannel': '<O2P>',
            'preview_channels': ['O2A', 'O2P', 'Z C'],
            'height_indicator': '<Z>',
            'amp_indicator': '<A>',
            'phase_indicator': '<P>',
            'backwards_indicator': '<R->',
            'real_indicator': '<Re>',
            'imag_indicator': '<Im>',
            'optical_indicator': '<O>',
            'mechanical_indicator': '<M>',
            'channel_prefix_default': '< >',
            'channel_prefix_custom': '< >',
            'channel_suffix_default': '<>',
            'channel_suffix_custom': '<>',
            'channel_suffix_synccorrected_phase': '<_synccorrected>',
            'channel_suffix_manipulated': '<_manipulated>',
            'channel_suffix_overlain': '<_overlain>',
            'parameters_name': '<.txt>', # measurement_directory + parameters_name
            'parameters_header_indicator': '<# >',
            'parameters_separator': '<:>',
            'file_ending': '<.gsf>',
            'phase_offset_default': np.pi, # shift raw data to the interval [0, 2pi]
            'phase_offset_custom': 0, # assume custom data is already in the interval [0, 2pi]
            'rounding_decimal_amp_default': 5,
            'rounding_decimal_amp_custom': 5,
            'rounding_decimal_phase_default': 5,
            'rounding_decimal_phase_custom': 5,
            'rounding_decimal_complex_default': 5,
            'rounding_decimal_complex_custom': 5,
            'rounding_decimal_height_default': 2, # when in nm
            'rounding_decimal_height_custom': 2, # when in nm
            'height_scaling_default': 10**9, # data is in m convert to nm
            'height_scaling_custom': 10**9, # data is in m convert to nm
            'measurement_tags': {
                # carful the keys will be used to create enums, so they should be unique and uppercase, they also must be identical for all filetypes
                # the values are the tags in the file so they should match the file format
                # 'SCAN': 'Scan', # scan type, afm, snom, approach curve, 2d/3d, PsHet...
                'PROJECT': 'Project',
                'DESCRIPTION': 'Description',
                'DATE': 'Date',
                'SCANNERCENTERPOSITION': 'Scanner Center Position (X, Y)',
                'ROTATION': 'Rotation',
                'SCANAREA': 'Scan Area (X, Y, Z)',
                'PIXELAREA': 'Pixel Area (X, Y, Z)',
                'AVERAGING': 'Averaging',
                'INTEGRATIONTIME': 'Integration time',
                'LASERSOURCE': 'Laser Source',
                # 'DETECTOR': 'Detector',
                'TARGETWAVELENGTH': 'Target Wavelength',
                'DEMODULATIONMODE': 'Demodulation Mode',
                'TIPFREQUENCY': 'Tip Frequency',
                'TIPAMPLITUTDE': 'Tip Amplitude',
                'TAPPINGAMPLITUDE': 'Tapping Amplitude',
                'MODULATIONFREQUENCY': 'Modulation Frequency',
                'MODULATIONAMPLITUDE': 'Modulation Amplitude',
                'MODULATIONOFFSET': 'Modulation Offset',
                'SETPOINT': 'Setpoint',
                'REGULATOR': 'Regulator (P, I, D)',
                'TIPPOTENTIAL': 'Tip Potential',
                'M1ASCALING': 'M1A Scaling',
                # 'Q-FACTOR': 'Q-Factor',
                'VERSION': 'Version',
            },
            'channel_tags': {
                'PIXELAREA': ['XRes', 'YRes'],
                'YINCOMPLETE': 'YResIncomplete',
                # 'ROTATION': 'Neaspec_Angle',
                'SCANAREA': ['XReal', 'YReal'],
                'SCANNERCENTERPOSITION': ['XOffset', 'YOffset'],
                'XYUNIT': 'XYUnits',
                'ZUNIT': 'ZUnits',
                'WAVENUMBERSCALING': 'Neaspec_WavenumberScaling',
            },
        }
        config['FILETYPE3'] = {
            'filetype': '<aachen_ascii>',
            'parametertype': '<new_parameters_txt>',
            'phase_channels': ['O1-F-arg','O2-F-arg','O3-F-arg','O4-F-arg', 'O1-B-arg','O2-B-arg','O3-B-arg','O4-B-arg'],
            'amp_channels': ['O1-F-abs','O2-F-abs','O3-F-abs','O4-F-abs', 'O1-B-abs','O2-B-abs','O3-B-abs','O4-B-abs'],
            'real_channels': ['O1-F-Re','O2-F-Re','O3-F-Re','O4-F-Re','O1-B-Re','O2-B-Re','O3-B-Re','O4-B-Re'],
            'imag_channels': ['O1-F-Im','O2-F-Im','O3-F-Im','O4-F-Im','O1-B-Im','O2-B-Im','O3-B-Im','O4-B-Im'],
            'height_channel': '<MT-F-abs>',
            'height_channels': ['MT-F-abs', 'MT-B-abs'],
            'mechanical_channels': ['M0-F-abs', 'M0-F-arg', 'M1-F-abs', 'M1-F-arg', 'M2-F-abs', 'M2-F-arg', 'M3-F-abs', 'M3-F-arg', 'M4-F-abs', 'M4-F-arg', 'M5-F-abs', 'M5-F-arg', 'M0-B-abs', 'M0-B-arg', 'M1-B-abs', 'M1-B-arg', 'M2-B-abs', 'M2-B-arg', 'M3-B-abs', 'M3-B-arg', 'M4-B-abs', 'M4-B-arg', 'M5-B-abs', 'M5-B-arg'],
            'preview_ampchannel': '<O2-F-abs>',
            'preview_phasechannel': '<O2-F-arg>',
            'preview_channels': ['O2-F-abs', 'O2-F-arg', 'MT-F-abs'],
            'height_indicator': '<MT>',
            'amp_indicator': '<abs>',
            'phase_indicator': '<arg>',
            'real_indicator': '<Re>',#not used
            'imag_indicator': '<Im>',#not used
            'optical_indicator': '<O>',
            'mechanical_indicator': '<M>',
            'backwards_indicator': '<-B->',
            'channel_prefix_default': '<_>',
            'channel_prefix_custom': '<_>',
            'channel_suffix_default': '<>',
            'channel_suffix_custom': '<>',
            'channel_suffix_synccorrected_phase': '<_synccorrected>',
            'channel_suffix_manipulated': '<_manipulated>',
            'channel_suffix_overlain': '<_overlain>',
            'parameters_name': '<.parameters.txt>', # measurement_directory + parameters_name
            'parameters_header_indicator': '<>',
            'parameters_separator': '<: >',
            'file_ending': '<.ascii>',
            # definitions for data loading:
            # todo the detector voltages should be handeled here, the following values are just placeholders
            # also gsf file reading for the gwyddion dump format is not implemented yet but ascii somewhat works
            'phase_offset_default': np.pi, # shift raw data to the interval [0, 2pi]
            'phase_offset_custom': 0, # assume custom data is already in the interval [0, 2pi]
            'rounding_decimal_amp_default': 5,
            'rounding_decimal_amp_custom': 5,
            'rounding_decimal_phase_default': 5,
            'rounding_decimal_phase_custom': 5,
            'rounding_decimal_complex_default': 5,
            'rounding_decimal_complex_custom': 5,
            'rounding_decimal_height_default': 2, # when in nm
            'rounding_decimal_height_custom': 2, # when in nm
            'height_scaling_default': 10**9, # data is in m convert to nm
            'height_scaling_custom': 10**9, # data is in m convert to nm
            'measurement_tags': {
                # carful the keys will be used to create enums, so they should be unique and uppercase, they also must be identical for all filetypes
                # the values are the tags in the file so they should match the file format
                'SCANAREA': ['scan_size_f (um)', 'scan_size_s (um)', 'scan_size_v (um)'],
                'PIXELAREA': ['resolution_f (pt)', 'resolution_s (pt)', 'resolution_v (pt)'],
                'INTEGRATIONTIME': 'pixel_time (ms)',
                'SCANNERCENTERPOSITION': ['offset_x (um)', 'offset_y (um)'],
                'ROTATION': 'rotation_a (deg)',
                'TIPFREQUENCY': 'probe_frequency (Hz)',
                'MODULATIONFREQUENCY': 'modulation_frequency (Hz)',
                'TAPPINGAMPLITUDE': 'probe_amplitude (V)',
                'MODULATIONAMPLITUDE': 'modulation_amplitude (V)',
                'MODULATIONOFFSET': 'modulation_offset (V)',
                'SETPOINT': 'setpoint (V)',
            },


        }
        # this filetype is not supported yet
        config['FILETYPE4'] = {
            'filetype': '<aachen_dumb>',
            'parametertype': '<new_parameters_txt>',
            'phase_channels': ['O1-F-arg','O2-F-arg','O3-F-arg','O4-F-arg', 'O1-B-arg','O2-B-arg','O3-B-arg','O4-B-arg'],
            'amp_channels': ['O1-F-abs','O2-F-abs','O3-F-abs','O4-F-abs', 'O1-B-abs','O2-B-abs','O3-B-abs','O4-B-abs'],
            'real_channels': ['O1-F-Re','O2-F-Re','O3-F-Re','O4-F-Re','O1-B-Re','O2-B-Re','O3-B-Re','O4-B-Re'],
            'imag_channels': ['O1-F-Im','O2-F-Im','O3-F-Im','O4-F-Im','O1-B-Im','O2-B-Im','O3-B-Im','O4-B-Im'],
            'mechanical_channels': ['M0-F-abs', 'M0-F-arg', 'M1-F-abs', 'M1-F-arg', 'M2-F-abs', 'M2-F-arg', 'M3-F-abs', 'M3-F-arg', 'M4-F-abs', 'M4-F-arg', 'M5-F-abs', 'M5-F-arg', 'M0-B-abs', 'M0-B-arg', 'M1-B-abs', 'M1-B-arg', 'M2-B-abs', 'M2-B-arg', 'M3-B-abs', 'M3-B-arg', 'M4-B-abs', 'M4-B-arg', 'M5-B-abs', 'M5-B-arg'],
            'height_channel': '<MT-F-abs>',
            'height_channels': ['MT-F-abs', 'MT-B-abs'],
            'preview_ampchannel': '<O2-F-abs>',
            'preview_phasechannel': '<O2-F-arg>',
            'preview_channels': ['O2-F-abs', 'O2-F-arg', 'MT-F-abs'],
            'height_indicator': '<MT>',
            'amp_indicator': '<abs>',
            'phase_indicator': '<arg>',
            'real_indicator': '<Re>',#not used
            'imag_indicator': '<Im>',#not used
            'optical_indicator': '<O>',
            'mechanical_indicator': '<M>',
            'backwards_indicator': '<-B->',
            'channel_prefix_default': '<_>',
            'channel_prefix_custom': '<_>',
            'channel_suffix_default': '<>',
            'channel_suffix_custom': '<>',
            'channel_suffix_synccorrected_phase': '<_synccorrected>',
            'channel_suffix_manipulated': '<_manipulated>',
            'channel_suffix_overlain': '<_overlain>',
            'parameters_name': '<.parameters.txt>', # measurement_directory + parameters_name
            'parameters_header_indicator': '<>',
            'parameters_separator': '<: >',
            'file_ending': '<.dumb>',
            # definitions for data loading:
            # todo the detector voltages should be handeled here, the following values are just placeholders
            # also gsf file reading for the gwyddion dump format is not implemented yet but ascii somewhat works
            'phase_offset_default': np.pi, # shift raw data to the interval [0, 2pi]
            'phase_offset_custom': 0, # assume custom data is already in the interval [0, 2pi]
            'rounding_decimal_amp_default': 5,
            'rounding_decimal_amp_custom': 5,
            'rounding_decimal_phase_default': 5,
            'rounding_decimal_phase_custom': 5,
            'rounding_decimal_complex_default': 5,
            'rounding_decimal_complex_custom': 5,
            'rounding_decimal_height_default': 2, # when in nm
            'rounding_decimal_height_custom': 2, # when in nm
            'height_scaling_default': 10**9, # data is in m convert to nm
            'height_scaling_custom': 10**9, # data is in m convert to nm
            'measurement_tags': {
                # carful the keys will be used to create enums, so they should be unique and uppercase, they also must be identical for all filetypes
                # the values are the tags in the file so they should match the file format
                'SCANAREA': ['scan_size_f (um)', 'scan_size_s (um)', 'scan_size_v (um)'],
                'PIXELAREA': ['resolution_f (pt)', 'resolution_s (pt)', 'resolution_v (pt)'],
                'INTEGRATIONTIME': 'pixel_time (ms)',
                'SCANNERCENTERPOSITION': ['offset_x (um)', 'offset_y (um)'],
                'ROTATION': 'rotation_a (deg)',
                'TIPFREQUENCY': 'probe_frequency (Hz)',
                'MODULATIONFREQUENCY': 'modulation_frequency (Hz)',
                'TAPPINGAMPLITUDE': 'probe_amplitude (V)',
                'MODULATIONAMPLITUDE': 'modulation_amplitude (V)',
                'MODULATIONOFFSET': 'modulation_offset (V)',
                'SETPOINT': 'setpoint (V)',
            },
        }
        
        config['FILETYPE5'] = {
            'filetype': '<standard>',
            'parametertype': '<new_parameters_txt>',
            'phase_channels': ['O1P','O2P','O3P','O4P','O5P', 'R-O1P','R-O2P','R-O3P','R-O4P','R-O5P'],
            'amp_channels': ['O1A','O2A','O3A','O4A','O5A', 'R-O1A','R-O2A','R-O3A','R-O4A','R-O5A'],
            'real_channels': ['O1Re', 'O2Re', 'O3Re', 'O4Re', 'R-O5Re', 'R-O1Re', 'R-O2Re', 'R-O3Re', 'R-O4Re', 'R-O5Re'],
            'imag_channels': ['O1Im', 'O2Im', 'O3Im', 'O4Im', 'R-O5Im', 'R-O1Im', 'R-O2Im', 'R-O3Im', 'R-O4Im', 'R-O5Im'],
            'height_channel': '<Z C>',
            'height_channels': ['Z C', 'R-Z C'],
            'mechanical_channels': ['M0A', 'M0P', 'M1A', 'M1P', 'M2A', 'M2P', 'M3A', 'M3P', 'M4A', 'M4P', 'M5A', 'M5P', 'R-M0A', 'R-M0P', 'R-M1A', 'R-M1P', 'R-M2A', 'R-M2P', 'R-M3A', 'R-M3P', 'R-M4A', 'R-M4P', 'R-M5A', 'R-M5P'],
            'preview_ampchannel': '<O2A>',
            'preview_phasechannel': '<O2P>',
            'preview_channels': ['O2A', 'O2P', 'Z C'],
            'height_indicator': '<Z>',
            'amp_indicator': '<A>',
            'phase_indicator': '<P>',
            'backwards_indicator': '<R->',
            'real_indicator': '<Re>',
            'imag_indicator': '<Im>',
            'optical_indicator': '<O>',
            'mechanical_indicator': '<M>',
            'channel_prefix_default': '< >',
            'channel_prefix_custom': '< >',
            'channel_suffix_default': '< raw>',
            'channel_suffix_custom': '<>',
            'channel_suffix_synccorrected_phase': '<_synccorrected>',
            'channel_suffix_manipulated': '<_manipulated>',
            'channel_suffix_overlain': '<_overlain>',
            'parameters_name': '<.txt>', # measurement_directory + parameters_name
            'parameters_header_indicator': '<>',
            'parameters_separator': '<>',
            'file_ending': '<.gsf>',
            'phase_offset_default': np.pi, # shift raw data to the interval [0, 2pi]
            'phase_offset_custom': 0, # assume custom data is already in the interval [0, 2pi]
            'rounding_decimal_amp_default': 5,
            'rounding_decimal_amp_custom': 5,
            'rounding_decimal_phase_default': 5,
            'rounding_decimal_phase_custom': 5,
            'rounding_decimal_complex_default': 5,
            'rounding_decimal_complex_custom': 5,
            'rounding_decimal_height_default': 2, # when in nm
            'rounding_decimal_height_custom': 2, # when in nm
            'height_scaling_default': 10**9, # data is in m convert to nm
            'height_scaling_custom': 10**9, # data is in m convert to nm
            'measurement_tags': {
                # carful the keys will be used to create enums, so they should be unique and uppercase, they also must be identical for all filetypes
                # the values are the tags in the file so they should match the file format
                # 'SCAN': 'Scan', # scan type, afm, snom, approach curve, 2d/3d, PsHet...
                'PROJECT': 'Project',
                'DESCRIPTION': 'Description',
                'DATE': 'Date',
                'SCANNERCENTERPOSITION': 'Scanner Center Position (X, Y)',
                'ROTATION': 'Rotation',
                'SCANAREA': 'Scan Size (X, Y, Z)',
                'PIXELAREA': 'Resolution (X, Y, Z)',
                'AVERAGING': 'Number of samples',
                'INTEGRATIONTIME': 'Pixel time',
                'LASERSOURCE': 'Laser Source',
                # 'DETECTOR': 'Detector',
                'TARGETWAVELENGTH': 'Target Wavelength',
                # 'DEMODULATIONMODE': 'Demodulation Mode',
                'TIPFREQUENCY': 'Tip Frequency',
                'TIPAMPLITUTDE': 'Tip Amplitude',
                'TAPPINGAMPLITUDE': 'Tapping Amplitude',
                'MODULATIONFREQUENCY': 'Modulation Frequency',
                'MODULATIONAMPLITUDE': 'Modulation Amplitude',
                'MODULATIONOFFSET': 'Modulation Offset',
                'SETPOINT': 'Setpoint',
                'REGULATOR': 'Regulator (P, I, D)',
                'TIPPOTENTIAL': 'Tip Potential',
                'M1ASCALING': 'M1A Scaling',
                # 'Q-FACTOR': 'Q-Factor',
                'VERSION': 'Version',
            },
            'channel_tags': {
                'PIXELAREA': ['XRes', 'YRes'],
                # 'YINCOMPLETE': 'YResIncomplete',
                # 'ROTATION': 'Neaspec_Angle',
                'SCANAREA': ['XReal', 'YReal'],
                'SCANNERCENTERPOSITION': ['XOffset', 'YOffset'],
                'XYUNIT': 'XYUnits',
                'ZUNIT': 'ZUnits',
                'WAVENUMBERSCALING': 'Neaspec_WavenumberScaling',
            },
        }
        config['FILETYPE6'] = {
            'filetype': '<comsol_gsf>',
            'parametertype': '<comsol_txt>',
            'all_channels_default': ['abs', 'arg', 'real', 'imag', 'Z'], # Z is not a standard channel, but the user might create it manually to show the simulation design
            'phase_channels': ['arg'],
            'amp_channels': ['abs'],
            'real_channels': ['real'],
            'imag_channels': ['imag'],
            'height_channel': '<Z>',
            'height_channels': ['Z'],
            'mechanical_channels': [],
            'preview_ampchannel': '<abs>',
            'preview_phasechannel': '<arg>',
            'preview_channels': ['abs', 'arg'],
            'height_indicator': '<Z>',
            'amp_indicator': '<abs>',
            'phase_indicator': '<arg>',
            'backwards_indicator': '<>',
            'real_indicator': '<real>',
            'imag_indicator': '<imag>',
            'optical_indicator': '<None>',
            'mechanical_indicator': '<None>',
            'channel_prefix_default': '<_>',
            'channel_prefix_custom': '<_>',
            'channel_suffix_default': '<>',
            'channel_suffix_custom': '<>',
            'channel_suffix_synccorrected_phase': '<_synccorrected>',
            'channel_suffix_manipulated': '<_manipulated>',
            'channel_suffix_overlain': '<_overlain>',
            'parameters_name': '<.txt>', # measurement_directory + parameters_name
            'parameters_header_indicator': '<# >',
            'parameters_separator': '<:>',
            'file_ending': '<.gsf>',

            # definitions for data loading:
            'phase_offset_default': 0, # assume default data is already in the interval [0, 2pi]
            'phase_offset_custom': 0, # assume custom data is already in the interval [0, 2pi]
            'rounding_decimal_amp_default': 5,
            'rounding_decimal_amp_custom': 5,
            'rounding_decimal_phase_default': 5,
            'rounding_decimal_phase_custom': 5,
            'rounding_decimal_complex_default': 5,
            'rounding_decimal_complex_custom': 5,
            'rounding_decimal_height_default': 2, # when in nm
            'rounding_decimal_height_custom': 2, # when in nm
            'height_scaling_default': 10**9, # data is in m convert to nm
            'height_scaling_custom': 10**9, # data is in m convert to nm
            'measurement_tags': {
                # carful the keys will be used to create enums, so they should be unique and uppercase, they also must be identical for all filetypes
                # the values are the tags in the file so they should match the file format
                # 'SCAN': 'Scan', # scan type, afm, snom, approach curve, 2d/3d, PsHet...
                # 'PROJECT': 'Project',
                # 'DESCRIPTION': 'Description',
                # 'DATE': 'Date',
                # 'SCANNERCENTERPOSITION': 'Scanner Center Position (X, Y)',
                # 'ROTATION': 'Rotation',
                'SCANAREA': 'Scan Area (X, Y)',
                'PIXELAREA': 'Pixel Area (X, Y)',
                # 'AVERAGING': 'Number of samples',
                # 'INTEGRATIONTIME': 'Pixel Time',
                # 'LASERSOURCE': 'Laser Source',
                # 'DETECTOR': 'Detector',
                # 'TARGETWAVELENGTH': 'Target Wavelength',
                # 'DEMODULATIONMODE': 'Demodulation Mode',
                # 'TIPFREQUENCY': 'Tip Frequency',
                # 'TIPAMPLITUTDE': 'Tip Amplitude',
                # 'TAPPINGAMPLITUDE': 'Tapping Amplitude',
                # 'MODULATIONFREQUENCY': 'Modulation Frequency',
                # 'MODULATIONAMPLITUDE': 'Modulation Amplitude',
                # 'MODULATIONOFFSET': 'Modulation Offset',
                # 'SETPOINT': 'Setpoint',
                # 'REGULATOR': 'Regulator (P, I, D)',
                # 'TIPPOTENTIAL': 'Tip Potential',
                # 'M1ASCALING': 'M1A Scaling',
                # 'Q-FACTOR': 'Q-Factor',
                'VERSION': 'Version',
            },
            'channel_tags': {
                'PIXELAREA': ['XRes', 'YRes'],
                # 'YINCOMPLETE': 'YResIncomplete',
                # 'ROTATION': 'Neaspec_Angle',
                'SCANAREA': ['XReal', 'YReal'],
                'SCANNERCENTERPOSITION': ['XOffset', 'YOffset'],
                'XYUNIT': 'XYUnits',
                # 'ZUNIT': 'ZUnits',
                # 'WAVENUMBERSCALING': 'Neaspec_WavenumberScaling',
            },
        }
        with open(self.config_path, 'w') as configfile:
            config.write(configfile)
        self.config = config

    def _load_config(self):
        """This function loads the config file and makes the config available throu self.config.
        """
        self.config = ConfigParser()
        with open(self.config_path, 'r') as f:
            self.config.read_file(f)

    def _load_mpl_style(self):
        if not Path.exists(self.mpl_style_path):
            # generate default mpl style file
            with open(self.mpl_style_path, 'w') as f:
                f.write('axes.grid: False\n')
                f.write('axes.grid.axis: both\n')
                f.write('axes.grid.which: major\n')
                f.write('grid.linestyle: -\n')
                f.write('grid.linewidth: 0.5\n')
                f.write('grid.color: black\n')
                f.write('xtick.direction: in\n')
                f.write('ytick.direction: in\n')
                f.write('xtick.minor.visible: True\n')
                f.write('ytick.minor.visible: True\n')
                f.write('xtick.major.size: 5\n')
                f.write('ytick.major.size: 5\n')
                f.write('xtick.minor.size: 3\n')
                f.write('ytick.minor.size: 3\n')
                f.write('xtick.major.width: 0.5\n')
                f.write('ytick.major.width: 0.5\n')
                f.write('xtick.minor.width: 0.5\n')
                f.write('ytick.minor.width: 0.5\n')
                f.write('xtick.major.pad: 5\n')
                f.write('ytick.major.pad: 5\n')
                f.write('xtick.minor.pad: 5\n')
                f.write('ytick.minor.pad: 5\n')
                f.write('xtick.major.top: True\n')
                f.write('ytick.major.right: True\n')
                f.write('xtick.minor.top: True\n')
                f.write('ytick.minor.right: True\n')
                f.write('axes.linewidth: 0.5\n')
                f.write('axes.edgecolor: black\n')
                f.write('axes.labelcolor: black\n')
                f.write('axes.labelsize: 12\n')
                f.write('axes.labelweight: normal\n')
                f.write('axes.labelpad: 4.0\n')
                f.write('axes.formatter.limits: -7, 7\n')
                f.write('axes.formatter.use_locale: False\n')
                f.write('axes.formatter.use_mathtext: False\n')
                f.write('axes.formatter.useoffset: True\n')
                f.write('axes.formatter.offset_threshold: 4\n')
                f.write('axes.formatter.min_exponent: 0\n')
        plt.style.use(self.mpl_style_path)

    def print_config(self):
        """This function prints the config file.
        """
        for section in self.config.sections():
            print(section)
            for option in self.config.options(section):
                print(f'{option} = {self.config.get(section, option)}')

    def _change_config(self, section:str, option:str, value):
        """This function changes the config file.

        Args:
            section (str): The section in the config file. Corresponds to the filetype, e.g. 'FILETYPE1'.
            option (str): The option in the section, e.g. amplitude_channels.
            value: The value to change, could be a string, int, float, list, dict, bool.
        """
        # if value is a string add quotes
        if isinstance(value, str):
            value = f'<{value}>'
        try:
            self.config[section][option] = value
        except:
            print('The specified section or option does not exist in the config file!')
            try:
                print('The available options are: ', self.config.options(section))
            except:	
                print('The available sections are: ', self.config.sections())
        # update the config file        
        with open(self.config_path, 'w') as configfile:
            self.config.write(configfile)
    
    def _get_from_config(self, option:str=None, section:str=None):
        """This function gets the value of an option in a section of the config file.
        If no option is specified the whole section is returned."""
        if section is None:
            # set the section to the file type if it is not specified, but only if file_type is defined
            try: section = self.file_type
            except: print('Filetype unknown, please specify the section! (In _get_from_config)')
        if option is None:
            return dict(self.config[section])
        else:
            value = self.config[section][option]
            # replace < and > with empty string if value is a string
            if isinstance(value, str):
                if value[0] == '<':
                    value = value.replace('<', '').replace('>', '')
                # convert string to list if it is a list
                # elif value[0] == '[':
                #     value = ast.literal_eval(value)
                # # convert string to dictionary if it is a dictionary
                # elif value[0] == '{':
                #     value = ast.literal_eval(value)
                # replace string with boolean if it is a boolean
                if value == 'True':
                    value = True
                elif value == 'False':
                    value = False
                elif value == 'None':
                    value = None
                else:
                    # try to convert string to float or int or list or dict
                    try:
                        value = ast.literal_eval(value)
                    except:
                        pass

            return value

    def _print_measurement_tags(self):
        """This function prints the measurement tags.
        """
        # print the content of the measurement tags class
        print('All measurement tags: ', list(MeasurementTags))

    def _find_filetype(self) -> bool:
        """This function tries to find the correct filetype for the given file.
        It will iterate through all filetypes in the config file and try to create the measurement tag dict.
        If the filetype is found the function returns True otherwise False.
        """
        filetypes = self._get_from_config(section='FILETYPES')
        for key in filetypes:
            filetype = self._get_from_config(key, 'FILETYPES')
            parameters_name = self._get_from_config('parameters_name', filetype)
            parameters_path = self.directory_name / Path(self.filename.name + parameters_name)
            # try to create the measurement tag dict
            succsess = self._create_measurement_tag_dict(parameters_path, filetype)
            # if succsess:
            #     print('measurement tag dict: ', self.measurement_tag_dict)
            # the correct creation of teh measurement tag dict is not enough to determine the filetype
            # try to also to create the channel tag dict for one arbitrary channel
            # self._initialize_file_type()
            self.file_type = filetype
            self._initialize_measurement_channel_indicators()
            # amp_channel = self._get_from_config('amp_channels', filetype)[0]
            # try to create the channel tag dict, if it fails the filetype is not correct
            # print('Trying to create channel tag dict')
            # print('all_channels_default[0]: ', self.all_channels_default[0])
            # print('filetype: ', filetype)
            # print('succsess: ', succsess)
            # this approach does not work for comsol files, approach curves and 3d scans
            # print('measurement_type: ', self.measurement_type)
            # in case the Filehandler was called directly the measurement type is not set yet
            # try to find the measurement type
            if self.measurement_type == MeasurementTypes.NONE:
                self._find_measurement_type()
            if self.measurement_type == MeasurementTypes.SNOM:
                try: self._create_channel_tag_dict([self.all_channels_default[0]])
                except: 
                    succsess = False
            self.file_type = None
            if succsess:
                # the correct filetype has been found
                # print(f'Filetype found: {filetype}')
                self.file_type = filetype
                # print('parameter dict was created successfully')                
                return True

        # if no filetype was found return False
        print('No filetype was found!')
        exit()
        return False

    def _find_measurement_type(self) -> None:
        # print('Trying to find the measurement type')
        if self.file_type != None:
            try:
                # not every filetype has a scan type
                scan_type = self.measurement_tag_dict[MeasurementTags.SCAN]
            except:
                # scan_type = None
                # self.plotting_mode = MeasurementTypes.NONE
                # todo, not all filetypes have a scan type, use additional ways to identify the measurement type
                # for now assume, that all files without a scan type are standard snom measurements
                self.measurement_type = MeasurementTypes.SNOM
            else:
                if 'Approach Curve' in scan_type:
                    self.measurement_type = MeasurementTypes.APPROACHCURVE
                elif '3D' in scan_type:
                    self.measurement_type = MeasurementTypes.SCAN3D
                elif 'Spectrum' in scan_type: # todo, not implemented yet
                    self.measurement_type = MeasurementTypes.SPECTRUM
                else:
                    self.measurement_type = MeasurementTypes.SNOM
        else:
            print('Could not identify the measurement type!')
            self.measurement_type = MeasurementTypes.NONE
        # print('Measurement type: ', self.measurement_type)

    def _initialize_logfile(self) -> str:
        # logfile_path = self.directory_name + '/python_manipulation_log.txt'
        logfile_path = self.directory_name / Path('python_manipulation_log.txt')
        file = open(logfile_path, 'a') # the new logdata will be appended to the existing file
        now = datetime.now()
        current_datetime = now.strftime("%d/%m/%Y %H:%M:%S")
        file.write(current_datetime + '\n' + 'filename = ' + self.filename.name + '\n')
        file.close()
        return logfile_path

    def _write_to_logfile(self, parameter_name:str, parameter):
        file = open(self.logfile_path, 'a')
        file.write(f'{parameter_name} = {parameter}\n')
        file.close()
 
    def _create_measurement_tag_dict(self, parameters_path:Path, filetype:str) -> bool:
        """This function creates a dictionary containing the measurement tags. The tags are extracted from the parameters file.
        If the tag dict cannot be created the function will return False otherwise True.

        Args:
            parameters_path (Path): The path to the parameters file.
            filetype (str): The filetype to use.
        """
        # first check if the file exists
        # print('trying to load parameters')
        # print('filetype: ', filetype)
        try:
            with open(parameters_path, 'r') as file:
                pass
        except:
            return False
        separator = self._get_from_config('parameters_separator', filetype)
        header_indicator = self._get_from_config('parameters_header_indicator', filetype)
        measurement_tags = self._get_from_config('measurement_tags', filetype)
        tags_list = list(measurement_tags.values())
        # print('tags_list: ', tags_list)
        # if tags contains a list of values flatten the list
        flattened_tags_list = []
        list_items = [] # keep track of list items to reverse the flattening after the creation of the parameters dict
        for i in range(len(tags_list)):
            tag = tags_list[i]
            if isinstance(tag, list):
                for item in tag:
                    flattened_tags_list.append(item)
            else:
                flattened_tags_list.append(tag)
            list_items.append(i)
        # if any(isinstance(i, list) for i in tags_list):
        #     tags_list = [item for sublist in tags_list for item in sublist]
        # print('flattened tags_list: ', tags_list)
        # print('trying to convert header to dict')  
        # print('flattenend_tags_list: ', flattened_tags_list)
        parameters_dict = convert_header_to_dict(parameters_path, separator=separator, header_indicator=header_indicator, tags_list=flattened_tags_list)
        # print('parameters_dict: ', parameters_dict)
        if parameters_dict is None:
            return False
        # reverse the flattening of the tags list and translate file tags to measurement tags
        new_parameters_dict = {}
        for i in range(len(tags_list)):
            tag = tags_list[i]
            if isinstance(tag, list):
                val_list = []
                count = 0
                for item in tag:
                    val_list.append(parameters_dict[tag[count]])
                    count += 1
                measurement_tag = list(measurement_tags.keys())[list(measurement_tags.values()).index(tag)]
                new_parameters_dict[measurement_tag] = val_list
                    # flattened_tags_list.append(item)
            else:
                val = parameters_dict[tag]
                measurement_tag = list(measurement_tags.keys())[list(measurement_tags.values()).index(tag)]
                new_parameters_dict[measurement_tag] = val
                # flattened_tags_list.append(tag)
            # list_items.append(i)

        # print('parameters_dict: ', parameters_dict)
        # print('new_parameters_dict: ', new_parameters_dict)
        # now create the measurement tag dict
        self.measurement_tag_dict = {}

        '''
        # SCAN = auto()   # scan type, afm, snom, approach curve, 2d/3d, PsHet...
        # PROJECT = auto()   
        # DESCRIPTION = auto()   
        # DATE = auto()   
        # SCANNERCENTERPOSITION = auto()   
        # ROTATION = auto()   
        # SCANAREA = auto()  
        # PIXELAREA = auto()   
        # AVERAGING = auto()   
        # INTEGRATIONTIME = auto()   
        # LASERSOURCE = auto()   
        # DETECTOR = auto()   
        # TARGETWAVELENGTH = auto()    
        # DEMODULATIONMODE = auto()   
        # TIPFREQUENCY = auto()   
        # TIPAMPLITUTDE = auto()   
        # TAPPINGAMPLITUDE = auto()   
        # MODULATIONFREQUENCY = auto()   
        # MODULATIONAMPLITUDE = auto()   
        # MODULATIONOFFSET = auto()   
        # SETPOINT = auto()   
        # REGULATOR = auto()   
        # TIPPOTENTIAL = auto()   
        # M1ASCALING = auto()  
        # QFACTOR = auto()   
        # VERSION = auto()
        # '''
        for key, value in new_parameters_dict.items():
            is_unit = False
            is_list = False
            if isinstance(value, list):
                is_list = True
                value = [item.replace(',', '') for item in value]
                # check if first value is a is_unit
                try: float(value[0])
                except: is_unit = True
                else: is_unit = False
            else: # sometimes only the is_unit is given
                try: float(value)
                except: is_unit = True
                else: is_unit = False
            if value == '':
                continue
            # if key in measurement_tags.values():
                # tag_key = list(measurement_tags.keys())[list(measurement_tags.values()).index(key)]
            # else:
            #     continue
            tag_key = key
            # print(f'tag_key: <{tag_key}>, value: <{value}>')
            # print(f'is_unit: {is_unit}, is_list: {is_list}')
            if tag_key == 'SCAN':
                self.measurement_tag_dict[MeasurementTags.SCAN] = value
            elif tag_key == 'PROJECT':
                self.measurement_tag_dict[MeasurementTags.PROJECT] = value
            elif tag_key == 'DESCRIPTION':
                self.measurement_tag_dict[MeasurementTags.DESCRIPTION] = value
            elif tag_key == 'DATE':
                self.measurement_tag_dict[MeasurementTags.DATE] = value
            elif tag_key == 'SCANNERCENTERPOSITION': # is_unit, x, y
                if is_unit:
                    try: self.measurement_tag_dict[MeasurementTags.SCANNERCENTERPOSITION] = [value[0], float(value[1]), float(value[2])]
                    except: self.measurement_tag_dict[MeasurementTags.SCANNERCENTERPOSITION] = [float(value[0]), float(value[1])]
                else:
                    self.measurement_tag_dict[MeasurementTags.SCANNERCENTERPOSITION] = [float(value[0]), float(value[1])]
            elif tag_key == 'ROTATION': # is_unit, angle
                if is_unit: self.measurement_tag_dict[MeasurementTags.ROTATION] = [value[0], float(value[1])]
                else: self.measurement_tag_dict[MeasurementTags.ROTATION] = float(value)
            elif tag_key == 'SCANAREA': # is_unit, x, y, z
                # some files have only 2 values for the scan area with or without is_unit
                # check if first value is a is_unit
                if is_unit:
                    try: self.measurement_tag_dict[MeasurementTags.SCANAREA] = [value[0], float(value[1]), float(value[2]), float(value[3])]
                    except: self.measurement_tag_dict[MeasurementTags.SCANAREA] = [value[0], float(value[1]), float(value[2])]
                else:
                    try: self.measurement_tag_dict[MeasurementTags.SCANAREA] = [float(value[0]), float(value[1]), float(value[2])]
                    except: self.measurement_tag_dict[MeasurementTags.SCANAREA] = [float(value[0]), float(value[1])]
            elif tag_key == 'PIXELAREA': # is_unit, x, y, z
                # print('PixelArea value: ', value)
                if is_unit:
                    try: self.measurement_tag_dict[MeasurementTags.PIXELAREA] = [value[0], int(value[1]), int(value[2]), int(value[3])]
                    except: self.measurement_tag_dict[MeasurementTags.PIXELAREA] = [int(value[0]), int(value[1]), int(value[2])]
                else:
                    try: self.measurement_tag_dict[MeasurementTags.PIXELAREA] = [int(value[0]), int(value[1]), int(value[2])]
                    except: self.measurement_tag_dict[MeasurementTags.PIXELAREA] = [int(value[0]), int(value[1])]
            elif tag_key == 'AVERAGING':
                self.measurement_tag_dict[MeasurementTags.AVERAGING] = int(value)  
            elif tag_key == 'INTEGRATIONTIME':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.INTEGRATIONTIME] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.INTEGRATIONTIME] = value
                    else: self.measurement_tag_dict[MeasurementTags.INTEGRATIONTIME] = float(value)
            elif tag_key == 'LASERSOURCE':
                self.measurement_tag_dict[MeasurementTags.LASERSOURCE] = value
            elif tag_key == 'DETECTOR':
                self.measurement_tag_dict[MeasurementTags.DETECTOR] = value
            elif tag_key == 'TARGETWAVELENGTH': # wavelength is usually not given...
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.TARGETWAVELENGTH] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.TARGETWAVELENGTH] = value
                    else: self.measurement_tag_dict[MeasurementTags.TARGETWAVELENGTH] = float(value)
            elif tag_key == 'DEMODULATIONMODE':
                self.measurement_tag_dict[MeasurementTags.DEMODULATIONMODE] = value
            elif tag_key == 'TIPFREQUENCY':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.TIPFREQUENCY] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.TIPFREQUENCY] = value
                    else: self.measurement_tag_dict[MeasurementTags.TIPFREQUENCY] = float(value)
            elif tag_key == 'TIPAMPLITUTDE':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.TIPAMPLITUTDE] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.TIPAMPLITUTDE] = value
                    else: self.measurement_tag_dict[MeasurementTags.TIPAMPLITUTDE] = float(value)
            elif tag_key == 'TAPPINGAMPLITUDE':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.TAPPINGAMPLITUDE] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.TAPPINGAMPLITUDE] = value
                    else: self.measurement_tag_dict[MeasurementTags.TAPPINGAMPLITUDE] = float(value)
            elif tag_key == 'MODULATIONFREQUENCY':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.MODULATIONFREQUENCY] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.MODULATIONFREQUENCY] = value
                    else: self.measurement_tag_dict[MeasurementTags.MODULATIONFREQUENCY] = float(value)
            elif tag_key == 'MODULATIONAMPLITUDE':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.MODULATIONAMPLITUDE] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.MODULATIONAMPLITUDE] = value
                    else: self.measurement_tag_dict[MeasurementTags.MODULATIONAMPLITUDE] = float(value)
            elif tag_key == 'MODULATIONOFFSET':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.MODULATIONOFFSET] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.MODULATIONOFFSET] = value
                    else: self.measurement_tag_dict[MeasurementTags.MODULATIONOFFSET] = float(value)
            elif tag_key == 'SETPOINT':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.SETPOINT] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.SETPOINT] = value
                    else: self.measurement_tag_dict[MeasurementTags.SETPOINT] = float(value)
            elif tag_key == 'REGULATOR':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.REGULATOR] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.REGULATOR] = value
                    else: self.measurement_tag_dict[MeasurementTags.REGULATOR] = float(value)
            elif tag_key == 'TIPPOTENTIAL':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.TIPPOTENTIAL] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.TIPPOTENTIAL] = value
                    else: self.measurement_tag_dict[MeasurementTags.TIPPOTENTIAL] = float(value)
            elif tag_key == 'M1ASCALING':
                if is_list:
                    self.measurement_tag_dict[MeasurementTags.M1ASCALING] = [value[0], float(value[1])]
                else:
                    if is_unit: self.measurement_tag_dict[MeasurementTags.M1ASCALING] = value
                    else: self.measurement_tag_dict[MeasurementTags.M1ASCALING] = float(value)
            elif tag_key == 'QFACTOR':
                self.measurement_tag_dict[MeasurementTags.QFACTOR] = float(value)
            elif tag_key == 'VERSION':
                self.measurement_tag_dict[MeasurementTags.VERSION] = value

        # only used by synccorrection, every other function should use the channels tag dict version, as pixel resolution could vary
        pixelarea = self._get_measurement_tag_dict_value(MeasurementTags.PIXELAREA)
        # scanarea = self.measurement_tag_dict[MeasurementTags.SCANAREA]
        scanarea = self._get_measurement_tag_dict_value(MeasurementTags.SCANAREA)
        if len(pixelarea) >= 3 and isinstance(pixelarea[0], str):
            self.XRes, self.YRes = pixelarea[1], pixelarea[2] # [unit, x, y, (z)]
            self.XReal, self.YReal = scanarea[1], scanarea[2]
        else:
            # print('pixelarea:', pixelarea)
            self.XRes, self.YRes = pixelarea[0], pixelarea[1] # [x, y]
            self.XReal, self.YReal = scanarea[0], scanarea[1]

        # if everything went well return True
        return True
        
    def print_measurement_tag_dict(self):
        """This function prints the measurement tag dict.
        """
        print('-------------------------------')
        print('Measurement tag dict:')
        print('-------------------------------')
        for key, value in self.measurement_tag_dict.items():
            print(f'{key} = {value}')

    def print_channel_tag_dict(self, channel=None):
        """This function prints the channel tag dict.

        Args:
            channel (str, optional): The channel to print. If None all channels will be printed. Defaults to None.
        """
        if channel is not None:
            print('-------------------------------')
            print(f'{channel} channel tag dict:')
            print('-------------------------------')
            for key, value in self.channel_tag_dict[self.channels.index(channel)].items():
                print(f'{key} = {value}')
        else:
            for channel in self.channels:
                print('-------------------------------')
                print(f'{channel} channel tag dict:')
                print('-------------------------------')
                for key, value in self.channel_tag_dict[self.channels.index(channel)].items():
                    print(f'{key} = {value}')

    def _replace_plotting_parameter_placeholders(self, dictionary:dict, placeholders:dict) -> dict:
        """This function replaces the placeholders in the plotting parameters dictionary with the actual values. 
        Afterwards it replaces the colormap placeholders with the actual colormaps.

        Args:
            dictionary (dict): plotting parameters dictionary
            placeholders (dict): dictionary containing the string definition of the placeholder and its value

        Returns:
            dict: the updated plotting parameters dictionary
        """
        # colormaps = {"<SNOM_amplitude>": SNOM_amplitude,
        #             "<SNOM_height>": SNOM_height,
        #             "<SNOM_phase>": SNOM_phase,
        #             "<SNOM_realpart>": SNOM_realpart}
        
        # first iterate through all placeholders and replace them in the dictionary
        for placeholder in placeholders:
            value = placeholders[placeholder]
            for key in dictionary:
                if placeholder in dictionary[key]:
                    dictionary[key] = dictionary[key].replace(placeholder, value)
                    # print('replaced channel!')
        # replace colormaps
        for key in dictionary:
            for colormap in all_colormaps:
                if colormap in dictionary[key]:
                    dictionary[key] = all_colormaps[colormap]
                    break
        return dictionary

    def _get_plotting_parameters(self) -> dict:
        """This will load the plotting parameters dictionary from the plotting_parameters.json file. If the file does not exist, it will be created with default values.
        The dictionary contains definitions for the colormaps, the colormap labels and the titles of the subplots. It also contains placeholders, which can be replaced by the actual values.
        The user can change the values in the plotting_parameters.json file to customize the plotting.

        Returns:
            dict: plotting parameters dictionary
        """
        try:
            with open(self.plotting_parameters_path, 'r') as file:
                plotting_parameters = json.load(file)
        except:
            self._generate_default_plotting_parameters()
            with open(self.plotting_parameters_path, 'r') as file:
                plotting_parameters = json.load(file)
        return plotting_parameters
    
    def _generate_default_plotting_parameters(self):
        dictionary = {
            "amplitude_cmap": "<SNOM_amplitude>",
            "amplitude_cbar_label": "Amplitude (arb.u.)",
            "amplitude_title": "<channel>",
            "phase_cmap": "<SNOM_phase>",
            "phase_cbar_label": "Phase (rad)",
            "phase_title": "<channel>",
            "phase_positive_title": "Positively corrected phase <channel>",
            "phase_negative_title": "Negatively corrected phase <channel>",
            "height_cmap": "<SNOM_height>",
            "height_cbar_label": "Height (nm)",
            "height_title": "<channel>",
            "real_cmap": "<SNOM_realpart>",
            "real_cbar_label": "E (arb.u.)",
            "real_title_real": "<channel>",
            "real_title_imag": "<channel>",
            "fourier_cmap": "viridis",
            "fourier_cbar_label": "Intensity (arb.u.)",
            "fourier_title": "Fourier transform",
            "gauss_blurred_title": "Blurred <channel>"
        }
        # Todo: add more parameters to the dictionary
        # make a similar file for the snom plotter app and overwrite the defaults from the snom anlaysis package
        # make it possible to add mutliple sets of parameters, each for a different filetype
        '''
        channel indicators
        channel labels
        channel prefixes
        channel suffixes
        file endings (.gsf, .txt, .ascii, ...)
        synccorrected channel indicator
        manipulated channel indicator
        filetype indicator? (standard, aachen, comsol, ...)
        parameters type indicator? (txt, html, gsf)
        add all plotting parameters
        enable/disable logfiles
        standard channels
        also add the default values for the loading of the data like:
            phaseoffset
            rounding_decimal (amp, phase, height, ...)
            scaling
        allow to add a list of custom channels which will be added to all_channels_custom
        '''
        with open(self.plotting_parameters_path, 'w') as file:
            json.dump(dictionary, file, indent=4)
   
    def _get_channel_tag_dict_value(self, channel:str, tag:ChannelTags) -> list:
        """This function returns the value of the specified tag for the specified channel. If the tag is not found, it will return None.

        Args:
            channel (str): channel name
            tag (ChannelTags): tag name

        Returns:
            list: tag value or values as a list
        """
        # check if the tag is in the channel tag dict
        if tag not in self.channel_tag_dict[self.channels.index(channel)]:
            return [None]
        value = self.channel_tag_dict[self.channels.index(channel)][tag]
        # check if a unit is part of the value
        if isinstance(value, list):
            for element in value:
                # if a unit is part of the value it must be in first place
                if isinstance(element, str):
                    return value[1:] 
                else:
                    return value
        else:
            if isinstance(value, str):
                return [None]
            else:
                return [value]

    def _get_measurement_tag_dict_value(self, tag:MeasurementTags) -> list:
        """This function returns the value of the specified tag for the current measurement. If the tag is not found, it will return None.

        Args:
            channel (str): channel name
            tag (MeasurementTags): tag name

        Returns:
            list: tag value or values as a list
        """
        # check if the tag is in the measurement tag dict
        if tag not in self.measurement_tag_dict:
            return [None]
        value = self.measurement_tag_dict[tag]
        # check if a unit is part of the value
        if isinstance(value, list):
            for element in value:
                # if a unit is part of the value it must be in first place
                if isinstance(element, str):
                    return value[1:] 
                else:
                    return value
        else:
            if isinstance(value, str):
                return [None]
            else:
                return [value]

    def _get_channel_tag_dict_unit(self, channel:str, tag:ChannelTags) -> str:
        """This function returns the value of the specified tag for the specified channel. If the tag is not found, it will return None.

        Args:
            channel (str): channel name
            tag (ChannelTags): tag name

        Returns:
            float: tag unit if there is one
        """
        # check if the tag is in the channel tag dict
        if tag not in self.channel_tag_dict[self.channels.index(channel)]: 
            return None
        value = self.channel_tag_dict[self.channels.index(channel)][tag]
        # check if a unit is part of the value
        if isinstance(value, list):
            for element in value:
                # if a unit is part of the value it must be in first place
                if isinstance(element, str):
                    return value[0] 
                else:
                    return None
        else:
            if isinstance(value, str):
                return value
            else:
                return None

    def _get_measurement_tag_dict_unit(self, tag:MeasurementTags) -> str:
        """This function returns the value of the specified tag for the current measurement. If the tag is not found, it will return None.

        Args:
            channel (str): channel name
            tag (MeasurementTags): tag name

        Returns:
            float: tag unit if there is one
        """
        # check if the tag is in the measurement tag dict
        if tag not in self.measurement_tag_dict:
            return None
        value = self.measurement_tag_dict[tag]
        # check if a unit is part of the value
        if isinstance(value, list):
            for element in value:
                # if a unit is part of the value it must be in first place
                if isinstance(element, str):
                    return value[0] 
                else:
                    return None
        else:
            if isinstance(value, str):
                return value
            else:
                return None

    def _set_channel_tag_dict_value(self, channel:str, tag:ChannelTags, value) -> None:
        """This function sets the value of the specified tag for the specified channel.
        It automatically tries to keep the unit of the value if there is one. 

        Args:
            channel (str): channel name
            tag (ChannelTags): tag name
            value (list): tag values as a list, or single value
        """
        # ckeck if value is a list
        if isinstance(value, list):
            # check that no strings are in the list
            if isinstance(value[0], str):
                print('One of the provided values is a string, use set_channel_tag_dict_unit to change the unit!')
                return 0
            else:
                # try to get old unit
                unit = self._get_channel_tag_dict_unit(channel, tag)
                if unit is not None:
                    new_value = [unit] + value
                else:
                    new_value = value
                # set the new values
                self.channel_tag_dict[self.channels.index(channel)][tag] = new_value
        else:
            # check if unit is provided
            if isinstance(value, str):
                # dont add the str value to the channel dict, if a unit should be changed use the set_channel_tag_dict_unit function
                print('Provided value is a string, use set_channel_tag_dict_unit to change the unit!')
                return 0
            else:
                # set new value
                self.channel_tag_dict[self.channels.index(channel)][tag] = value

    def _set_measurement_tag_dict_value(self, tag:MeasurementTags, value) -> None:
        """This function sets the value of the specified tag for the current measurement.
        It automatically tries to keep the unit of the value if there is one. 

        Args:
            tag (MeasurementTags): tag name
            value (list): tag values as a list, or single value
        """
        # ckeck if value is a list
        if isinstance(value, list):
            # check that no strings are in the list
            if isinstance(value[0], str):
                print('One of the provided values is a string, use set_channel_tag_dict_unit to change the unit!')
                return 0
            else:
                # try to get old unit
                unit = self._get_measurement_tag_dict_unit(tag)
                if unit is not None:
                    new_value = [unit] + value
                else:
                    new_value = value
                # set the new values
                self.measurement_tag_dict[tag] = new_value
        else:
            # check if unit is provided
            if isinstance(value, str):
                # dont add the str value to the channel dict, if a unit should be changed use the set_channel_tag_dict_unit function
                print('Provided value is a string, use set_channel_tag_dict_unit to change the unit!')
                return 0
            else:
                # set new value
                self.measurement_tag_dict[tag] = value

    def _set_channel_tag_dict_unit(self, channel:str, tag:ChannelTags, value:str) -> None:
        """This function sets the unit of the specified tag for the specified channel.

        Args:
            channel (str): channel name
            tag (ChannelTags): tag name
            value (str): unit of the specified tag
        """
        # check if old unit exists
        old_unit = self._get_channel_tag_dict_unit(channel, tag)
        if old_unit is None:
            print('This filtype has no unit for the specified tag!\nSetting the unit anayways...')
            old_values = self._get_channel_tag_dict_value(channel, tag) # shift the old values to the right
            # check if old values are in a list
            if isinstance(old_values, list):
                new_values = [value] + old_values
            else:
                new_values = [value, old_values]
            self.channel_tag_dict[self.channels.index(channel)][tag] = new_values
        else:
            self.channel_tag_dict[self.channels.index(channel)][tag][0] = value

    def _set_measurement_tag_dict_unit(self, tag:MeasurementTags, value:str) -> None:
        """This function sets the unit of the specified tag for the current measurement.

        Args:
            tag (MeasurementTags): tag name
            value (str): unit of the specified tag
        """
        # check if old unit exists
        old_unit = self._get_measurement_tag_dict_unit(tag)
        if old_unit is None:
            print('This filtype has no unit for the specified tag!\nSetting the unit anayways...')
            old_values = self._get_measurement_tag_dict_value(tag) # shift the old values to the right
            # check if old values are in a list
            if isinstance(old_values, list):
                new_values = [value] + old_values
            else:
                new_values = [value, old_values]
            self.measurement_tag_dict[tag] = new_values
        else:
            self.measurement_tag_dict[tag][0] = value

    def _get_tagval(self, content, tag):
        """This function gets the value of the tag listed in the file header"""
        content_array = content.split('\n')
        tag_array = []
        tagval = 0 # if no tag val can be found return 0
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
            is_unit = False
            try: float(tag_array[i][1])
            except: is_unit = True
            if tag_array[i][0] == tag:
                if is_unit:
                    tagval = tag_array[i][1]
                else:
                    tagval = float(tag_array[i][1])
        return tagval

    def _get_optomechanical_indicator(self, channel) -> tuple:
        """This function returns the optomechanical indicator of the channel and its index in the channel name.
        Meaning it tries to find out wether the cannel is an optical or mechanical channel."""
        channel_list = list(channel)
        indicator = None
        indicator_index = None
        if self.optical_indicator != None and self.mechanical_indicator != None:
            for i in range(len(channel_list)):
                opto_len = len(self.optical_indicator)
                mech_len = len(self.mechanical_indicator)
                if channel[i:i+opto_len][0] == self.optical_indicator:
                    indicator = self.optical_indicator
                    indicator_index = i # return the first occurence of the indicator
                    break
                elif channel[i:i+mech_len][0] == self.mechanical_indicator:	
                    indicator = self.mechanical_indicator
                    indicator_index = i
                    break
        else:
            indicator = None
            indicator_index = None
            # print('optomechanical indicator for this filetype is not yet implemented')
        # check that the channel is not a height channel
        if self.height_indicator in channel:
            indicator = None
            indicator_index = None
        return indicator, indicator_index

    def _is_amp_channel(self, channel) -> bool:
        """This function returns True if the channel is an amplitude channel, False otherwise."""
        optomechanical_indicator, indicator_index = self._get_optomechanical_indicator(channel)
        if optomechanical_indicator == self.optical_indicator and self.amp_indicator in channel:
            return True
        else:
            return False
            
    def _is_phase_channel(self, channel) -> bool:
        """This function returns True if the channel is a phase channel, False otherwise."""
        optomechanical_indicator, indicator_index = self._get_optomechanical_indicator(channel)
        if optomechanical_indicator == self.optical_indicator and self.phase_indicator in channel:
            return True
        else:
            return False

    def _is_complex_channel(self, channel) -> bool:
        """This function returns True if the channel is a complex channel, False otherwise."""
        optomechanical_indicator, indicator_index = self._get_optomechanical_indicator(channel)
        if optomechanical_indicator == self.optical_indicator and (self.real_indicator in channel or self.imag_indicator in channel):
            return True
        else:
            return False

    def _is_height_channel(self, channel) -> bool:
        """This function returns True if the channel is a height channel, False otherwise."""
        optomechanical_indicator, indicator_index = self._get_optomechanical_indicator(channel)
        if optomechanical_indicator == None and self.height_indicator in channel:
            return True
        else:
            return False

    def _channel_has_demod_num(self, channel) -> bool:
        """This function returns True if the channel has a demodulation number, False otherwise.

        Args:
            channel (str): channel name

        Returns:
            bool: _description_
        """
        # only amplitude, phase, complex and mechanical (amp, phase) channels can have a demodulation number not the height channels
        if self._is_amp_channel(channel) or self._is_phase_channel(channel) or self._is_complex_channel(channel):
            return True
        elif self._is_height_channel(channel):
            return False
        else:
            try:
                if channel in self.mechanical_channels:
                    return True
            except:
                print('unknown channel encountered in _channel_has_demod_num')
                return False

    def _get_demodulation_num(self, channel) -> int:
        """This function returns the demodulation number of the channel.
        So far for all known filetypes the demodulation number is the number behind the optomechanical indicator (O or M) in the channel name."""
        optomechanical_indicator, indicator_index = self._get_optomechanical_indicator(channel)
        demodulation_num = None
        if indicator_index != None: # if the index is None the channel is a height channel and has no demodulation number
            demodulation_num = int(channel[indicator_index +1 : indicator_index +2])

        if demodulation_num == None and self._channel_has_demod_num(channel):
            # height channel for example has no demodulation number but should not cause an error
            print('demodulation number could not be found')
        return demodulation_num
    
    def _initialize_measurement_channel_indicators(self):
        """This function initializes the channel indicators for the measurement channels.
        More precisely it loades all the parameters from the config file.
        """
        # the cannel prefix and suffix are characters surrounding the channel name in the filename, they will be used when loading and saving the data
        # filename = directory_name + channel_prefix + channel + channel_suffix + appendix + '.gsf' (or '.txt') 
        # appendix is just a standard appendix when saving to not overwrite the original files, can be changed by the user default is '_manipulated'
        # new approach based on cofigfile
        self.phase_channels = self._get_from_config('phase_channels')
        self.amp_channels = self._get_from_config('amp_channels')
        self.real_channels = self._get_from_config('real_channels')
        self.imag_channels = self._get_from_config('imag_channels')
        self.complex_channels = self.imag_channels + self.real_channels
        self.height_channel = self._get_from_config('height_channel')
        self.height_channels = self._get_from_config('height_channels')
        self.mechanical_channels = self._get_from_config('mechanical_channels')
        self.all_channels_default = self.phase_channels + self.amp_channels + self.mechanical_channels
        self.preview_ampchannel = self._get_from_config('preview_ampchannel')
        self.preview_phasechannel = self._get_from_config('preview_phasechannel')
        self.preview_channels = self._get_from_config('preview_channels')
        self.height_indicator = self._get_from_config('height_indicator')
        self.amp_indicator = self._get_from_config('amp_indicator')
        self.phase_indicator = self._get_from_config('phase_indicator')
        self.backwards_indicator = self._get_from_config('backwards_indicator')
        self.real_indicator = self._get_from_config('real_indicator')
        self.imag_indicator = self._get_from_config('imag_indicator')
        self.optical_indicator = self._get_from_config('optical_indicator')
        self.mechanical_indicator = self._get_from_config('mechanical_indicator')
        self.channel_prefix_default = self._get_from_config('channel_prefix_default')
        self.channel_prefix_custom = self._get_from_config('channel_prefix_custom')
        self.channel_suffix_default = self._get_from_config('channel_suffix_default')
        self.channel_suffix_custom = self._get_from_config('channel_suffix_custom')
        self.channel_suffix_synccorrected_phase = self._get_from_config('channel_suffix_synccorrected_phase')
        self.channel_suffix_manipulated = self._get_from_config('channel_suffix_manipulated')
        self.channel_suffix_overlain = self._get_from_config('channel_suffix_overlain')
        self.file_ending = self._get_from_config('file_ending')
        self.phase_offset_default = self._get_from_config('phase_offset_default')
        self.phase_offset_custom = self._get_from_config('phase_offset_custom')
        self.rounding_decimal_amp_default = self._get_from_config('rounding_decimal_amp_default')
        self.rounding_decimal_amp_custom = self._get_from_config('rounding_decimal_amp_custom')
        self.rounding_decimal_phase_default = self._get_from_config('rounding_decimal_phase_default')
        self.rounding_decimal_phase_custom = self._get_from_config('rounding_decimal_phase_custom')
        self.rounding_decimal_complex_default = self._get_from_config('rounding_decimal_complex_default')
        self.rounding_decimal_complex_custom = self._get_from_config('rounding_decimal_complex_custom')
        self.rounding_decimal_height_default = self._get_from_config('rounding_decimal_height_default')
        self.rounding_decimal_height_custom = self._get_from_config('rounding_decimal_height_custom')
        self.height_scaling_default = self._get_from_config('height_scaling_default')
        self.height_scaling_custom = self._get_from_config('height_scaling_custom')
        
        # additional definitions independent of filetype:
        self.filter_gauss_indicator = 'gauss'
        self.filter_fourier_indicator = 'fft'

        #create also lists for the overlain channels
        self.overlain_phase_channels = [channel+'_overlain' for channel in self.phase_channels]
        self.corrected_phase_channels = [channel+'_corrected' for channel in self.phase_channels]
        self.corrected_overlain_phase_channels = [channel+'_corrected_overlain' for channel in self.phase_channels]
        self.overlain_amp_channels = [channel+'_overlain' for channel in self.amp_channels]
        
        self.all_channels_custom = self.height_channels + self.complex_channels + self.overlain_phase_channels + self.overlain_amp_channels + self.corrected_phase_channels + self.corrected_overlain_phase_channels
  
    def _create_channel_tag_dict(self, channels:list=None):
        """This function reads in the header of the gsf file for the specified channel and extracts the tag values. The tag values are stored in a dictionary for each channel.
        This tag dict is very similar to the measurement_tag_dict, but the measurement_tag_dict is only created on the basis of the parameter file.
        If individual channels have been modified this will only be stored in the channel_tag_dict.

        Args:
            channel (str): channel name for which the tag values should be extracted
        """
        if channels == None:
            channels = self.channels
        # create a list containing the tag dictionary for each channel
        self.channel_tag_dict = []
        for channel in channels:
            channel_dict = {}
            if channel in self.all_channels_default:
                suffix = self.channel_suffix_default
                prefix = self.channel_prefix_default
                channel_type = 'default'
            elif channel in self.all_channels_custom:
                suffix = self.channel_suffix_custom
                prefix = self.channel_prefix_custom
                channel_type = 'custom'
            else:
                print(f'channel {channel} not found in default or custom channels!')
                # assume it is a custom channel and try loading anyways
                suffix = self.channel_suffix_custom
                prefix = self.channel_prefix_custom
                channel_type = 'custom'
                # exit()
            # we want to read the non binary part of the datafile
            if self.file_ending == '.gsf':
                encod = 'latin1'
            elif self.file_ending == '.ascii':
                encod = 'latin1'
            else:
                pass
                # not necessarily a problem, since the creation of the channel tag dict is also a test if the correct filetype was found
                # print('file ending not supported')
                # print('in _create_channel_tag_dict')
            with open(self.directory_name / Path(self.filename.name + f'{prefix}{channel}{suffix}{self.file_ending}'), 'r', encoding=encod) as f:
                content=f.read()


            try: channel_tags = self._get_from_config('channel_tags')
            except:
                # seem like there are no channel tags in the config file
                # so we will just use the measurement tags to initialize the channel tags
                measurement_tags = self._get_from_config('measurement_tags', self.file_type)
                measurement_tag_enums = list(MeasurementTags)
                channel_tag_enums = list(ChannelTags)
                for key, value in measurement_tags.items():
                    # get the data from the measurement tag dict
                    for i in range(len(measurement_tag_enums)):
                        if key == measurement_tag_enums[i].name:
                            data = self.measurement_tag_dict[measurement_tag_enums[i]]
                    # insert the data into the channel tag dict with the corresponding key wich is an enum of the channel tags class but has the same name as the measurement tag
                    for i in range(len(channel_tag_enums)):
                        if key == channel_tag_enums[i].name:
                            channel_dict[channel_tag_enums[i]] = data

            else:
                for key, tag in channel_tags.items():
                    is_list = False
                    tag_value_found = False
                    value = None
                    values = [None]
                    if isinstance(tag, list):
                        is_list = True
                    # so far each tag contains a maximum of 2 values
                    if is_list:
                        values = []
                        for element in tag:
                            try: value = self._get_tagval(content, element)
                            except: 
                                values.append(None)
                                tag_value_found = False
                            else: 
                                values.append(value)
                                tag_value_found = True
                    else:
                        try: value = self._get_tagval(content, tag)
                        except: value = None
                        else: tag_value_found = True
                        # try to find out if the value is a number or a unit
                        try: float(value)
                        except: pass
                    # check if tag value was found
                    if not tag_value_found:
                        print(f'Could not find the tag value for {tag} in channel {channel}. You should probably check the config file.')
                        continue
                    if key == 'PIXELAREA':
                        try: channel_dict[ChannelTags.PIXELAREA] = [int(values[0]), int(values[1]), int(values[2])]
                        except: channel_dict[ChannelTags.PIXELAREA] = [int(values[0]), int(values[1])]
                    elif key == 'YINCOMPLETE':
                        channel_dict[ChannelTags.YINCOMPLETE] = int(value)
                    elif key == 'SCANNERCENTERPOSITION':
                        try: channel_dict[ChannelTags.SCANNERCENTERPOSITION] = [float(values[0]), float(values[1]), float(values[2])]
                        except: channel_dict[ChannelTags.SCANNERCENTERPOSITION] = [float(values[0]), float(values[1])]
                    elif key == 'ROTATION':
                        channel_dict[ChannelTags.ROTATION] = float(value)
                    elif key == 'SCANAREA':
                        try: channel_dict[ChannelTags.SCANAREA] = [float(values[0]), float(values[1]), float(values[2])]
                        except: channel_dict[ChannelTags.SCANAREA] = [float(values[0]), float(values[1])]
                    elif key == 'XYUNIT':
                        channel_dict[ChannelTags.XYUNIT] = value
                    elif key == 'ZUNIT':
                        channel_dict[ChannelTags.ZUNIT] = value
                    elif key == 'WAVENUMBERSCALING':
                        channel_dict[ChannelTags.WAVENUMBERSCALING] = float(value)
            # add pixel scaling to the channel dict, initially this is always 1
            channel_dict[ChannelTags.PIXELSCALING] = 1
                    
            self.channel_tag_dict.append(channel_dict)


# this could be split in AFM and SNOM measurement classes where AFM has all the base functions and SNOM inherits from it
# make it easier for AFM users to finde the functions they need
class SnomMeasurement(FileHandler):
    """This class opens a snom measurement and handels all the snom related functions.
    
    Args:
        directory_name (str): path to the directory containing the measurement
        channels (list, optional): list of channels to load. Defaults to None.
        title (str, optional): title of the measurement. Defaults to None.
        autoscale (bool, optional): if True the data will be scaled to quadratic pixels. Defaults to True.
    """
    def __init__(self, directory_name:str, channels:list=None, title:str=None, autoscale:bool=True) -> None:
        self.all_subplots = [] # list containing all subplots
        self.measurement_type = MeasurementTypes.SNOM
        super().__init__(directory_name, title)
        self._initialize_measurement_channel_indicators()
        if channels == None: # the standard channels which will be used if no channels are specified
            channels = self.preview_channels
        self.channels = channels.copy() # make sure to copy the list to avoid changing the original list     
        self.autoscale = autoscale
        self._initialize_data(self.channels)
        if PlotDefinitions.autodelete_all_subplots: self._delete_all_subplots() # automatically delete old subplots
        # get the plotting style from the mpl style file
        self._load_mpl_style()
    
    def _initialize_data(self, channels=None) -> None:
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
            self._create_channel_tag_dict()
            self.all_data, self.channels_label = self._load_data(channels) # could be changed to a single dictionary containing the data and the channel names
            xres = len(self.all_data[0][0])
            yres = len(self.all_data[0])
            # reset all the instance variables dependent on the data, but not the ones responsible for plotting
            if self.autoscale == True:
                self.quadratic_pixels()
            # initialize instance variables:
            self.mask_array = [] # not shure if it's best to reset the mask...
            self.align_points = None
            self.scalebar_channels = []    

    def initialize_channels(self, channels:list) -> None:
        """This function will load the data from the specified channels and replace the ones in memory.
        
        Args:
            channels [list]: a list containing the channels you want to initialize
        """
        self._initialize_data(channels)

    def add_channels(self, channels:list) -> None:
        """This function will add the specified channels to memory without changing the already existing ones.

        Args:
            channels (list): Channels to add to memory.
        """
        self.channels += channels
        # update the channel tag dictionary, makes the program compatible with differrently sized datasets, like original data plus manipulated, eg. cut data
        self._create_channel_tag_dict(channels)
        all_data, channels_label = self._load_data(channels)
        for i in range(len(channels)):
            self.all_data.append(all_data[i])
            self.channels_label.append(channels_label[i])
        # reset all the instance variables dependent on the data, but nor the ones responsible for plotting
        # self.scaling_factor = 1
        if self.autoscale == True:
            self.quadratic_pixels(channels)

    def _load_all_subplots(self) -> None:
        """Load all subplots from memory (located under home/SNOM_Analysis/all_subplots.p).
        """
        try:
            with open(self.all_subplots_path, 'rb') as file:
                self.all_subplots = pkl.load(file)
        except: self.all_subplots = []
         
    def _export_all_subplots(self) -> None:
        """Export all subplots to memory.
        """
        with open(self.all_subplots_path, 'wb') as file:
            pkl.dump(self.all_subplots, file)
        self.all_subplots = []

    def _delete_all_subplots(self):
        """Delete the subplot memory. Should be done always if new measurement row is investigated.
        """
        try:
            os.remove(self.all_subplots_path)
        except: pass
        self.all_subplots = []
        
    def _scale_array(self, array, scaling) -> np.array:
        """This function scales a given 2D Array, it thus creates 'scaling'**2 subpixels per pixel.
        The scaled array is returned."""
        yres = len(array)
        xres = len(array[0])
        scaled_array = np.zeros((yres*scaling, xres*scaling))
        for i in range(len(array)):
            for j in range(len(array[0])):
                for k in range(scaling):
                    for l in range(scaling):
                        scaled_array[i*scaling + k][j*scaling + l] = array[i][j]
        return scaled_array

    def scale_channels(self, channels:list=None, scaling:int=4) -> None:
        """This function scales all the data in memory or the specified channels.
                
        Args:
            channels (list, optional): List of channels to scale. If not specified all channels in memory will be scaled. Defaults to None.
            scaling (int, optional): Defines scaling factor. Each pixel will be scaled to scaling**2 subpixels. Defaults to 4.
        """
        if channels is None:
            channels = self.channels
        self._write_to_logfile('scaling', scaling)
        for channel in channels:
            if channel in self.channels:
                self.all_data[self.channels.index(channel)] = self._scale_array(self.all_data[self.channels.index(channel)], scaling)
                XReal, YReal, *args = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELAREA)
                self._set_channel_tag_dict_value(channel, ChannelTags.PIXELAREA, [XReal*scaling, YReal*scaling])
                self._set_channel_tag_dict_value(channel, ChannelTags.PIXELSCALING, scaling)
            else:
                print(f'Channel {channel} is not in memory! Please initiate the channels you want to use first!')

    def _load_data(self, channels:list) -> list:
        """Loads all binary data of the specified channels and returns them in a list plus the dictionary with the channel information.
        Height data is automatically converted to nm. 
        
        Args:
            channels (list): list of channels to load
        """

        data_dict = []
        all_data = []
        for channel in channels:
            # check if channel is a default channel or something user made
            # if default use the standard naming convention
            # if user made dont use the '_raw' suffix
            if channel in self.all_channels_default:
                suffix = self.channel_suffix_default
                prefix = self.channel_prefix_default
                channel_type = 'default'
            elif channel in self.all_channels_custom:
                suffix = self.channel_suffix_custom
                prefix = self.channel_prefix_custom
                channel_type = 'custom'
            else:
                print(f'channel {channel} not found in default or custom channels!')
                # assume it is a custom channel and try loading anyways
                suffix = self.channel_suffix_custom
                prefix = self.channel_prefix_custom
                channel_type = 'custom'
                # exit()
            # check the readmode depending on the filetype
            # this also affects the way the data is read and processed
            if self.file_ending == '.gsf':
                read_mode = 'br'
            elif self.file_ending == '.ascii':
                read_mode = 'r'
            else:
                print('file ending not supported')
            with open(self.directory_name / Path(self.filename.name + f'{prefix}{channel}{suffix}{self.file_ending}'), read_mode) as f:
                data=f.read()

            if read_mode == 'br':
                binarydata = data
            elif read_mode == 'r':
                datalist = data.split('\n')
                datalist = [element.split(' ') for element in datalist]
                datalist = np.array(datalist[:-1], dtype=float)#, dtype=np.float convert list to np.array and strings to float
            
            # get the resolution of the channel 
            XRes, YRes, *args = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELAREA)
            datasize=int(XRes*YRes*4)
            channel_data = np.zeros((YRes, XRes))
            # we knwo the resolution of the data from the header or parameter file
            # we use that to read the data from the end of the file until the end of the file minus the datasize
            # in this way we ignore the header and read only the data
            if read_mode == 'br':
                reduced_binarydata=binarydata[-datasize:]

            # depending on the channel type set the scaling, phase_offset and rounding_decimal
            scaling = 1 # default scaling, not every channel needs scaling
            phase_offset = 0 # default phase offset, not every channel needs a phase offset
            if self._is_amp_channel(channel):
                if channel_type == 'default':
                    rounding_decimal = self.rounding_decimal_amp_default
                elif channel_type == 'custom':
                    rounding_decimal = self.rounding_decimal_amp_custom
            if self._is_height_channel(channel):
                if channel_type == 'default':
                    scaling = self.height_scaling_default
                    rounding_decimal = self.rounding_decimal_height_default
                elif channel_type == 'custom':
                    scaling = self.height_scaling_custom
                    rounding_decimal = self.rounding_decimal_height_custom
            if self._is_phase_channel(channel):
                if channel_type == 'default':
                    phase_offset = self.phase_offset_default
                    rounding_decimal = self.rounding_decimal_phase_default
                elif channel_type == 'custom':
                    phase_offset = self.phase_offset_custom
                    rounding_decimal = self.rounding_decimal_phase_custom
            if self._is_complex_channel(channel):
                if channel_type == 'default':
                    rounding_decimal = self.rounding_decimal_complex_default
                elif channel_type == 'custom':
                    rounding_decimal = self.rounding_decimal_complex_custom
            # print(f'channel: {channel} is a {channel_type} channel')
            # print(f'channel: {channel} is a amp channel ', self._is_amp_channel(channel))
            # print(f'channel: {channel} is a phase channel ', self._is_phase_channel(channel))
            # print(f'channel: {channel} is a height channel ', self._is_height_channel(channel))
            # print(f'channel: {channel}, scaling: {scaling}, phase_offset: {phase_offset}, rounding_decimal: {rounding_decimal}')

            # now read the data and apply the scaling, phase offset and rounding
            for y in range(0,YRes):
                for x in range(0,XRes):
                    if read_mode == 'br':
                        pixval = unpack("f",reduced_binarydata[4*(y*XRes+x):4*(y*XRes+x+1)])[0]
                        channel_data[y][x] = round(pixval*scaling + phase_offset, rounding_decimal)
                    elif read_mode == 'r':
                        channel_data[y][x] = round(datalist[y][x]*scaling + phase_offset, rounding_decimal)
            
            all_data.append(channel_data)
            data_dict.append(channel)
        # data_dict currently is just a list of the channels, this list is not equivalent to self.channels as the data_dict
        # or later self.channels_label contains the names of the channels which are used as the plot title, they will change depending on the functions applied, eg. 'channel_blurred' or channel_manipulated'...
        # but self.channels will always contain the original channel name as this is used for internal referencing
        return all_data, data_dict

    def _load_data_binary(self, channels:list) -> list:
        """Loads all binary data of the specified channels and returns them in a list plus the dictionary for access.
        
        Args:
            channels (list): list of channels to load
        """
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

    def set_min_to_zero(self, channels:list=None) -> None:
        """This function sets the min value of the specified channels to zero.
                
        Args:
            channels (list, optional): List of channels to set min value to zero. If not specified this will apply to all height channels in memory. Defaults to None.
        """
        if channels is None:
            channels = []
            for channel in self.channels:
                if self.height_indicator in channel:
                    channels.append(channel)

        self._write_to_logfile('set_min_to_zero', True)
        for channel in channels:
            if channel in self.channels:
                data = self.all_data[self.channels.index(channel)]
                flattened_data = data.flatten()
                data_min = min(flattened_data)
                self.all_data[self.channels.index(channel)] = data - data_min
            else:
                print('At least one of the specified channels is not in memory! You probably should initialize the channels first.')

    def _get_plotting_values(self, channel:str) -> tuple:
        """This function returns the colormap, the colormap label and the title for the specified channel.

        Args:
            channel (str): channel name
        """
        # import plotting_parameters.json, here the user can tweek some options for the plotting, like automatic titles and colormap choices
        plotting_parameters = self._get_plotting_parameters()

        # update the placeholders in the dictionary
        # the dictionary contains certain placeholders, which are now being replaced with the actual values
        # until now only the channel placeholder is used but more could be added
        # placeholders are indicated by the '<' and '>' characters
        # this step insures, that for example the title contains the correct channel name
        placeholders = {'<channel>': channel}
        plotting_parameters = self._replace_plotting_parameter_placeholders(plotting_parameters, placeholders)
    
        if self.amp_indicator in channel and self.height_indicator not in channel:
            cmap = plotting_parameters["amplitude_cmap"]
            label = plotting_parameters["amplitude_cbar_label"]
            title = plotting_parameters["amplitude_title"]
        elif self.phase_indicator in channel:
            cmap = plotting_parameters["phase_cmap"]
            if 'positive' in channel:
                title = plotting_parameters["phase_positive_title"]
            elif 'negative' in channel:
                title = plotting_parameters["phase_negative_title"]
            else:
                title = plotting_parameters["phase_title"]
            label = plotting_parameters["phase_cbar_label"]
        elif self.height_indicator in channel:
            cmap = plotting_parameters["height_cmap"]
            label = plotting_parameters["height_cbar_label"]
            title = plotting_parameters["height_title"]
        elif self.real_indicator in channel or self.imag_indicator in channel:
            cmap = plotting_parameters["real_cmap"]
            label = plotting_parameters["real_cbar_label"]
            if self.real_indicator in channel:
                title = plotting_parameters["real_title_real"]
            else:
                title = plotting_parameters["real_title_imag"]
        elif self.filter_fourier_indicator in channel:
            cmap = plotting_parameters["fourier_cmap"]
            label = plotting_parameters["fourier_cbar_label"]
            title =  plotting_parameters["fourier_title"]
        elif self.filter_gauss_indicator in channel:
            title = plotting_parameters["gauss_blurred_title"]
        
        else:
            print('channel: ', channel)
            print('self.amp_indicator: ', self.amp_indicator)
            print('self.phase_indicator: ', self.phase_indicator)
            print('self.height_indicator: ', self.height_indicator)
            print('self.real_indicator: ', self.real_indicator)
            print('self.imag_indicator: ', self.imag_indicator)
            print('In _add_subplot(), encountered unknown channel')
            exit()
        return cmap, label, title

    def _add_subplot(self, data:np.array, channel:str, scalebar:list=None) -> list:
        """This function adds the specified data to the list of subplots. The list of subplots contains the data, the colormap,
        the colormap label and a title, which are generated from the channel information. The same array is also returned,
        so it can also be iterated by an other function to only plot the data of interest.
        
        Args:
            data (np.array): data to plot
            channel (str): channel name
            scalebar (list, optional): list of scalebar parameters. Defaults to None.
        
        Returns:
            list: [data, cmap, label, title, scalebar]
        """
        cmap, label, title = self._get_plotting_values(channel)
        # subplots.append([data, cmap, label, title])
        if self.measurement_title != None:
            title = self.measurement_title + title
        '''
        if scalebar != None:
            self.all_subplots.append([np.copy(data), cmap, label, title, scalebar])
            return [data, cmap, label, title, scalebar]
        else:
            self.all_subplots.append([np.copy(data), cmap, label, title])
            return [data, cmap, label, title]
        '''
        supplot = {'data': np.copy(data), 'cmap': cmap, 'label': label, 'title': title, 'scalebar': scalebar}
        self._load_all_subplots()
        self.all_subplots.append(supplot)
        self._export_all_subplots()
        return supplot
    
    def remove_subplots(self, index_array:list) -> None:
        """This function removes the specified subplot from the memory.
        
        Args:
            index_array (list): The indices of the subplots to remove from the plot list
        """
        #sort the index array in descending order and delete the corresponding plots from the memory
        index_array.sort(reverse=True)
        self._load_all_subplots()
        for index in index_array:
            del self.all_subplots[index]
        self._export_all_subplots()

    def remove_last_subplots(self, times:int=1) -> None:
        """This function removes the last added subplots from the memory.
        Times specifies how often the last subplot should be removed.
        Times=1 means only the last, times=2 means the two last, ...
        
        Args:
            times (int): how many subplots should be removed from the end of the list?
        """
        self._load_all_subplots()
        for i in range(times):
            self.all_subplots.pop()
        self._export_all_subplots()

    def _plot_subplots(self, subplots:list) -> None:
        """This function plots the subplots. The plots are created in a grid, by default the grid is optimized for 3 by 3.
        The layout changes dependent on the number of subplots of subplots and also the dimensions.
        Wider subplots are prefferably created vertically, otherwise they are plotted horizontally. Probably subject to future changes...
        
        Args:
            subplots (list): list of subplots to plot
        """
        number_of_axis = 9
        number_of_subplots = len(subplots)
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
        data = subplots[0]['data']
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
        
        # get the plotting style from the mpl style file and apply it
        self._load_mpl_style()

        #start the plotting process
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
                    data = subplots[counter]['data']
                    cmap = subplots[counter]['cmap']
                    label = subplots[counter]['label']
                    title = subplots[counter]['title']
                    scalebar = subplots[counter]['scalebar']
                    if scalebar is not None:
                        dx, units, dimension, scalebar_label, length_fraction, height_fraction, width_fraction, location, loc, pad, border_pad, sep, frameon, color, box_color, box_alpha, scale_loc, label_loc, font_properties, label_formatter, scale_formatter, fixed_value, fixed_units, animated, rotation = scalebar
                        scalebar = ScaleBar(dx, units, dimension, scalebar_label, length_fraction, height_fraction, width_fraction,
                            location, loc, pad, border_pad, sep, frameon, color, box_color, box_alpha, scale_loc,
                            label_loc, font_properties, label_formatter, scale_formatter, fixed_value, fixed_units, animated, rotation) 
                        axis.add_artist(scalebar)

                    #center the colorscale for real data around 0
                    # get minima and maxima from data:
                    flattened_data = data.flatten()
                    min_data = np.min(flattened_data)
                    max_data = np.max(flattened_data)
                    if self.real_indicator in title or self.imag_indicator in title: # for real part or imaginary part data
                        if self.file_type == 'FILETYPE6':
                            data = set_nan_to_zero(data) #comsol data can contain nan values which are problematic for min and max
                        data_limit = get_largest_abs(min_data, max_data)
                        if PlotDefinitions.vlimit_real is None: PlotDefinitions.vlimit_real = data_limit
                        if PlotDefinitions.real_cbar_range is True:
                            if PlotDefinitions.vlimit_real < data_limit: PlotDefinitions.vlimit_real = data_limit
                            img = axis.pcolormesh(data, cmap=cmap, vmin=-PlotDefinitions.vlimit_real, vmax=PlotDefinitions.vlimit_real, rasterized=True)
                        else:
                            img = axis.pcolormesh(data, cmap=cmap, vmin=-data_limit, vmax=data_limit, rasterized=True)
                    else:
                        if cmap == SNOM_phase and PlotDefinitions.full_phase_range is True: # for phase data
                            vmin = 0
                            vmax = 2*np.pi
                            img = axis.pcolormesh(data, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
                        elif cmap == SNOM_phase and PlotDefinitions.full_phase_range is False:
                            if PlotDefinitions.vmin_phase is None: PlotDefinitions.vmin_phase = min_data
                            if PlotDefinitions.vmax_phase is None: PlotDefinitions.vmax_phase = max_data
                            if PlotDefinitions.shared_phase_range is True:
                                if PlotDefinitions.vmin_phase > min_data: PlotDefinitions.vmin_phase = min_data
                                if PlotDefinitions.vmax_phase < max_data: PlotDefinitions.vmax_phase = max_data
                            else:
                                PlotDefinitions.vmin_phase = min_data
                                PlotDefinitions.vmax_phase = max_data
                            img = axis.pcolormesh(data, cmap=cmap, vmin=PlotDefinitions.vmin_phase, vmax=PlotDefinitions.vmax_phase, rasterized=True)
                            
                        elif cmap == SNOM_amplitude and PlotDefinitions.amp_cbar_range is True:
                            if PlotDefinitions.vmin_amp is None: PlotDefinitions.vmin_amp = min_data
                            if PlotDefinitions.vmax_amp is None: PlotDefinitions.vmax_amp = max_data
                            if min_data < PlotDefinitions.vmin_amp: PlotDefinitions.vmin_amp = min_data # update the min and max values in PlotDefinitions if new values are outside of range
                            if max_data > PlotDefinitions.vmax_amp: PlotDefinitions.vmax_amp = max_data
                            vmin = PlotDefinitions.vmin_amp
                            vmax = PlotDefinitions.vmax_amp
                            img = axis.pcolormesh(data, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
                        elif cmap == SNOM_height and PlotDefinitions.height_cbar_range is True:
                            if PlotDefinitions.vmin_height is None: PlotDefinitions.vmin_height = min_data # initialize for the first time
                            if PlotDefinitions.vmax_height is None: PlotDefinitions.vmax_height = max_data
                            if min_data < PlotDefinitions.vmin_height: PlotDefinitions.vmin_height = min_data # update the min and max values in PlotDefinitions if new values are outside of range
                            if max_data > PlotDefinitions.vmax_height: PlotDefinitions.vmax_height = max_data
                            vmin = PlotDefinitions.vmin_height
                            vmax = PlotDefinitions.vmax_height
                            img = axis.pcolormesh(data, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
                        else:
                            # print('not plotting full range phase')
                            img = axis.pcolormesh(data, cmap=cmap, rasterized=True)
                    
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
                                mean = self._get_mean_value(self.mask_array, x, y, border_width)
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
                    
                    # invert y axis to fit to the scanning procedure which starts in the top left corner
                    axis.invert_yaxis()
                    divider = make_axes_locatable(axis)
                    # cax = divider.append_axes("right", size=f"{self.colorbar_width}%", pad=0.05) # size is the size of colorbar relative to original axis, 100% means same size, 10% means 10% of original
                    cax = divider.append_axes("right", size=f"{calculate_colorbar_size(fig, axis, self.colorbar_width)}%", pad=0.05) # size is the size of colorbar relative to original axis, 100% means same size, 10% means 10% of original
                    cbar = plt.colorbar(img, aspect=1, cax=cax)
                    cbar.ax.get_yaxis().labelpad = 15
                    cbar.ax.set_ylabel(label, rotation=270)
                    if self.hide_ticks == True:
                        # remove ticks on x and y axis, they only show pixelnumber anyways, better to add a scalebar
                        axis.set_xticks([])
                        axis.set_yticks([])
                    if self.show_titles == True:
                        axis.set_title(title)
                    axis.axis('scaled')
                    counter += 1

        #turn off all unneeded axes
        counter = 0
        for row in range(nrows):
            for col in range(int(np.sqrt(number_of_axis))):
                if  counter >= number_of_subplots > 3 and changed_orientation == False and number_of_subplots != 4: 
                    axis = ax[row, col]
                    axis.axis('off')
                counter += 1

        plt.subplots_adjust(hspace=PlotDefinitions.hspace)
        if self.tight_layout is True:
            plt.tight_layout()
        if PlotDefinitions.show_plot is True:
            plt.show()
        gc.collect()
    
    def switch_supplots(self, first_id:int=None, second_id:int=None) -> None:
        """
        This function changes the position of the subplots.
        The first and second id corresponds to the positions of the two subplots which should be switched.
        This function will be repea you are satisfied.

        Args:
            first_id (int): the first id of the two subplots which should be switched
            second_id (int): the second id of the two subplots which should be switched
        """
        if (first_id == None) or (second_id == None):
            first_id = int(input('Please enter the id of the first image: '))
            second_id = int(input('Please enter the id of the second image: '))
        self._load_all_subplots()
        first_subplot = self.all_subplots[first_id]
        self.all_subplots[first_id] = self.all_subplots[second_id]
        self.all_subplots[second_id] = first_subplot
        self._export_all_subplots()
        self.display_all_subplots()
        print('Are you happy with the new positioning?')
        user_input = self._user_input_bool()
        if user_input == False:
            print('Do you want to change the order again?')
            user_input = self._user_input_bool()
            if user_input == False:
                exit()
            else:
                self.switch_supplots()

    def _display_dataset(self, dataset, channels) -> None:
        """Add all data contained in dataset as subplots to one figure.
        The data has to be shaped beforehand!
        channels should contain the information which channel is stored at which position in the dataset.

        Args:
            dataset (list): list of data to plot
            channels (list): list of channel names
        """
        subplots = []
        for i in range(len(dataset)):
            scalebar = None
            for j in range(len(self.scalebar_channels)):
                if self.channels[i] == self.scalebar_channels[j][0]:
                    scalebar = self.scalebar_channels[j][1]
            subplots.append(self._add_subplot(dataset[i], channels[i], scalebar))
        self._plot_subplots(subplots)

    def display_all_subplots(self) -> None:
        """
        This function displays all the subplots which have been created until this point.
        """
        self._load_all_subplots()
        self._plot_subplots(self.all_subplots)
        self.all_subplots = []
        gc.collect()

    def display_channels(self, channels:list=None) -> None: 
        """This function displays the channels in memory or the specified ones.
                
        Args:
            channels (list, optional): List of channels to display. If not specified all channels from memory will be plotted. Defaults to None.

        """
        if channels == None:
            dataset = self.all_data
            plot_channels = self.channels
        else:
            dataset = []
            plot_channels = []
            for channel in channels:
                if channel in self.channels:
                    dataset.append(self.all_data[self.channels.index(channel)])
                    plot_channels.append(channel)
                else: 
                    print(f'Channel {channel} is not in memory! Please initiate the channels you want to display first!')
                    print(self.channels)

        self._display_dataset(dataset, plot_channels)
        gc.collect()

    def display_overlay(self, channel1:str, channel2:str, alpha=0.5) -> None:
        """This function displays an overlay of two channels. The first channel is displayed in full color, the second channel is displayed width a specified alpha.

        Args:
            channel1 (str): channel name of the first channel
            channel2 (str): channel name of the second channel
            alpha (float, optional): alpha value of the second channel. Defaults to 0.5.
        """
        # get the colormaps
        cmap1, label1, title1 = self._get_plotting_values(channel1)
        cmap2, label2, title2 = self._get_plotting_values(channel2)
        # get the data
        data1 = self.all_data[self.channels.index(channel1)]
        data2 = self.all_data[self.channels.index(channel2)]
        # create the figure
        fig, ax = plt.subplots()
        fig.set_figheight(self.figsizey)
        fig.set_figwidth(self.figsizex)
        # plot the data
        img1 = ax.imshow(data1, cmap=cmap1)
        img2 = ax.imshow(data2, cmap=cmap2, alpha=alpha)
        # add the colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=f"{self.colorbar_width}%", pad=0.05)
        cbar = plt.colorbar(img1, aspect=1, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(label1, rotation=270)
        # invert y axis to fit to the scanning procedure which starts in the top left corner
        ax.invert_yaxis()
        # add the title
        # ax.set_title(title)
        # remove ticks on x and y axis, they only show pixelnumber anyways, better to add a scalebar
        if self.hide_ticks == True:
            ax.set_xticks([])
            ax.set_yticks([])
        plt.show()
        gc.collect()

    def _gauss_blurr_data(self, array, sigma) -> np.array:
        """Applies a gaussian blurr to the specified array, with a specified sigma. The blurred data is returned as a np.array."""
        return gaussian_filter(array, sigma)

    def gauss_filter_channels(self, channels:list=None, sigma=2):
        """This function will gauss filter the specified channels. If no channels are specified, the ones in memory will be used.
        Only for amplitude and height data, phase data will be ignored. Works fine, but the gauss_filter_channels_complex() function is more versatile.

        Args:
            channels (list, optional): List of channels to blurr, if not specified all channels will be blurred. Should not be used for phase. Defaults to None.
            sigma (int, optional): The 'width' of the gauss blurr in pixels, you should scale the data before blurring. Defaults to 2.
        """
        if channels is None:
            channels = self.channels
        self._write_to_logfile('gaussian_filter_sigma', sigma)
        
        # start the blurring:
        for channel in channels:
            if channel in self.channels:
                channel_index = self.channels.index(channel)
                # check pixel scaling from channel tag dict for each channel
                pixel_scaling = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELSCALING)[0]
                if pixel_scaling == 1:
                    if PlotDefinitions.show_plot:
                        print(f'The data in channel {channel} is not yet scaled! Do you want to scale the data?')
                        user_input = self._user_input_bool()
                        if user_input == True:
                            self.scale_channels([channel])
                self.all_data[channel_index] = self._gauss_blurr_data(self.all_data[channel_index], sigma)
                self.channels_label[channel_index] += '_' + self.filter_gauss_indicator
            else: 
                print(f'Channel {channel} is not in memory! Please initiate the channels you want to use first!')

    def _find_gauss_compatible_channels(self) -> list:
        """This function goes through all channels in memory and tries to find compatible pairs of amplitude and phase channels.
        The function returns a list of lists, where each sublist contains the indices of the amplitude and phase channel.
        """
        channel_pairs = [] # list of lists, where each sublist contains the indices of the amplitude and phase channel relative to the self.channels list
        phase_channels = [] # sort the phase channels in a separate list
        amp_channels = [] # sort the amplitude channels in a separate list e.g. [[demod, channel_index, channel_name]]
        for i in range(len(self.channels)):
            demod = self._get_demodulation_num(self.channels[i])
            if self._is_amp_channel(self.channels[i]):
                amp_channels.append([demod, i])
            elif self._is_phase_channel(self.channels[i]):
                phase_channels.append([demod, i])

        # now try to find a partner for each phase channel, if there are amp channels without a partner they will be blurred ignoring the phase
        for i in range(len(phase_channels)):
            possible_amp_partners = []
            for j in range(len(amp_channels)):
                if phase_channels[i][0] == amp_channels[j][0]: # check if the demodulation number is the same
                    if self.all_data[phase_channels[i][1]].shape == self.all_data[amp_channels[j][1]].shape: # check if the data shape is the same
                        possible_amp_partners.append(amp_channels[j][1])
            if len(possible_amp_partners) == 1:
                channel_pairs.append([possible_amp_partners[0], phase_channels[i][1]])
            elif len(possible_amp_partners) > 1:
                print(f'Found more than one possible amplitude channel for phase channel {self.channels[phase_channels[i][1]]}!')
                print('Please specify the correct one! This channel will be ignored for now.')
        
        return channel_pairs

    def gauss_filter_channels_complex(self, channels:list=None, scaling:int=4, sigma:int=2) -> None:
        """This fucton gauss filters the specified channels. If no channels are specified, all channels in memory will be used.
        The function is designed to work with complex data, where amplitude and phase are stored in separate channels.
        It will also blurr height, real part and imaginary part channels and amplitude channels without phase partner and phase channels without amplitude partner if you want to.
        If the data is not scaled already the function will do it automatically, the default scaling factor is 4, works good with sigma=2.
                
        Args:
            channels [list]: list of channels to blurr, must contain amplitude and phase of same channels.
            scaling [int]: the scaling factor used for scaling the data, default is 4
            sigma [int]: the sigma used for blurring the data, bigger sigma means bigger blurr radius

        """
        self._write_to_logfile('gaussian_filter_complex_sigma', sigma)
        if channels is None:
            channels = self.channels
        for channel in channels:
            if channel not in self.channels:
                print(f'Channel {channel} is not in memory! Please initiate the channels you want to use first!')

        # get pairs of amplitude and phase channels
        channel_pairs = self._find_gauss_compatible_channels()
        # make a list of the remaining channels
        remaining_channels = []
        for i in range(len(self.channels)):
            if i not in [pair[0] for pair in channel_pairs] and i not in [pair[1] for pair in channel_pairs]:
                if self._is_phase_channel(self.channels[i]) == False: # ignore phase channels
                    remaining_channels.append(i)
                else:
                    print(f'Channel {self.channels[i]} is a phase channel and does not have a compatible amplitude channel!')
                    print('For phase data without amplitude please use the gauss_filter_channels() function!')
                    # get user input if the phase channel should be blurred without amplitude, might be useful in some cases when the phase is flat
                    print('Do you want to blur this channel without amplitude anyways?')
                    user_input = self._user_input_bool()
                    if user_input == True:
                        remaining_channels.append(i)
        
        # check if the data is scaled, if not scale it
        for i in range(len(channel_pairs)):
            if self._get_channel_tag_dict_value(self.channels[channel_pairs[i][0]], ChannelTags.PIXELSCALING)[0] == 1:
                # scale the data
                self.scale_channels([self.channels[channel_pairs[i][0]]], scaling)
            if self._get_channel_tag_dict_value(self.channels[channel_pairs[i][1]], ChannelTags.PIXELSCALING)[0] == 1:
                # scale the data
                self.scale_channels([self.channels[channel_pairs[i][1]]], scaling)
        
        for i in range(len(remaining_channels)):
            if self._get_channel_tag_dict_value(self.channels[remaining_channels[i]], ChannelTags.PIXELSCALING)[0] == 1:
                # scale the data
                self.scale_channels([self.channels[remaining_channels[i]]], scaling)

        # now start the blurring process for the amplitude and phase channel pairs
        print('Starting the blurring process, this might take a while...')
        for i in range(len(channel_pairs)):
            amp = self.all_data[channel_pairs[i][0]]
            phase = self.all_data[channel_pairs[i][1]]
            real = amp*np.cos(phase)
            imag = amp*np.sin(phase)

            real_blurred = self._gauss_blurr_data(real, sigma)
            imag_blurred = self._gauss_blurr_data(imag, sigma)
            compl_blurred = np.add(real_blurred, 1J*imag_blurred)
            amp_blurred = np.abs(compl_blurred)
            phase_blurred = self._get_compl_angle(compl_blurred)

            # update the data in memory and the labels used for plotting but not the channel names
            self.all_data[channel_pairs[i][0]] = amp_blurred
            self.channels_label[channel_pairs[i][0]] = self.channels_label[channel_pairs[i][0]] + '_' + self.filter_gauss_indicator
            self.all_data[channel_pairs[i][1]] = phase_blurred
            self.channels_label[channel_pairs[i][1]] = self.channels_label[channel_pairs[i][1]] + '_' + self.filter_gauss_indicator

        # now start the blurring process for the remaining channels
        # this will blurr height, real part, imaginary part channels and amplitude channels without phase partner and phase channels without amplitude partner if the user wants to
        for i in range(len(remaining_channels)):
            data = self.all_data[remaining_channels[i]]
            data_blurred = self._gauss_blurr_data(data, sigma)
            self.all_data[remaining_channels[i]] = data_blurred
            self.channels_label[remaining_channels[i]] = self.channels_label[remaining_channels[i]] + '_' + self.filter_gauss_indicator
        print('Blurring process finished!')
                   
    def _get_compl_angle(self, compl_number_array) -> np.array:
        """This function returns the angles of a clomplex number array.
        
        Args:
            compl_number_array (np.array): complex number array
        """
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

    def _fourier_filter_array(self, complex_array) -> np.array:
        '''
        Takes a complex array and returns the fourier transformed complex array.

        Args:
            complex_array (np.array): complex array to fourier transform
        '''
        FS_compl = np.fft.fftn(complex_array)
        return FS_compl
    
    def fourier_filter_channels(self, channels:list=None) -> None:
        """This function applies the Fourier filter to all data in memory or specified channels.
                
        Args:
            channels [list]: list of channels, will override the already existing channels
        """
        self._initialize_data(channels)
        self._write_to_logfile('fourier_filter', True)
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
            FS_compl = self._fourier_filter_array(compl)
            FS_compl_abs = np.absolute(FS_compl)
            FS_compl_angle = self._get_compl_angle(FS_compl)
            self.all_data[channels_to_filter[i]] = np.log(np.abs(np.fft.fftshift(FS_compl_abs))**2)
            self.channels_label[channels_to_filter[i]] = self.channels_label[channels_to_filter[i]] + '_fft'
            self.all_data[channels_to_filter[i+1]] = FS_compl_angle
            self.channels_label[channels_to_filter[i+1]] = self.channels_label[channels_to_filter[i+1]] + '_fft'

    def fourier_filter_channels_V2(self, channels:list=None) -> None:
        """This function applies the Fourier filter to all data in memory or specified channels
                
        Args:
            channels [list]: list of channels, will override the already existing channels
        """
        self._write_to_logfile('fourier_filter', True)
        if channels is None:
            channels = self.channels
        
        for i in range(len(channels)):
            FS = self._fourier_filter_array(self.all_data[self.channels.index(channels[i])])
            self.all_data[channels[i]] = np.log(np.abs(np.fft.fftshift(FS))**2)
            self.channels_label[channels[i]] = self.channels_label[channels[i]] + '_fft'

    def _create_header(self, channel, data=None, filetype='gsf'):
        """This function creates the header for the gsf file. The header contains all necessary information for the gsf file.

        Args:
            channel (str): channel name
            data (np.array, optional): the data to save, if not specified the data will be loaded from the file. Defaults to None.
            filetype (str, optional): the filetype to save the data. Defaults to 'gsf'.
        """
        # todo XOffset, YOffset dont work properly, also if the measurement is rotated or cut this is not considered so far
        # actually not shure if that isn't fixed by now...
        if data is None:
            # channel is not in memory, so the standard values will be used
            data = self._load_data([channel])[0][0]
            try: XReal, YReal = self._get_measurement_tag_dict_value(MeasurementTags.SCANAREA)
            except: XReal, YReal, ZReal = self._get_measurement_tag_dict_value(MeasurementTags.SCANAREA)
            rotation = self._get_measurement_tag_dict_value(MeasurementTags.ROTATION)[0]
            XOffset, YOffset = self._get_measurement_tag_dict_value(MeasurementTags.SCANNERCENTERPOSITION)
        else: 
            # if channel is in memory it has to have a channel dict, where all necessary infos are stored
            XReal, YReal, *args = self._get_channel_tag_dict_value(channel, ChannelTags.SCANAREA)
            rotation = self._get_channel_tag_dict_value(channel, ChannelTags.ROTATION)[0]
            XOffset, YOffset = self._get_channel_tag_dict_value(channel, ChannelTags.SCANNERCENTERPOSITION)
        if rotation is None:
            # try to get the rotation from the measurement tags
            rotation = self._get_measurement_tag_dict_value(MeasurementTags.ROTATION)[0]
            if rotation is None: rotation = ''
        XRes = len(data[0])
        YRes  = len(data)
        if filetype=='gsf':
            header = f'Gwyddion Simple Field 1.0\n'
        elif filetype=='txt':
            header = 'Simple Textfile, floats seperated by spaces\n'
        else:
            header = ''
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

    def save_to_gsf(self, channels:list=None, appendix:str='default'):
        """This function is ment to save all specified channels to external .gsf files.
        
        Args:
            channels (list, optional):    list of the channels to be saved, if not specified, all channels in memory are saved.
                                Careful! The data will be saved as it is right now, so with all the manipulations.
                                Therefor the data will have an '_manipulated' appendix in the filename.
            appendix (str, optional):     appendix to add to the filename, default is the default specified in the config of the current filetype.
        """
        if appendix == 'default':
            appendix = self.channel_suffix_manipulated
        if channels == None:
            channels = self.channels
        for channel in channels:
            # find out if channel is default or not
            if channel in self.all_channels_default:
                suffix = self.channel_suffix_default
                prefix = self.channel_prefix_default
                channel_type = 'default'
            elif channel in self.all_channels_custom:
                suffix = self.channel_suffix_custom
                prefix = self.channel_prefix_custom
                channel_type = 'custom'
                # ignore the default appendix if the channel is not a default channel 
                if self.channel_suffix_overlain in channel:
                    appendix = ''
                elif self.channel_suffix_synccorrected_phase in channel:
                    appendix = ''
            else:
                print('channel not found in default or custom channels')
                exit()
            filepath = self.directory_name / Path(self.filename.name + f'{prefix}{channel}{suffix}{appendix}.gsf')
            data = self.all_data[self.channels.index(channel)]
            XRes = len(data[0])
            YRes  = len(data)
            header, NUL = self._create_header(channel, data)
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
        self._write_to_logfile('save_to_gsf_appendix', appendix)

    def save_to_txt(self, channels:list=None, appendix:str='default'):
        """This function is ment to save all specified channels to external .txt files.
        
        Args:
            channels (list, optional):    list of the channels to be saved, if not specified, all channels in memory are saved
                                Careful! The data will be saved as it is right now, so with all the manipulations.
                                Therefor the data will have an '_manipulated' appendix in the filename.
            appendix (str, optional):     appendix to add to the filename, default is the default specified in the config of the current filetype.
        """
        if appendix == 'default':
            appendix = self.channel_suffix_manipulated
        if channels == None:
            channels = self.channels
        for channel in channels:
            if channel in self.all_channels_default:
                suffix = self.channel_suffix_default
                prefix = self.channel_prefix_default
                channel_type = 'default'
            elif channel in self.all_channels_custom:
                suffix = self.channel_suffix_custom
                prefix = self.channel_prefix_custom
                channel_type = 'custom'
                # ignore the default appendix if the channel is not a default channel 
                if self.channel_suffix_overlain in channel:
                    appendix = ''
                elif self.channel_suffix_synccorrected_phase in channel:
                    appendix = ''
            else:
                print('channel not found in default or custom channels')
                exit()
            
            filepath = self.directory_name / Path(self.filename.name + f'{prefix}{channel}{suffix}{appendix}.txt')
            data = self.all_data[self.channels.index(channel)]
            XRes = len(data[0])
            YRes  = len(data)
            header, NUL = self._create_header(channel, data, 'txt')
            file = open(filepath, 'w')
            file.write(header)
            # file.write(NUL) # the NUL marks the end of the header and konsists of 0 characters in the first dataline
            for y in range(YRes):
                for x in range(XRes):
                    file.write(f'{round(data[y][x], 5)} ')
            file.close()
            print(f'successfully saved channel {channel} to .txt')
        self._write_to_logfile('save_to_txt_appendix', appendix)
    
    def _create_synccorr_preview(self, channel, wavelength, nouserinput=False) -> None:
        """
        This function is part of the synccorrection and creates a preview of the corrected data.

        Args:
            channel (str): channel to create the preview from
            wavelength (float): wavelength in µm
            nouserinput (bool, optional): if True, the function will not ask for user input. Defaults to False.
        """
        scanangle = self._get_measurement_tag_dict_value(MeasurementTags.ROTATION)[0]*np.pi/180
        phasedir_positive = 1
        phasedir_negative = -1
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
                phase_no_correction[y][x] = phase_data[y][x]
                phase_positive[y][x] = np.mod(phase_data[y][x] - phasedir_positive*(np.cos(-scanangle)*xreal + np.sin(-scanangle)*yreal)/wavelength*2*np.pi, 2*np.pi)
                phase_negative[y][x] = np.mod(phase_data[y][x] - phasedir_negative*(np.cos(-scanangle)*xreal + np.sin(-scanangle)*yreal)/wavelength*2*np.pi, 2*np.pi)
        #create plots of the uncorrected and corrected images
        subplots = []
        subplots.append(self._add_subplot(phase_no_correction, channel))
        subplots.append(self._add_subplot(phase_positive, channel + '_positive'))
        subplots.append(self._add_subplot(phase_negative, channel + '_negative'))
        self._plot_subplots(subplots)
        # remove the preview subplots from the subplot memory after plotting
        self.remove_last_subplots(3)
        #ask the user to chose a correction direction
        if nouserinput is False:
            phasedir = self._gen_from_input_phasedir()
            return phasedir

    def synccorrection(self, wavelength:float, phasedir:int=None) -> None:
        """This function corrects all the phase channels for the linear phase gradient which stems from the synchronized measurement mode.
        The wavelength must be given in µm. The phasedir is either 1 or -1. If you are unshure about the direction just leave the parameter out.
        You will be shown a preview for both directions and then you must choose the correct one.
        The synccorrection will then be applied to all phase channels in memory.
        The corrected channels will then be saved as new files with the synccorrection appendix specified in the config.ini file.
        Afterwards the original channels and data will be reloaded in memory.
                
        Args:
            wavelenght (float): please enter the wavelength in µm.
            phasedir (int, optional): the phase direction, leave out if not known and you will be prompted with a preview and can select the appropriate direction.

        """
        if self.autoscale == True:
            print('careful! The synccorretion does not work when autoscale is enabled.')
            exit()
        # now load all channels in memory for the synccorrection, but save the original data and channels and reinitialize the data lateron
        old_channels = self.channels.copy()
        old_data = self.all_data.copy()
        old_channel_tag_dict = self.channel_tag_dict.copy()
        old_channels_label = self.channels_label.copy()
        old_measurement_tag_dict = self.measurement_tag_dict.copy()
        # load new channels for synccorrection
        all_channels = self.phase_channels + self.amp_channels
        self._initialize_data(all_channels)
        scanangle = self._get_measurement_tag_dict_value(MeasurementTags.ROTATION)[0]*np.pi/180
        if phasedir == None:
            phasedir = self._create_synccorr_preview(self.preview_phasechannel, wavelength)
        self._write_to_logfile('synccorrection_wavelength', wavelength)
        self._write_to_logfile('synccorrection_phasedir', phasedir)
        header, NUL = self._create_header(self.preview_phasechannel) # channel for header just important to distinguish z axis unit either m or nothing
        for channel in self.phase_channels:
            i = self.phase_channels.index(channel)
            phasef = open(self.directory_name / Path(self.filename.name + f' {channel}_corrected.gsf'), 'bw')
            realf = open(self.directory_name / Path(self.filename.name + f' {self.real_channels[i]}_corrected.gsf'), 'bw')
            phasef.write(header.encode('utf-8'))
            realf.write(header.encode('utf-8'))
            phasef.write(NUL) # add NUL terminator
            realf.write(NUL)
            for y in range(0,self.YRes):
                for x in range(0,self.XRes):
                    #convert pixel number to realspace coordinates in µm
                    xreal=x*self.XReal/self.XRes
                    yreal=y*self.YReal/self.YRes
                    #open the phase, add pi to change the range from 0 to 2 pi and then substract the linear phase gradient, which depends on the scanangle!
                    amppixval = self.all_data[self.channels.index(self.amp_channels[i])][y][x]
                    phasepixval = self.all_data[self.channels.index(self.phase_channels[i])][y][x]
                    phasepixval_corr = np.mod(phasepixval + np.pi - phasedir*(np.cos(-scanangle)*xreal + np.sin(-scanangle)*yreal)/wavelength*2*np.pi, 2*np.pi)
                    realpixval = amppixval*np.cos(phasepixval_corr)
                    phasef.write(pack("f",phasepixval_corr))
                    realf.write(pack("f",realpixval))
            phasef.close()
            realf.close()
        # reinitialize the old data
        self.channels = old_channels
        self.all_data = old_data
        self.channel_tag_dict = old_channel_tag_dict
        self.channels_label = old_channels_label
        self.measurement_tag_dict = old_measurement_tag_dict
        gc.collect()

    def _gen_from_input_phasedir(self) -> int:
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
            self._gen_from_input_phasedir()
    
    def _get_channel_scaling(self, channel_id) -> int :
        """This function checks if an instance channel is scaled and returns the scaling factor.
        
        Args:
            channel_id (int): the channel index
        """
        channel_yres = len(self.all_data[channel_id])
        return int(channel_yres/self.YRes)

    def _create_height_mask_preview(self, mask_array) -> None:
        """This function creates a preview of the height masking.
        The preview is based on all channels in the instance
        
        Args:
            mask_array (np.array): the mask array to preview
        """
        channels = self.channels
        dataset = self.all_data
        subplots = []
        for i in range(len(dataset)):
            masked_array = np.multiply(dataset[i], mask_array)
            subplots.append(self._add_subplot(np.copy(masked_array), channels[i]))
        self._plot_subplots(subplots)
        # remove the preview subplots from the memory
        self.remove_last_subplots(3)
        
    def _user_input_bool(self) -> bool: 
        """This function asks the user to input yes or no and returns a boolean value."""
        user_input = input('Please type y for yes or n for no. \nInput: ')
        if user_input == 'y':
            user_bool = True
        elif user_input == 'n':
            user_bool = False
        return user_bool

    def _user_input(self, message:str):
        """This function confronts the user with the specified message and returns the user input

        Args:
            message (str): the message to display
        """
        return input(message)

    def _create_mask_array(self, height_data:np.array, threshold:float) -> np.array:
        """This function takes the height data and a threshold value to create a mask array containing 0 and 1 values.

        Args:
            height_data (np.array): the height data
            threshold (float): the threshold value
        
        Returns:
            np.array: the mask array
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

    def _get_height_treshold(self, height_data:np.array) -> float:
        """This function returns the height threshold value. The user is prompted with a preview of the mask array and can adjust the threshold using a slider.
        
        Args:
            height_data (np.array): the height data
        
        Returns:
            float: the new threshold [0-1]
        """
        threshold = round(get_height_treshold(height_data), 2)
        '''self._create_height_mask_preview(mask_array)
        print('Do you want to use these parameters to mask the data?')
        mask_data = self._user_input_bool()
        if mask_data == False:
            print('Do you want to change the treshold?')
            change_treshold = self._user_input_bool()
            if change_treshold == True:
                print(f'The old threshold was {threshold}')
                threshold = float(input('Please enter the new treshold value: '))
                mask_array = self._create_mask_array(height_data, threshold)
                self._get_height_treshold(height_data, mask_array, threshold)
            else:
                print('Do you want to abort the masking procedure?')
                abort = self._user_input_bool()
                if abort == True:
                    exit()'''
        return threshold

    def heigth_mask_channels(self, channels:list=None, mask_channel:str=None, threshold:float=None) -> None:
        """
        The treshold factor should be between 0 and 1. It sets the threshold for the height pixels.
        Every pixel below threshold will be set to 0. This also applies for all other channels. 
        You can either specify specific channels to mask or if you don't specify channels,
        all standard channels will be masked. If export is False only the channels in self.channels will be masked
        and nothing will be exported. 
        For this function to also work with scaled data the height channel has to be specified and scaled as well!
                
        Args:
            channels (list): list of channels, will override the already existing channels
            mask_channel (str): The channel to use for the height mask, if not specified the height channel will be used
            threshold (float): Threshold value to create the height mask from. Default is None, the user can select the threshold with a slider.
        """
        if channels == None:
            channels = self.channels
        if (mask_channel == None) or (mask_channel not in self.channels):
            if self.height_channel in self.channels:
                mask_channel = self.height_channel
            else:
                print('Please specify a mask channel!')
                exit()
        if self.height_indicator not in mask_channel:
            print('Please specify a height channel!')
            exit()
        else:
            height_data = self.all_data[self.channels.index(mask_channel)]                

        if threshold is None:
            threshold = self._get_height_treshold(height_data)

        mask_array = self._create_mask_array(height_data, threshold)
        self.mask_array = mask_array # todo, mask array must be saved as part of the image, otherwise multiple measurement creations will use the same mask

        self._write_to_logfile('height_masking_threshold', threshold)
        for channel in channels:
            if channel not in self.channels:
                print(f'Channel {channel} is not in memory! Please initiate the channels you want to use first!')
            self.all_data[self.channels.index(channel)] = np.multiply(self.all_data[self.channels.index(channel)], mask_array)
            self.channels_label[self.channels.index(channel)] = self.channels_label[self.channels.index(channel)] + '_masked'
            
        # dataset = self.all_data
        # for i in range(len(dataset)):
        #     if self.height_channel not in self.channels_label[i]:
        #         self.all_data[i] = np.multiply(dataset[i], mask_array)
        #     self.channels_label[i] = self.channels_label[i] + '_masked'
        print('Channels in memory have been masked!')

    def _check_pixel_position(self, xres:int, yres:int, x:int, y:int) -> bool:
        """This function checks if the pixel position is within the bounds.
        
        Args:
            xres (int): x resolution
            yres (int): y resolution
            x (int): x coordinate
            y (int): y coordinate
        
        Returns:
            bool: True if the pixel position is within the bounds, False otherwise
        """
        if x < 0 or x > xres:
            return False
        elif y < 0 or y > yres:
            return False
        else: return True

    def _get_mean_value(self, data:np.array, x_coord:int, y_coord:int, zone:int) -> float:
        """This function returns the mean value of the pixel and its nearest neighbors.
        The zone specifies the number of neighbors. 1 means the pixel and the 8 nearest pixels.
        2 means zone 1 plus the next 16, so a total of 25 with the pixel in the middle. 

        Args:
            data (np.array): the data
            x_coord (int): x coordinate
            y_coord (int): y coordinate
            zone (int): the number of neighbors

        Returns:
            float: the mean value
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
                if self._check_pixel_position(xres, yres, x_pixel, y_pixel) == True:
                    mean += data[y_pixel][x_pixel]
                    count += 1
        return mean/count

    def get_pixel_coordinates(self, channel) -> list:
        """This function returns the pixel coordinates of the clicked pixel.
        
        Args:
            channel (str): the channel to display
            
        Returns:
            list: the pixel coordinates
        """
        data = self.all_data[self.channels.index(channel)]
        # identify the colormap
        if self.height_indicator in channel:
            cmap = SNOM_height
        elif self.phase_indicator in channel:
            cmap = SNOM_phase
        elif self.amp_indicator in channel:
            cmap = SNOM_amplitude
        else:
            cmap = 'viridis'
        fig, ax = plt.subplots()
        ax.pcolormesh(data, cmap=cmap)
        klicker = clicker(ax, ["event"], markers=["x"])
        ax.legend()
        ax.axis('scaled')
        ax.invert_yaxis()
        plt.title('Please click on the pixel you want to get the coordinates from.')
        if PlotDefinitions.show_plot:
            plt.show()
        klicker_coords = klicker.get_positions()['event'] #klicker returns a dictionary for the events
        coordinates = [[round(element[0]), round(element[1])] for element in klicker_coords]
        # display image with the clicked pixel
        fig, ax = plt.subplots()
        ax.pcolormesh(data, cmap=cmap)
        ax.plot(coordinates[0][0], coordinates[0][1], 'rx')
        ax.legend()
        ax.axis('scaled')
        ax.invert_yaxis()
        plt.title('You clicked on the following pixel.')
        if PlotDefinitions.show_plot:
            plt.show()
        return coordinates

    def get_pixel_value(self, channel, coordinates:list=None, zone:int=1) -> float:
        """This function returns the pixel value of a channel at the specified coordinates.
        The zone specifies the number of neighbors. 0 means only the pixel itself. 1 means the pixel and the 8 nearest pixels.
        2 means zone 1 plus the next 16, so a total of 25 with the pixel in the middle.
        If the channel is scaled the zone will be scaled as well.
        
        Args:
            channel (str): the channel to display
            coordinates (list, optional): the pixel coordinates. Defaults to None.
            zone (int, optional): the number of neighbors. Defaults to 1.
        
        Returns:
            float: the pixel value
        """
        # adjust the zone if the data is scaled
        zone = zone*self._get_channel_scaling(self.channels.index(channel))
        # display the channel
        data = self.all_data[self.channels.index(channel)]
        if coordinates == None:
            coordinates = self.get_pixel_coordinates(channel)
        if len(coordinates) != 1:
            print('You need to specify one pixel coordinate! \nDo you want to try again?')
            user_input = self._user_input_bool()
            if user_input == True:
                self.get_pixel_value(channel, zone)
            else:
                exit()
        x = coordinates[0][0]
        y = coordinates[0][1]
        # get the mean value of the pixel and its neighbors
        pixel_value = self._get_mean_value(data, x, y, zone)
        return pixel_value

    def _height_levelling_3point(self, height_data:np.array, coords:list=None, zone:int=1) -> np.array:
        """This function levels the height data with a 3 point leveling.
        The user has to click on three points to specify the underground plane.
        The function returns the leveled height data.
        
        Args:
            height_data (np.array): the height data
            zone (int, optional): the number of neighbors. Defaults to 1.

        Returns:
            np.array: the leveled height data
        """
        # check if coordinates are given, then we don't need to display the image
        if coords == None:
            fig, ax = plt.subplots()
            ax.pcolormesh(height_data, cmap=SNOM_height)
            klicker = clicker(ax, ["event"], markers=["x"])
            ax.legend()
            ax.axis('scaled')
            plt.title('3 Point leveling: please click on three points\nto specify the underground plane.')
            if PlotDefinitions.show_plot:
                plt.show()
            klicker_coords = klicker.get_positions()['event'] #klicker returns a dictionary for the events
            coords = [[round(element[0]), round(element[1])] for element in klicker_coords]
        if len(coords) != 3:
            print('You need to specify 3 point coordinates! \nDo you want to try again?')
            user_input = self._user_input_bool()
            if user_input == True:
                self._height_levelling_3point(height_data, coords, zone)
            else:
                exit()
        self._write_to_logfile('height_leveling_coordinates', coords)
        # for the 3 point coordinates the height data is calculated over a small area around the clicked pixels to reduce deviations due to noise
        mean_values = [self._get_mean_value(height_data, coords[i][0], coords[i][1], zone) for i in range(len(coords))]
        matrix = [[coords[i][0], coords[i][1], mean_values[i]] for i in range(3)]
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
    
    def _level_height_data(self, height_data:np.array, klick_coordinates:list, zone:int):
        """This function levels the height data with a 3 point leveling.
        The user has to click on three points to specify the underground plane.
        The function returns the leveled height data. This version is just for the gui.

        Args:
            height_data (np.array): the height data
            klick_coordinates (list): the pixel coordinates
            zone (int): the number of neighbors

        Returns:
            np.array: the leveled height data
        """

        mean_values = [self._get_mean_value(height_data, klick_coordinates[i][0], klick_coordinates[i][1], zone) for i in range(len(klick_coordinates))]
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

    def _get_klicker_coordinates(data, cmap, message:str):
        """This function returns the pixel coordinates of the clicked pixel.

        Args:
            data (np.array): the data
            cmap (str): the colormap
            message (str): the message to display as the title
        """
        fig, ax = plt.subplots()
        ax.pcolormesh(data, cmap=cmap)
        klicker = clicker(ax, ["event"], markers=["x"])
        ax.legend()
        ax.axis('scaled')
        plt.title(message)
        plt.show()
        klicker_coords = klicker.get_positions()['event'] #klicker returns a dictionary for the events
        klick_coordinates = [[round(element[0]), round(element[1])] for element in klicker_coords]
        return klick_coordinates

    def _height_levelling_3point_forGui(self, height_data, zone=1) -> np.array:
        klick_coordinates = self._get_klicker_coordinates(height_data, SNOM_height, '3 Point leveling: please click on three points\nto specify the underground plane.')
        if len(klick_coordinates) != 3:
            print('You need to specify 3 point coordinates! Data was not leveled!')
            return height_data
        # for the 3 point coordinates the height data is calculated over a small area around the clicked pixels to reduce deviations due to noise
        self._write_to_logfile('height_leveling_coordinates', klick_coordinates)
        return self._level_height_data(klick_coordinates, zone)

    def _level_phase_slope(self, data:np.array, slope:float) -> np.array:
        """This function substracts a linear phase gradient in y direction from the specified phase data.
        The data is then also shifted by 0 to ensure that the phase data is still in the range of 0 to 2pi.

        Args:
            data (np.array): the phase data
            slope (float): the slope

        Returns:
            np.array: the leveled phase data
        """
        yres = len(data)
        xres = len(data[0])
        for y in range(yres):
            for x in range(xres):
                data[y][x] -= y*slope
        return self._shift_phase_data(data, 0)

    def correct_phase_drift(self, channels:list=None, export:bool=False, phase_slope:float=None, zone:int=1) -> None:
        """This function asks the user to click on two points which should have the same phase value.
        Only the slow drift in y-direction will be compensated. Could in future be extended to include a percentual drift compensation along the x-direction.
        But should usually not be necessary.
                
        Args:
            channels (list, optional): list of channels, will override the already existing channels
            export (bool, optional): do you want to aply the correction to all phase channels and export them?
            phase_slope (float, optional): if you already now the phase slope you can enter it, otherwise leave it out
                                and it will prompt you with a preview to select two points to calculate the slope from
            zone (int, optional): defines the area which is used to calculate the mean around the click position in the preview,
                        0 means only the click position, 1 means the nearest 9 ...
        """
        self._initialize_data(channels)
        phase_data = None
        if self.preview_phasechannel in self.channels:
            phase_data = np.copy(self.all_data[self.channels.index(self.preview_phasechannel)])
            phase_channel = self.preview_phasechannel
        else:
            phase_data = self._load_data([self.preview_phasechannel])[0][0]
            phase_channel = self.preview_phasechannel
        if export == True:
            # ToDo
            # do something with the phase slope...
            print('You want to export a phase slope correction, but nothing happens!')
            pass
        else:
            if phase_slope != None:
                #level all phase channels in memory...
                self._write_to_logfile('phase_driftcomp_slope', phase_slope)
                for i in range(len(self.channels)):
                    if 'P' in self.channels[i]:
                        self.all_data[i] = self._level_phase_slope(self.all_data[i], phase_slope)
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
                    user_input = self._user_input_bool()
                    if user_input == True:
                        self.correct_phase_drift(channels, export, None)
                    else: 
                        exit()
                mean_values = [self._get_mean_value(phase_data, klick_coordinates[i][0], klick_coordinates[i][1], zone) for i in range(len(klick_coordinates))]
                #order points from top to bottom
                if klick_coordinates[0][1] > klick_coordinates[1][1]:
                    second_corrd = klick_coordinates[0]
                    second_mean = mean_values[0]
                    klick_coordinates[0] = klick_coordinates[1]
                    klick_coordinates[1] = second_corrd
                    mean_values[0] = mean_values[1]
                    mean_values[1] = second_mean
                phase_slope = (mean_values[1] - mean_values[0])/(klick_coordinates[1][1] - klick_coordinates[0][1])
                leveled_phase_data = self._level_phase_slope(phase_data, phase_slope)
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
                user_input = self._user_input_bool()
                if user_input == True:
                    #use the phase slope to level all phase channels in memory
                    self.correct_phase_drift(None, False, phase_slope)
                else:
                    print('Do you want to repeat the leveling?')
                    user_input = self._user_input_bool()
                    if user_input == True:
                        #start the leveling process again
                        self.correct_phase_drift()
                    else:
                        exit()
        gc.collect()

    def correct_phase_drift_nonlinear(self, channels:list=None, reference_area:list = [None, None]) -> None:
        """This function corrects the phase drift in the y-direction by using a reference area across the full length of the scan.	
        The reference area is used to calculate the average phase value per row.
        This value is then substracted from the phase data to level the phase.
        The reference area is specified by two coordinates, the left and right border. If no area is specified the whole image will be used.
        Make shure not to rotate the image prior to this function, since the reference area is defined in y-direction.
        This function is somewhat redundant to the level_data_columnwise function, which works for all channels (amplitude, height and phase).

        Args:
            channels (list, optional): list of channels, will override the already existing channels
            reference_area (list, optional): The reference area to calculate the phase offset, specify as reference_area=[left-border, right-border].
                If not specified the whole image will be used. Defaults to [None, None].
        """

        # if a list of channels is specified those will be loaded and the old ones will be overwritten
        self._initialize_data(channels)
        # define local list of channels to use for leveling
        channels = self.channels
        phase_data = None
        if self.preview_phasechannel in self.channels:
            phase_data = np.copy(self.all_data[self.channels.index(self.preview_phasechannel)])
            phase_channel = self.preview_phasechannel
        else:
            phase_data = self._load_data([self.preview_phasechannel])[0][0]
            phase_channel = self.preview_phasechannel
        
        # cut out the reference area
        # if no area is specified just use the whole data
        if reference_area[0] == None:
            reference_area[0] = 0 # left border
        if reference_area[1] == None:
            reference_area[1] = len(phase_data[0]) # right border

        # get the phase values per column of the reference area, then flatten each column 
        flattened_phase_profiles = []
        for j in range(reference_area[0], reference_area[1]):
            reference_values = [phase_data[i][j] for i in range(len(phase_data))]
            reference_values_flattened = phase_analysis.flatten_phase_profile(reference_values, 1)
            # reference_values_flattened = np.unwrap(reference_values)
            flattened_phase_profiles.append(reference_values_flattened)

        # average all flattened profiles
        reference_values_flattened = np.mean(flattened_phase_profiles, axis=0)

        # remove the averaged reference data per line from the phase data
        leveled_phase_data = np.copy(phase_data)
        for i in range(len(phase_data)):
            leveled_phase_data[i] = (leveled_phase_data[i] - reference_values_flattened[i] + np.pi) %(2*np.pi)

        # display the leveled phase data
        fig, ax = plt.subplots()
        img = ax.pcolormesh(leveled_phase_data, cmap=SNOM_phase)
        ax.invert_yaxis()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="10%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('phase', rotation=270)
        # ax.legend()
        ax.axis('scaled')
        plt.title('Leveled Pase: ' + phase_channel)
        plt.show()

        print('Are you satisfied with the phase leveling?')
        user_input = self._user_input_bool()
        if user_input == True:
            # write to logfile
            self._write_to_logfile('phase_driftcomp_nonlinear_reference_area', reference_area)
            # do the leveling for all channels but use always the same reference data, channels should only differ in phase offset
            for i in range(len(channels)):
                if 'P' in channels[i]:
                    self.all_data[self.channels.index(channels[i])] = np.array([(self.all_data[self.channels.index(channels[i])][j] - reference_values_flattened[j] + np.pi) %(2*np.pi) for j in range(len(reference_values_flattened))])
                    # also apply a phase shift to ensure that the phase is between 0 and 2pi
                    # for now take the average phase an shift it to pi/2 should be white on the colormap
                    phase_shift = np.pi/2 - np.mean(self.all_data[self.channels.index(channels[i])])
                    self.all_data[self.channels.index(channels[i])] = self._shift_phase_data(self.all_data[self.channels.index(channels[i])], phase_shift)
        gc.collect()

    def match_phase_offset(self, channels:list=None, reference_channel:str=None, reference_area=None, manual_width=5) -> None:
        """This function matches the phase offset of all phase channels in memory to the reference channel.
        The reference channel is the first phase channel in memory if not specified.

        Args:
            channels (list, optional): list of channels, will override the already existing channels
            reference_channel (str, optional): The reference channel to which all other phase channels will be matched.
                If not specified the first phase channel in memory will be used. Defaults to None.
            reference_area (list or str, optional): The area in the reference channel which will be used to calculate the phase offset. If not specified the whole image will be used.
                You can also specify 'manual' then you will be asked to click on a point in the image. The area around that pixel will then be used as reference
                You can also specify a list like in the logfile to use a specific area. Defaults to None.
            manual_width (int, optional): The width of the manual reference area. Only applies if reference_area='manual'. Defaults to 5.
        """
        # if a list of channels is specified those will be loaded and the old ones will be overwritten
        self._initialize_data(channels)
        # define local list of channels to use for leveling
        channels = self.channels
        if reference_channel == None:
            for channel in channels:
                if self.phase_indicator in channel:
                    reference_channel = channel
                    break
        if reference_area is None:
            # reference_area = [[xmin, xmax][ymin, ymax]]
            reference_area = [[0, len(self.all_data[self.channels.index(reference_channel)][0])],[0, len(self.all_data[self.channels.index(reference_channel)])]]
        elif reference_area == 'manual':
            # use pointcklicker to get the reference area
            fig, ax = plt.subplots()
            ax.pcolormesh(self.all_data[self.channels.index(reference_channel)], cmap=SNOM_phase)
            klicker = clicker(ax, ["event"], markers=["x"])
            ax.legend()
            ax.axis('scaled')
            ax.invert_yaxis()
            plt.title('Please click in the area to use as reference.')
            plt.show()
            klicker_coords = klicker.get_positions()['event']
            klick_coordinates = [[round(element[0]), round(element[1])] for element in klicker_coords]
            # make sure only one point is selected
            if len(klick_coordinates) != 1 and type(klick_coordinates[0]) != list:
                print('You must specify one point which should define the reference area!')
                print('Do you want to try again?')
                user_input = self._user_input_bool()
                if user_input == True:
                    self.match_phase_offset(channels, reference_channel, 'manual')
                else:
                    exit()
            reference_area = [[klick_coordinates[0][0] - manual_width,klick_coordinates[0][0] + manual_width],[klick_coordinates[0][1] - manual_width, klick_coordinates[0][1] + manual_width]]
        
        reference_data = self.all_data[self.channels.index(reference_channel)]
        reference_phase = np.mean([reference_data[i][reference_area[0][0]:reference_area[0][1]] for i in range(reference_area[1][0], reference_area[1][1])])
        
        # display the reference area
        fig, ax = plt.subplots()
        img = ax.pcolormesh(reference_data, cmap=SNOM_phase)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('phase', rotation=270)
        # ax.legend()
        ax.axis('scaled')  
        rect = patches.Rectangle((reference_area[0][0], reference_area[1][0]), reference_area[0][1]-reference_area[0][0], reference_area[1][1]-reference_area[1][0], linewidth=1, edgecolor='g', facecolor='none')
        ax.add_patch(rect)
        ax.invert_yaxis()
        plt.title('Reference Area: ' + reference_channel)
        plt.show()

        for channel in channels:
            if self.phase_indicator in channel:
                phase_data = self.all_data[self.channels.index(channel)]
                # phase_offset = np.mean(phase_data) - reference_phase
                phase_offset = np.mean([phase_data[i][reference_area[0][0]:reference_area[0][1]] for i in range(reference_area[1][0], reference_area[1][1])]) - reference_phase
                self.all_data[self.channels.index(channel)] = self._shift_phase_data(phase_data, -phase_offset)
        self._write_to_logfile('match_phase_offset_reference_area', reference_area)
        gc.collect()

    def correct_amplitude_drift_nonlinear(self, channels:list=None, reference_area:list = [None, None]) -> None:
        """This function corrects the amplitude drift in the y-direction by using a reference area across the full length of the scan.	
        The reference area is used to calculate the average amplitude value per row.
        This value is then divided from the amplitude data to level the amplitude.
        The reference area is specified by two coordinates, the left and right border. If no area is specified the whole image will be used.
        Make shure not to rotate the image prior to this function, since the reference area is defined in y-direction.
        This function is somewhat redundant to the level_data_columnwise function, which works for all channels (amplitude, height and phase).

        Args:
            channels (list, optional): list of channels, will override the already existing channels
            reference_area (list, optional): The reference area to calculate the amplitude offset, specify as reference_area=[left-border, right-border].
                If not specified the whole image will be used. Defaults to [None, None].
        """

        # if a list of channels is specified those will be loaded and the old ones will be overwritten
        self._initialize_data(channels)
        # define local list of channels to use for leveling
        channels = self.channels
        amplitude_data = None
        if self.preview_ampchannel in self.channels:
            amplitude_data = np.copy(self.all_data[self.channels.index(self.preview_ampchannel)])
            amplitude_channel = self.preview_ampchannel
        else:
            amplitude_data = self._load_data([self.preview_ampchannel])[0][0]
            amplitude_channel = self.preview_ampchannel
        
        # cut out the reference area
        # if no area is specified just use the whole data
        if reference_area[0] == None:
            reference_area[0] = 0
        if reference_area[1] == None:
            reference_area[1] = len(amplitude_data[0])
        
        # iterate through the reference area and get the average amplitude value per row
        reference_values = [np.mean(amplitude_data[i][reference_area[0]:reference_area[1]]) for i in range(len(amplitude_data))]

        # we assume the average amplitude should stay constant, so we divide the amplitude data by the reference values and multiply by the mean reference value
        leveled_amplitude_data = np.copy(amplitude_data)
        for i in range(len(amplitude_data)):
            leveled_amplitude_data[i] = amplitude_data[i] / reference_values[i] * np.mean(reference_values)
        
        # display the original data besides the leveled amplitude data
        fig, ax = plt.subplots(1, 2)
        img1 = ax[0].pcolormesh(amplitude_data, cmap=SNOM_amplitude)
        img2 = ax[1].pcolormesh(leveled_amplitude_data, cmap=SNOM_amplitude)
        ax[0].invert_yaxis()
        ax[1].invert_yaxis()
        divider = make_axes_locatable(ax[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img1, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('amplitude', rotation=270)
        divider = make_axes_locatable(ax[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img2, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('amplitude', rotation=270)
        # ax[0].legend()
        # ax[1].legend()
        ax[0].axis('scaled')
        ax[1].axis('scaled')
        ax[0].set_title('Original Amplitude: ' + amplitude_channel)
        ax[1].set_title('Leveled Amplitude: ' + amplitude_channel)
        plt.show()

        # ask the user if he is satisfied with the leveling
        print('Are you satisfied with the amplitude leveling?')
        user_input = self._user_input_bool()
        if user_input == True:
            # do the leveling for all channels, each channel should be referenced to itself since the amplitudes of the channels will be different
            for i in range(len(channels)):
                if self.amp_indicator in channels[i]:
                    # self.all_data[self.channels.index(channels[i])] = np.copy(self.all_data[self.channels.index(channels[i])])
                    reference_values = [np.mean(self.all_data[self.channels.index(channels[i])][j][reference_area[0]:reference_area[1]]) for j in range(len(self.all_data[self.channels.index(channels[i])]))]
                    self.all_data[self.channels.index(channels[i])] = [(self.all_data[self.channels.index(channels[i])][j] / reference_values[j] * np.mean(reference_values)) for j in range(len(reference_values))]
        else:
            print('Do you want to repeat the leveling?')
            user_input = self._user_input_bool()
            if user_input == True:
                # write to logfile
                self._write_to_logfile('amplitude_driftcomp_nonlinear_reference_area', reference_area)
                #start the leveling process again
                self.correct_amplitude_drift_nonlinear(channels, reference_area)
            else:
                exit()
        gc.collect()

    def correct_height_drift_nonlinear(self, channels:list=None, reference_area:list = [None, None]) -> None:
        """This function corrects the height drift in the y-direction by using a reference area across the full length of the scan.	
        The reference area is used to calculate the average height value per row.
        This value is then divided from the height data to level the height.
        The reference area is specified by two coordinates, the left and right border. If no area is specified the whole image will be used.
        Make shure not to rotate the image prior to this function, since the reference area is defined in y-direction.
        This function is somewhat redundant to the level_data_columnwise function, which works for all channels (amplitude, height and phase).

        Args:
            channels (list, optional): list of channels, will override the already existing channels
            reference_area (list, optional): The reference area to calculate the height offset, specify as reference_area=[left-border, right-border].
                If not specified the whole image will be used. Defaults to [None, None].
        """

        # zone = int(zone*self.scaling_factor/4) #automatically enlargen the zone if the data has been scaled by more than a factor of 4
        # if a list of channels is specified those will be loaded and the old ones will be overwritten
        self._initialize_data(channels)
        # define local list of channels to use for leveling
        channels = self.channels
        height_data = None
        if self.height_channel in self.channels:
            height_data = np.copy(self.all_data[self.channels.index(self.height_channel)])
            height_channel = self.height_channel
        else:
            height_data = self._load_data([self.height_channel])[0][0]
            height_channel = self.height_channel
        
        # cut out the reference area
        # new version: let the user specify the reference area by moving two borders in the preview
        # if no area is specified just use the whole data
        if reference_area[0] == None:
            reference_area[0] = 0
        if reference_area[1] == None:
            reference_area[1] = len(height_data[0])
        
        # iterate through the reference area and get the average height value per row
        reference_values = [np.mean(height_data[i][reference_area[0]:reference_area[1]]) for i in range(len(height_data))]

        # we assume the average height should stay constant, so we divide the height data by the reference values and multiply by the mean reference value
        leveled_height_data = np.copy(height_data)
        for i in range(len(height_data)):
            leveled_height_data[i] = height_data[i] / reference_values[i] * np.mean(reference_values)
        
        # display the original data besides the leveled height data
        fig, ax = plt.subplots(1, 2)
        img1 = ax[0].pcolormesh(height_data, cmap=SNOM_height)
        img2 = ax[1].pcolormesh(leveled_height_data, cmap=SNOM_height)
        ax[0].invert_yaxis()
        ax[1].invert_yaxis()
        divider = make_axes_locatable(ax[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img1, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('height', rotation=270)
        divider = make_axes_locatable(ax[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img2, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('height', rotation=270)
        # ax[0].legend()
        # ax[1].legend()
        ax[0].axis('scaled')
        ax[1].axis('scaled')
        ax[0].set_title('Original height: ' + height_channel)
        ax[1].set_title('Leveled height: ' + height_channel)
        plt.show()

        # ask the user if he is satisfied with the leveling
        print('Are you satisfied with the height leveling?')
        user_input = self._user_input_bool()
        if user_input == True:
            # do the leveling for all channels, each channel should be referenced to itself since the heights of the channels will be different
            for i in range(len(channels)):
                if self.height_indicator in channels[i]:
                    # self.all_data[self.channels.index(channels[i])] = np.copy(self.all_data[self.channels.index(channels[i])])
                    reference_values = [np.mean(self.all_data[self.channels.index(channels[i])][j][reference_area[0]:reference_area[1]]) for j in range(len(self.all_data[self.channels.index(channels[i])]))]
                    self.all_data[self.channels.index(channels[i])] = [(self.all_data[self.channels.index(channels[i])][j] / reference_values[j] * np.mean(reference_values)) for j in range(len(reference_values))]
        else:
            print('Do you want to repeat the leveling?')
            user_input = self._user_input_bool()
            if user_input == True:
                # write to logfile
                self._write_to_logfile('height_driftcomp_nonlinear_reference_area', reference_area)
                #start the leveling process again
                self.correct_height_drift_nonlinear(channels, reference_area)
            else:
                exit()
        gc.collect()

    def level_height_channels_3point(self, channels:list=None, coords:list=None) -> None:
        """This function levels all height channels which are either user specified or in the instance memory.
        The leveling will prompt the user with a preview to select 3 points for getting the coordinates of the leveling plane.
        
        Args:
            channels (list, optional): List of channels to level. If not specified all channels in memory will be used. Defaults to None.
            coords (list, optional): List of coordinates to use for the leveling. If not specified the user will be prompted to click on the points. Defaults to None.
        """
        if channels is None:
            channels = self.channels
        for channel in channels:
            if channel in self.channels and self.height_indicator in channel:
                self.all_data[self.channels.index(channel)] = self._height_levelling_3point(self.all_data[self.channels.index(channel)], coords)
                self.channels_label[self.channels.index(channel)] += '_leveled' 
        gc.collect()

    def level_height_channels_forGui(self, channels:list=None):# todo not used?
        """This function levels all height channels which are either user specified or in the instance memory.
        The leveling will prompt the user with a preview to select 3 points for getting the coordinates of the leveling plane.
        This function is specifically for use with GUI.
        
        Args:
            channels (list, optional): List of channels to level. If not specified all channels in memory will be used. Defaults to None.
        """
        if channels is None:
            channels = self.channels
        for channel in channels:
            if self.height_indicator in channel:
                self.all_data[self.channels.index(channel)] = self._height_levelling_3point_forGui(self.all_data[self.channels.index(channel)])
                self.channels_label[self.channels.index(channel)] += '_leveled' 
        gc.collect()

    def _shift_phase_data(self, data, shift) -> np.array:
        """This function adds a phaseshift to the specified phase data. The phase data is automatically kept in the 0 to 2 pi range.
        Could in future be extended to show a live view of the phase data while it can be modified by a slider...
        e.g. by shifting the colorscale in the preview rather than the actual data..."""
        yres = len(data)
        xres = len(data[0])
        for y in range(yres):
            for x in range(xres):
                data[y][x] = (data[y][x] + shift) % (2*np.pi)
        return data

    def shift_phase(self, shift:float=None, channels:list=None) -> None:
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
        if shift == None:
            shift_known = False
        else:
            shift_known = True
        if shift_known is False:
            if self.preview_phasechannel in channels:
                    phase_data = np.copy(self.all_data[self.channels.index(self.preview_phasechannel)])
            else:
                # check if corrected phase channel is present
                # just take the first phase channel in memory
                for channel in channels:
                    if self.phase_indicator in channel:
                        phase_data = np.copy(self.all_data[self.channels.index(channel)])
                        break
            shift = get_phase_offset(phase_data)
            print('The phase shift you chose is:', shift)
            shift_known = True

        # export shift value to logfile
        self._write_to_logfile('phase_shift', shift)
        # shift all phase channels in memory
        # could also be implemented to shift each channel individually...
        
        for channel in channels:
            # print(channel)
            if self.phase_indicator in channel:
                # print('Before phase shift: ', channel)
                # print('Min phase value:', np.min(self.all_data[self.channels.index(channel)]))
                # print('Max phase value:', np.max(self.all_data[self.channels.index(channel)]))
                self.all_data[self.channels.index(channel)] = self._shift_phase_data(self.all_data[self.channels.index(channel)], shift)
                # print('After phase shift: ', channel)
                # print('Min phase value:', np.min(self.all_data[self.channels.index(channel)]))
                # print('Max phase value:', np.max(self.all_data[self.channels.index(channel)]))
        gc.collect()

    def _fit_horizontal_wg(self, data):
        YRes = len(data)
        XRes = len(data[0])
        #just calculate the shift for each pixel for now
        number_align_points = XRes #the number of intersections fitted with gaussian to find waveguide center along the x direction
        align_points = np.arange(0, XRes, int((XRes)/number_align_points), int)
        cutline_data_sets = []
        for element in align_points:
            cutline = []
            for i in range(YRes):
                cutline.append(data[i][element])
            cutline_data_sets.append(cutline)
        list_of_coefficients = []
        p0 = [100, (YRes)/2, 5, 0]
        bounds = ([0, -YRes, 0, -1000], [1000, YRes, YRes/2, 1000])
        for cutline in cutline_data_sets:
            coeff, var_matrix = curve_fit(gauss_function, range(0, YRes), cutline, p0=p0, bounds=bounds)
            list_of_coefficients.append(coeff)
            p0 = coeff #set the starting parameters for the next fit
        return align_points, list_of_coefficients

    def _shift_data(self, data, axis, shifts) -> np.array:
        # if shifts are not int round them
        if not all(isinstance(n, int) for n in shifts):
            shifts = [round(element) for element in shifts]
        YRes = len(data)
        XRes = len(data[0])
        min_shift = round(min(shifts))
        max_shift = round(max(shifts))
        if axis == 1:
            new_YRes = YRes + int(abs(min_shift-max_shift))
            data_shifted = np.zeros((new_YRes, XRes))
            #create the realigned height
            for x in range(XRes):
                shift = int(-shifts[x] + abs(max_shift)) #the calculated shift has to be compensated by shifting the pixels
                # shift = round(-shifts[x] + abs(max_shift)) #the calculated shift has to be compensated by shifting the pixels
                for y in range(YRes):
                    data_shifted[y + shift][x] = data[y][x]
        elif axis == 0:
            YRes = len(data)
            XRes = len(data[0])
            min_shift = round(min(shifts))
            max_shift = round(max(shifts))
            new_XRes = XRes + int(abs(min_shift-max_shift))
            data_shifted = np.zeros((YRes, new_XRes))
            #create the realigned height
            for y in range(YRes):
                shift = int(-shifts[y] + abs(max_shift))
                # shift = round(-shifts[y] + abs(max_shift))
                for x in range(XRes):
                    data_shifted[y][x + shift] = data[y][x]
        return data_shifted

    def _get_mean_from_area(self, data, axis=1, threshold=0.5):
        """This function calculates the mean index of an array along a specified axis.
        The mean index is calculated by setting all values below a certain threshold to zero.

        Args:
            data (np.array): 2d array of data.
            axis (int): The axis along which the mean index should be calculated. 0 means x-axis, 1 means y-axis. Defaults to 1.
            threshold (float, optional): threshold, all values below will be set to zero to better estimate the mean index position. Defaults to 0.5.

        Returns:
            float: np.array of the mean position indices.
        """
        if axis == 1:
            res = len(data[0])
            sliced_data = [data[:,i] for i in range(res)]
        elif axis == 0:
            res = len(data)
            sliced_data = [data[i] for i in range(res)]
        #just calculate the shift for each pixel column for now
        # number_align_points = XRes
        shifts = np.zeros(res)
        for i in range(res):
            max_val = np.max(sliced_data[i])
            # set all values below threshold to zero
            sliced_data[i] = np.where(sliced_data[i] < threshold*max_val, 0, sliced_data[i])
            mean_index = mean_index_array(sliced_data[i])
            # plot the column data
            # if i%100 == 0:
            #     print('mean index:', mean_index)
            #     plt.plot(column_data)
            #     plt.vlines(mean_index, ymin=min(column_data), ymax=max(column_data), color='red')
            #     plt.show()
            shifts[i] = mean_index
        return shifts

    def realign(self, channels:list=None, bounds:list=None, axis=1, threshold=0.5):
        """This function corrects the drift of the piezo motor. As of now it needs a reference region of the sample which is assumed to be straight.
        In the future this could be implemented with a general map containing the distortion created by the piezo motor, if it turns out to be temporally constant...
        Anyways, you will be prompted with a preview of the height data, please select an area of the scan with only one 'straight' reference. 
        It will then calculate the index of the mean according to the specified axis. If you specify a threshold all values below this threshold will be set to zero.
        This makes the mean index calculation more robust.
        The bounds for the fitting routine are based on the lower and upper limit of this selection.
        
        Args:  
            channels (list): list of channels, will override the already existing channels
            bounds (list): The bounds for the fitting routine. If not specified you will be prompted with a window to select an area.
                Should be specified like this: [lower_bound, upper_bound] in px.
            axis (int): The axis along which the mean index should be calculated. 0 means x-axis, 1 means y-axis. Defaults to 1.
            threshold (float, optional): threshold, all values below will be set to zero to better estimate the mean index position. Defaults to 0.5.
        
        """
        self._initialize_data(channels)
        # store the bounds in the instance so the plotting algorithm can access them
        # get the bounds from drawing a rectangle:
        if self.height_channel in self.channels:
            data = self.all_data[self.channels.index(self.height_channel)]
        else:
            data, trash = self._load_data([self.height_channel])
        if bounds is None:
            coords = select_rectangle(data, self.height_channel)
            if axis == 1:
                lower = coords[0][1]
                upper = coords[1][1]
            elif axis == 0:
                lower = coords[0][0]
                upper = coords[1][0]
        else:
            lower = bounds[0]
            upper = bounds[1]
        self._write_to_logfile('realign_axis_bounds', [axis, [lower, upper]])
        if self.height_channel in self.channels:
            height_data = self.all_data[self.channels.index(self.height_channel)]
        else:
            height_data_array, trash = self._load_data([self.height_channel])
            height_data = height_data_array[0]
            # if the channels have been scaled, the height has to be scaled as well
            scaling = self._get_channel_scaling(0)
            if scaling != 1:
                height_data = self._scale_array(height_data, self.height_channel, scaling)
        YRes = len(height_data)
        XRes = len(height_data[0])
        if axis == 1:
            reduced_height_data = np.zeros((upper-lower +1,XRes))
            for y in range(YRes):
                if (lower <= y) and (y <= upper):
                    for x in range(XRes):
                        reduced_height_data[y-lower][x] = height_data[y][x]
        elif axis == 0:
            reduced_height_data = np.zeros((YRes, upper-lower +1))
            for y in range(YRes):
                for x in range(XRes):
                    if (lower <= x) and (x <= upper):
                        reduced_height_data[y][x-lower] = height_data[y][x]
        shifts = self._get_mean_from_area(reduced_height_data, axis, threshold)

        # plot 
        fig, axs = plt.subplots()    
        fig.set_figheight(self.figsizey)
        fig.set_figwidth(self.figsizex) 
        cmap = SNOM_height
        img = axs.pcolormesh(height_data, cmap=cmap)
        # axs.invert_yaxis()
        divider = make_axes_locatable(axs)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('height (nm)', rotation=270)
        axs.set_title('Realigned')
        axs.axis('scaled')
        if axis == 1:
            axs.plot(range(XRes), [element + lower for element in shifts], color='red')
            axs.hlines([upper, lower], xmin=0, xmax=XRes, color='white')
        elif axis == 0:
            axs.plot([element + lower for element in shifts], range(YRes), color='red')
            axs.vlines([upper, lower], ymin=0, ymax=YRes, color='white')
        plt.show()

        # reinitialize the instance data to fit the new bigger arrays
        min_shift = round(min(shifts))
        max_shift = round(max(shifts))
        new_YRes = YRes + int(abs(min_shift-max_shift))
        all_data = self.all_data
        self.all_data = []
        for i in range(len(self.channels)):
            shifted_data = self._shift_data(all_data[i], axis, shifts)
            
            self.all_data.append(shifted_data)
            self.channels_label[i] += '_shifted'
            # adjust the scan area and pixel area
            if axis == 1:
                xres, yres, *args = self._get_channel_tag_dict_value(self.channels[i], ChannelTags.PIXELAREA)
                yres_new = new_YRes
                # new_values = [xres, yres_new, *args]
                self._set_channel_tag_dict_value(self.channels[i], ChannelTags.PIXELAREA, [xres, yres_new, *args])
                xreal, yreal, *args = self._get_channel_tag_dict_value(self.channels[i], ChannelTags.SCANAREA)
                yreal_new = yres_new*yreal/yres
                self._set_channel_tag_dict_value(self.channels[i], ChannelTags.SCANAREA, [xreal, yreal_new, *args])
            elif axis == 0:
                xres, yres, *args = self._get_channel_tag_dict_value(self.channels[i], ChannelTags.PIXELAREA)
                xres_new = new_YRes
                self._set_channel_tag_dict_value(self.channels[i], ChannelTags.PIXELAREA, [xres_new, yres, *args])
                xreal, yreal, *args = self._get_channel_tag_dict_value(self.channels[i], ChannelTags.SCANAREA)
                xreal_new = xres_new*xreal/xres
                self._set_channel_tag_dict_value(self.channels[i], ChannelTags.SCANAREA, [xreal_new, yreal, *args])
        gc.collect()

    def cut_channels(self, channels:list=None, preview_channel:str=None, autocut:bool=False, coords:list=None, reset_mask:bool=True) -> None:
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
                print('There does not seem to be an old mask... ')
        # generate new mask by selecting a region in the preview channel
        elif autocut is False:
            data = self.all_data[self.channels.index(preview_channel)]
            # get the coordinates of the selection rectangle
            if coords is None:
                coords = select_rectangle(data, preview_channel)
            # check if coords are none, if so, the user has canceled the selection
            if coords is not None:
                self._write_to_logfile('cut_coords', coords)
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
                    # set all values outside of the mask to zero and then cut all zero away from the outside with _auto_cut_channels(channels)
                    self.all_data[index] = np.multiply(self.all_data[index], self.mask_array)
        # apply the auto cut function to remove masked areas around the data
        self._auto_cut_channels(channels)
        gc.collect()

    def _auto_cut_channels(self, channels:list=None) -> None:
        """This function automatically cuts away all rows and lines which are only filled with zeros.
        This function applies to all channels in memory.

        Args:
            channels (list, optional): List of channels to apply the cut to. If not specified all channels in memory will be used. Defaults to None.
        """
        if channels is None:
            channels = self.channels
        
        # get the new size of the reduced channels
        reduced_data = self._auto_cut_data(self.all_data[0])
        yres = len(reduced_data)
        xres = len(reduced_data[0])
        for channel in channels:
            index = self.channels.index(channel)
            # get the old size of the data
            xres, yres, *args = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELAREA)
            xreal, yreal, *args = self._get_channel_tag_dict_value(channel, ChannelTags.SCANAREA)
            self.all_data[index] = self._auto_cut_data(self.all_data[index])
            xres_new = len(self.all_data[index][0])
            yres_new = len(self.all_data[index])
            xreal_new = xreal*xres_new/xres
            yreal_new = yreal*yres_new/yres
            # save new resolution and scan area in channel tag dict:
            self._set_channel_tag_dict_value(channel, ChannelTags.PIXELAREA, [xres_new, yres_new])
            self._set_channel_tag_dict_value(channel, ChannelTags.SCANAREA, [xreal_new, yreal_new])
            # add new appendix to channel
            self.channels_label[index] += '_reduced'
        self._write_to_logfile('cut', 'autocut')

    def _auto_cut_data(self, data) -> np.array:
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

    def scalebar(self, channels:list=[], units="m", dimension="si-length", label=None, length_fraction=None, height_fraction=None, width_fraction=None,
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
            XRes, YRes, *args = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELAREA)
            XReal, YReal, *args = self._get_channel_tag_dict_value(channel, ChannelTags.SCANAREA)
            pixel_scaling = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELSCALING)
            dx = XReal/(XRes)
            scalebar_var = [dx, units, dimension, label, length_fraction, height_fraction, width_fraction,
                            location, loc, pad, border_pad, sep, frameon, color, box_color, box_alpha, scale_loc,
                            label_loc, font_properties, label_formatter, scale_formatter, fixed_value, fixed_units, animated, rotation]
            if (channel in channels) or (len(channels)==0):
                self.scalebar_channels.append([channel, scalebar_var])                
            else:
                self.scalebar_channels.append([channel, None])                
            count += 1

    def rotate_90_deg(self, orientation:str = 'right'):
        """This function will rotate all data in memory by 90 degrees.

        Args:
            orientation (str, optional): rotate clockwise ('right') or counter clockwise ('left'). Defaults to 'right'.
        """
        if orientation == 'right':
            axes=(1,0)
            self._write_to_logfile('rotation', +90)
        elif orientation == 'left':
            axes=(0,1)
            self._write_to_logfile('rotation', -90)
        #rotate data:
        all_data = self.all_data.copy()
        # initialize data array
        self.all_data = []
        for channel in self.channels:
            # flip pixelarea and scanarea as well
            XReal, YReal, *args = self._get_channel_tag_dict_value(channel, ChannelTags.SCANAREA)
            self._set_channel_tag_dict_value(channel, ChannelTags.SCANAREA, [YReal, XReal])
            XRes, YRes, *args = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELAREA)
            self._set_channel_tag_dict_value(channel, ChannelTags.PIXELAREA, [YRes, XRes])
            self.all_data.append(np.rot90(all_data[self.channels.index(channel)], axes=axes))

    def _get_positions_from_plot(self, channel, data, coordinates:list=None, orientation=None) -> list:
        # Todo redundant to the get clicker corrdinates function?!
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
            self._plot_profile_lines(data, ax, coordinates, orientation)
        plt.title('Please select one or more points to continue.')
        plt.tight_layout()
        plt.show()
        klicker_coords = klicker.get_positions()['event'] #klicker returns a dictionary for the events
        klick_coordinates = [[round(element[0]), round(element[1])] for element in klicker_coords]
        return klick_coordinates

    def _get_profile(self, data, coordinates:list, orientation:Definitions, width:int) -> list:
        YRes = len(data)
        XRes = len(data[0])
        all_profiles = []
        for coord in coordinates:
            profile = []
            if orientation == Definitions.vertical:
                for y in range(YRes):
                    value = 0
                    for x in range(int(coord[0] - width/2), int(coord[0] + width/2)):
                        value += data[y][x]
                    value = value/width
                    profile.append(value)
            if orientation == Definitions.horizontal:
                for x in range(XRes):
                    value = 0
                    for y in range(int(coord[1] - width/2), int(coord[1] + width/2)):
                        value += data[y][x]
                    value = value/width
                    profile.append(value)
            all_profiles.append(profile)
        return all_profiles

    def select_profile(self, profile_channel:str, preview_channel:str=None, orientation:Definitions=Definitions.vertical, width:int=10, phase_orientation:int=1, coordinates:list=None):
        # Todo
        """This function lets the user select a profile with given width in pixels and displays the data.
        This is quite unfinished and only allows for profiles which extend over the whole image in the x-direction or y-direction.

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
            coordinates = self._get_positions_from_plot(preview_channel, previewdata)

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

        profiles = self._get_profile(profiledata, coordinates, orientation, width)
        for profile in profiles:
            xvalues = np.linspace(0, 10, len(profile))
            plt.plot(xvalues, profile, 'x')
        plt.title('Phase profiles')
        plt.tight_layout()
        plt.show()

        flattened_profiles = [phase_analysis.flatten_phase_profile(profile, phase_orientation) for profile in profiles]
        for profile in flattened_profiles:
            xvalues = np.linspace(0, 10, len(profile))
            plt.plot(xvalues, profile)
        plt.title('Flattened phase profiles')
        plt.tight_layout()
        plt.show()

        difference_profile = phase_analysis.get_profile_difference(profiles[0], profiles[1])
        # difference_profile = get_profile_difference(flattened_profiles[0], flattened_profiles[1])
        xres, yres = self._get_channel_tag_dict_value(self.channels.index(profile_channel), ChannelTags.PIXELAREA)
        xreal, yreal = self._get_channel_tag_dict_value(self.channels.index(profile_channel), ChannelTags.SCANAREA)
        pixel_scaling = self._get_channel_tag_dict_value(self.channels.index(profile_channel), ChannelTags.PIXELSCALING)
        xvalues = [i*yreal/yres/pixel_scaling for i in range(yres*pixel_scaling)]
        plt.plot(xvalues, difference_profile)
        plt.xlabel('Y [µm]')
        plt.ylabel('Phase difference')
        plt.ylim(ymin=0, ymax=2*np.pi)
        plt.title('Phase difference')
        plt.tight_layout()
        plt.show()
        gc.collect()

    def _plot_data_and_profile_pos(self, channel, data, coordinates, orientation):
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
        self._plot_profile_lines(data, ax, coordinates, orientation)
        plt.title('You chose the following line profiles')
        plt.tight_layout()
        plt.show()

    def _plot_profile_lines(self, data, ax, coordinates, orientation):
        xcoord = [coord[0] for coord in coordinates]
        ycoord = [coord[1] for coord in coordinates]
        if orientation == Definitions.vertical:
            ax.vlines(xcoord, ymin=0, ymax=len(data))
        elif orientation == Definitions.horizontal:
            ax.hlines(ycoord, xmin=0, xmax=len(data[0]))

    def _get_profiles_Coordinates(self, profile_channel, profiledata, preview_channel, previewdata, orientation, redo:bool=False, coordinates=None, redo_coordinates=None):
        if redo == False:
            coordinates = self._get_positions_from_plot(preview_channel, previewdata)
        else:
            display_coordinates = [coordinates[i] for i in range(len(coordinates)) if i not in redo_coordinates]# remove coordinates to redo and plot the other ones while selecton is active
            redone_coordinates = self._get_positions_from_plot(preview_channel, previewdata, display_coordinates, orientation)
            count = 0
            for index in redo_coordinates:
                coordinates[index] = redone_coordinates[count]
                count += 1

        self._plot_data_and_profile_pos(profile_channel, profiledata, coordinates, orientation)
        print('Are you satisfied with the profile positions? Or would you like to change one ore more profile positions?')
        user_input_bool = self._user_input_bool() 
        if user_input_bool == False:
            user_input = self._user_input('Please enter the indices of the profiles you like to redo, separated by a space character e.g. (0 1 3 11 ...)\nYour indices: ') 
            redo_coordinates = user_input.split(' ')
            redo_coordinates = [int(coord) for coord in redo_coordinates]
            print('coordinates to redo: ', redo_coordinates)
            print('Please select the new positons only for the indices you selected and in the same ordering, those were: ', redo_coordinates)
            coordinates = self._get_profiles_Coordinates(profile_channel, profiledata, preview_channel, previewdata, orientation, redo=True, coordinates=coordinates, redo_coordinates=redo_coordinates)
        
        return coordinates

    def select_profiles(self, profile_channel:str, preview_channel:str=None, orientation:Definitions=Definitions.vertical, width:int=10, coordinates:list=None):
        # Todo
        """This function lets the user select multiple profiles with given width in pixels and displays the data.
        Also unfinished, but allows for the selection of multiple profiles.

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
            self._initialize_data([profile_channel, preview_channel])#this will negate any modifications done prior like blurr...
        profiledata = self.all_data[self.channels.index(profile_channel)]
        previewdata = self.all_data[self.channels.index(preview_channel)]

        if coordinates == None:
            coordinates = self._get_profiles_Coordinates(profile_channel, profiledata, preview_channel, previewdata, orientation)
        
        print('The final profiles are shown in this plot.')
        self._plot_data_and_profile_pos(profile_channel, profiledata, coordinates, orientation)
        # get the profile data and save to class variables
        # additional infos are also stored and can be used by plotting and analysis functions
        self.profiles = self._get_profile(profiledata, coordinates, orientation, width)
        self.profile_channel = profile_channel
        self.profile_orientation = orientation
        return self.profiles

    def select_profiles_SSH(self, profile_channel_amp:str, profile_channel_phase:str, preview_channel:str=None, orientation:Definitions=Definitions.vertical, width_amp:int=10, width_phase:int=1, coordinates:list=None):
        # Todo
        """This function lets the user select a profile with given width in pixels and displays the data.
        Specific function for ssh model measurements. This will create a plot of field per waveguide index for the topological array.
        The field is calculated from the amplitude profiles times the cosine of the phasedifference to the central waveguide. 
        Also a very specific function for me, will probably not make it into the final version of the software.

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
            self._initialize_data([profile_channel_amp, profile_channel_phase, preview_channel])#this will negate any modifications done prior like blurr...
        profiledata_amp = self.all_data[self.channels.index(profile_channel_amp)]
        profiledata_phase = self.all_data[self.channels.index(profile_channel_phase)]
        previewdata = self.all_data[self.channels.index(preview_channel)]
        # get the profile coordinates
        if coordinates == None:
            coordinates = self._get_profiles_Coordinates(profile_channel_phase, profiledata_phase, preview_channel, previewdata, orientation)
        print(f'You selected the following coordinates: ', coordinates)
        print('The final profiles are shown in this plot.')
        self._plot_data_and_profile_pos(profile_channel_phase, profiledata_phase, coordinates, orientation)
        self._plot_data_and_profile_pos(profile_channel_amp, profiledata_amp, coordinates, orientation)
        self.profile_channel = profile_channel_phase
        self.profile_orientation = orientation

        # get the profile data for amp and phase
        self.phase_profiles = self._get_profile(profiledata_phase, coordinates, orientation, width_phase)
        # test:
        self._display_profile([self.phase_profiles[6], self.phase_profiles[16]])

        self.amp_profiles = self._get_profile(profiledata_amp, coordinates, orientation, width_amp)
        mean_amp = [np.mean(amp) for amp in self.amp_profiles]
        reference_index = int((len(self.phase_profiles)-1)/2)
        # phase_difference_profiles = [Phase_Analysis.get_profile_difference(self.phase_profiles[reference_index], self.phase_profiles[i]) for i in range(len(self.phase_profiles))]
        flattened_profiles = [phase_analysis.flatten_phase_profile(profile, +1) for profile in self.phase_profiles]
        self._display_profile(flattened_profiles, linestyle='-', title='Flattened phase profiles') # display the flattened profiles
        # phase_difference_profiles = [Phase_Analysis.get_profile_difference_2(self.phase_profiles[reference_index], self.phase_profiles[i]) for i in range(len(self.phase_profiles))]
        phase_difference_profiles = [phase_analysis.get_profile_difference_2(flattened_profiles[reference_index], flattened_profiles[i]) for i in range(len(flattened_profiles))]
        self._display_profile(phase_difference_profiles, linestyle='-', title='Phase difference to center wg') # display the phase difference profiles, no jumps close to 2 pi should occure or the average will lead to false values!
        # mean_phase_differences = [np.mean(diff) for diff in phase_difference_profiles]# todo this does not work!
        mean_phase_differences = [np.mean(diff) if np.mean(diff)>0 else np.mean(diff) + np.pi*2 for diff in phase_difference_profiles]# todo this does not work!
        real_per_wg_index = [mean_amp[i]*np.cos(mean_phase_differences[i]) for i in range(len(self.phase_profiles))]
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
        
    def _display_profile(self, profiles, ylabel=None, labels=None, linestyle='x', title=None):
        if self.profile_orientation == Definitions.horizontal:
            xrange, yrange = self._get_channel_tag_dict_value(self.channels.index(self.profile_channel), ChannelTags.SCANAREA)
            x_center_pos, y_center_pos = self._get_channel_tag_dict_value(self.channels.index(self.profile_channel), ChannelTags.SCANNERCENTERPOSITION)
            xres, yres = self._get_channel_tag_dict_value(self.channels.index(self.profile_channel), ChannelTags.PIXELAREA)
            xvalues = [x_center_pos - xrange/2 + x*(xrange/xres) for x in range(xres)]
            xlabel = 'X [µm]'
            if title == None:
                title = 'Horizontal profiles of channel ' + self.profile_channel
        elif self.profile_orientation == Definitions.vertical:
            xrange, yrange = self._get_channel_tag_dict_value(self.channels.index(self.profile_channel), ChannelTags.SCANAREA)
            x_center_pos, y_center_pos = self._get_channel_tag_dict_value(self.channels.index(self.profile_channel), ChannelTags.SCANNERCENTERPOSITION)
            xres, yres = self._get_channel_tag_dict_value(self.channels.index(self.profile_channel), ChannelTags.PIXELAREA)
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
        plt.legend()
        plt.tight_layout()
        plt.show()

    def display_profiles(self, ylabel:str=None, labels:list=None):
        """This function will display all current profiles from memory.

        Args:
            ylabel (str, optional): label of the y axis. The x axis label is in µm per default. Defaults to None.
            labels (list, optional): the description of the profiles. Will be displayed in the legend. Defaults to None.
        """
        self._display_profile(self.profiles)
        gc.collect()

    def display_flattened_profile(self, phase_orientation:int):
        """This function will flatten all profiles in memory and display them. Only useful for phase profiles!

        Args:
            phase_orientation (int): direction of the phase, must be '1' or '-1'
        """
        flattened_profiles = [phase_analysis.flatten_phase_profile(profile, phase_orientation) for profile in self.profiles]
        self._display_profile(flattened_profiles)
        gc.collect()

    def display_phase_difference(self, reference_index:int):
        """This function will calculate the phase difference of all profiles relative to the profile specified by the reference index.

        Args:
            reference_index (int): index of the reference profile. Basically the nth-1 selected profile.
        """
        difference_profiles = [phase_analysis.get_profile_difference(self.profiles[reference_index], self.profiles[i]) for i in range(len(self.profiles)) if i != reference_index]
        labels = ['Wg index ' + str(i) for i in range(len(difference_profiles))]
        self._display_profile(difference_profiles, 'Phase difference', labels)
        gc.collect()

    def _get_mean_phase_difference(self, profiles, reference_index:int):
        difference_profiles = [phase_analysis.get_profile_difference(profiles[reference_index], profiles[i]) for i in range(len(profiles)) if i != reference_index]
        mean_differences = [np.mean(diff) for diff in difference_profiles]
        return mean_differences

    def _scale_data_xy(self, data:np.array, scale_x:int, scale_y:int) -> np.array:
        XRes = len(data[0])
        YRes = len(data)
        new_data = np.zeros((YRes*scale_y, XRes*scale_x))
        for y in range(YRes):
            for i in range(scale_y):
                for x in range(XRes):
                    for j in range(scale_x):
                        new_data[y*scale_y + i][x*scale_x + j]= data[y][x]
        return new_data

    def quadratic_pixels(self, channels:list=None):
        """This function scales the data such that each pixel is quadratic, eg. the physical dimensions are equal.
        This is important because the pixels will be set to quadratic in the plotting function.
        However make shure that the pixel scaling x relative to y is an integer, otherwise the scaling will not work properly.
        This function will be applied to all channels in memory automatically when creating a measurement instance if autoscale is set to True.
        
        Args:
            channels [list]: list of channels the scaling should be applied to. If not specified the scaling will be applied to all channels
        """
        self._write_to_logfile('quadratic_pixels', True)
        if channels == None:
            channels = self.channels
        for channel in channels:
            if channel in self.channels:
                XRes, YRes, *args = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELAREA)
                XReal, YReal, *args = self._get_channel_tag_dict_value(channel, ChannelTags.SCANAREA)
                pixel_size_x = round(XReal/XRes *1000000000) # pixel size in nm
                pixel_size_y = round(YReal/YRes *1000000000)
                scale_x = 1
                scale_y = 1
                # if pixel_size_x < pixel_size_y:
                #     scale_y = int(pixel_size_y/pixel_size_x)
                # elif pixel_size_x > pixel_size_y:
                #     scale_x = int(pixel_size_x/pixel_size_y)
                # if pixel_size_x/scale_x != pixel_size_y/scale_y:
                    # print('The pixel size does not fit perfectly, you probably chose weired resolution values. You should probably not use this function then...\nScaling the data anyways!')
                # self.all_data[self.channels.index(channel)] = self._scale_data_xy(self.all_data[self.channels.index(channel)], scale_x, scale_y)
                # self._set_channel_tag_dict_value(channel, ChannelTags.PIXELAREA, [XRes*scale_x, YRes*scale_y])
                ###### New method using pillow to scale the image with interpolation, better if the scaling is not an integer
                # one could also implement a method using pillow to scale the image with interpolation, better if the scaling is not an integer
                rescaling = False
                if pixel_size_x < pixel_size_y:
                    # scale_y, rest = divmod(pixel_size_y, pixel_size_x)
                    xres = XRes
                    yres = int(YRes*pixel_size_y/pixel_size_x)
                    rescaling = True
                elif pixel_size_x > pixel_size_y:
                    # scale_x, rest = divmod(pixel_size_x, pixel_size_y)
                    yres = YRes
                    xres = int(XRes*pixel_size_x/pixel_size_y)
                    rescaling = True
                if rescaling:
                    img = Image.fromarray(self.all_data[self.channels.index(channel)])
                    img = img.resize((xres, yres), Image.Resampling.NEAREST)
                    self.all_data[self.channels.index(channel)] = np.array(img)
                    self._set_channel_tag_dict_value(channel, ChannelTags.PIXELAREA, [xres, yres])

    def overlay_forward_and_backward_channels(self, height_channel_forward:str, height_channel_backward:str, channels:list=None):
        """This function is ment to overlay the backwards and forwards version of the specified channels.
        The function will create a mean version which can then be displayed and saved. Note that the new version will be larger then the previous ones.
        Also make shure to use leveled data if you want to apply to height data.
        The overlain data will be larger, because the programm automatically tries to shift the data to match the best.
        The data will also be gauss blurred for better overlap.
        This function is still quite experimental and might not work properly in all cases. 
        But if it does you can basically double the integration time of your measurement.
        Works best for amplitude data, height data is also ok, to monitor the quality of the overlay.
        Phase channels don't work well, because the phase is not continuous and the mean phase is not meaningful because there are typically some slight shifts
        between forward and backward channel.

        Args:
            height_channel_forward (str): usual corrected height channel
            height_channel_backward (str): backwards height channel
            channels (list, optional): a list of all channels to be overlain. Defaults to None.
        """
        all_channels = []
        for channel in channels:
            all_channels.extend([channel, self.backwards_indicator + channel])
        all_channels.extend([height_channel_forward, height_channel_backward])
        self._initialize_data(all_channels)

        self.set_min_to_zero([height_channel_forward, height_channel_backward])
        
        #scale and blurr channels for better overlap
        self.scale_channels()
        # self.gauss_filter_channels_complex()

        height_data_forward = self.all_data[self.channels.index(height_channel_forward)]
        height_data_backward = self.all_data[self.channels.index(height_channel_backward)]
        
        #gauss blurr the data used for the alignment, so it might be a litte more precise
        height_channel_forward_blurr = self._gauss_blurr_data(height_data_forward, 2)
        height_channel_backward_blurr = self._gauss_blurr_data(height_data_backward, 2)

        # array_1 = height_data_forward[0]
        # array_2 = height_data_backward[0]

        '''
        mean_deviation_array = realign.Calculate_Squared_Deviation(array_1, array_2)
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
        pixel_scaling = self._get_channel_tag_dict_value(self.channels[0], ChannelTags.PIXELSCALING)[0]
        N = 5*pixel_scaling #maximum iterations, scaled if pixelnumber was increased

        # realign.minimize_deviation_1d(array_1, array_2, n_tries=N)
        # realign.Minimize_Deviation_2D(height_data_forward, height_data_backward, n_tries=N)

        # get the index which minimized the deviation of the height channels
        # index = realign.Minimize_Deviation_2D(height_data_forward, height_data_backward, N, False)
        index = realign.minimize_deviation_2d(height_channel_forward_blurr, height_channel_backward_blurr, N, False)
        # self.all_data[self.channels.index(height_channel_forward)], self.all_data[self.channels.index(height_channel_backward)] = realign.Shift_Array_2D_by_Index(height_data_forward, height_data_backward, index)


        for channel in channels:
            if self.backwards_indicator not in channel:
                #test:
                if self.height_indicator in channel:
                    # get current res and size and add the additional res and size due to addition of zeros while shifting
                    XRes, YRes, *args = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELAREA)
                    XReal, YReal, *args = self._get_channel_tag_dict_value(channel, ChannelTags.SCANAREA)
                    XRes_new = XRes + abs(index)# absolute value? index can be negative, but resolution can only increase, same for real dimensions
                    XReal_new = XReal + XReal/XRes*abs(index)
                    
                    # create channel_dict for new mean data 
                    self.channel_tag_dict.append(self.channel_tag_dict[self.channels.index(channel)])

                    # also create data dict entry
                    self.channels_label.append(self.channels_label[self.channels.index(channel)] + '_overlain')

                    # add new channel to channels
                    self.channels.append(channel + '_overlain')

                    self._set_channel_tag_dict_value(channel + '_overlain', ChannelTags.PIXELAREA, [XRes_new, YRes])
                    self._set_channel_tag_dict_value(channel + '_overlain', ChannelTags.SCANAREA, [XReal_new, YReal])

                    #test realign (per scan) based on minimization of differences 
                    #not usable right now, drift compensation might lead to differently sized data
                    # self.all_data[self.channels.index(height_channel_forward)] = realign.Minimize_Drift(self.all_data[self.channels.index(height_channel_forward)], display=False)
                    # self.all_data[self.channels.index(height_channel_backward)] = realign.Minimize_Drift(self.all_data[self.channels.index(height_channel_backward)])

                    # shift the data of the forward and backwards channel to match
                    self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)] = realign.Shift_Array_2D_by_Index(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)], index)
        

                    # create mean data and append to all_data
                    self.all_data.append(realign.Create_Mean_Array(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)]))
                else:
                    # get current res and size and add the additional res and size due to addition of zeros while shifting
                    XRes, YRes, *args = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELAREA)
                    XReal, YReal, *args = self._get_channel_tag_dict_value(channel, ChannelTags.SCANAREA)
                    XRes_new = XRes + abs(index)# absolute value? index can be negative, but resolution can only increase, same for real dimensions
                    XReal_new = XReal + XReal/XRes*abs(index)
                    
                    # create channel_dict for new mean data 
                    self.channel_tag_dict.append(self.channel_tag_dict[self.channels.index(channel)])

                    # also create data dict entry
                    self.channels_label.append(self.channels_label[self.channels.index(channel)] + '_overlain')

                    # add new channel to channels
                    self.channels.append(channel + '_overlain')
                    
                    self._set_channel_tag_dict_value(channel + '_overlain', ChannelTags.PIXELAREA, [XRes_new, YRes])
                    self._set_channel_tag_dict_value(channel + '_overlain', ChannelTags.SCANAREA, [XReal_new, YReal])

                    #test realign (per scan) based on minimization of differences 
                    # self.all_data[self.channels.index(channel)] = realign.Minimize_Drift(self.all_data[self.channels.index(channel)], display=False)
                    # self.all_data[self.channels.index(self.backwards_indicator+ channel)] = realign.Minimize_Drift(self.all_data[self.channels.index(self.backwards_indicator+ channel)])

                    # shift the data of the forward and backwards channel to match
                    self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)] = realign.Shift_Array_2D_by_Index(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)], index)

                    # create mean data and append to all_data
                    self.all_data.append(realign.create_mean_array(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)]))

        gc.collect()

    def overlay_forward_and_backward_channels_v2(self, height_channel_forward:str, height_channel_backward:str, channels:list=None):
        """
        Caution! This variant is ment to keep the scan size identical!

        This function is ment to overlay the backwards and forwards version of the specified channels.
        You should only specify the forward version of the channels you want to overlay. The function will create a mean version
        which can then be displayed and saved.

        Args:
            height_channel_forward (str): Usual corrected height channel
            height_channel_backward (str): Backwards height channel
            channels (list, optional): List of all channels to be overlain. Only specify the forward direction. Defaults to None. If not specified only the amp channels and the height channel will be overlain.
        """
        if channels is None:
            channels = [channel for channel in self.amp_channels if self.backwards_indicator not in channel]
            channels.append(self.height_channel)
        all_channels = []
        for channel in channels:
            if self.backwards_indicator not in channel:
                all_channels.extend([channel, self.backwards_indicator + channel]) # this is not optimal, what if the indicator does not come first?
        if height_channel_forward not in channels:
            all_channels.extend([height_channel_forward, height_channel_backward])
        self.initialize_channels(all_channels)
        self.set_min_to_zero([height_channel_forward, height_channel_backward])
        
        #scale channels for more precise overlap
        self.scale_channels()
        height_data_forward = self.all_data[self.channels.index(height_channel_forward)]
        height_data_backward = self.all_data[self.channels.index(height_channel_backward)]
        
        #gauss blurr the data used for the alignment, so it might be a litte more precise
        height_channel_forward_blurr = self._gauss_blurr_data(height_data_forward, 2)
        height_channel_backward_blurr = self._gauss_blurr_data(height_data_backward, 2)

        # try to optimize by shifting second array and minimizing mean deviation
        pixel_scaling = self._get_channel_tag_dict_value(self.channels[0], ChannelTags.PIXELSCALING)[0]
        N = 5*pixel_scaling #maximum iterations, scaled if pixelnumber was increased

        # get the index which minimized the deviation of the height channels
        index = realign.minimize_deviation_2d(height_channel_forward_blurr, height_channel_backward_blurr, N, False)

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
                    self.all_data.append(realign.create_mean_array_v2(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)], index))
                else:
                    # create channel_dict for new mean data 
                    self.channel_tag_dict.append(self.channel_tag_dict[self.channels.index(channel)])

                    # also create data dict entry
                    self.channels_label.append(self.channels_label[self.channels.index(channel)] + '_overlain')

                    # add new channel to channels
                    self.channels.append(channel + '_overlain')
                    
                    # create mean data and append to all_data
                    self.all_data.append(realign.create_mean_array_v2(self.all_data[self.channels.index(channel)], self.all_data[self.channels.index(self.backwards_indicator+ channel)], index))
        gc.collect()

    def manually_create_complex_channel(self, amp_channel:str, phase_channel:str, complex_type:str=None) -> None:
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
        demodulation_amp = self._get_demodulation_num(amp_channel)
        demodulation_phase = self._get_demodulation_num(phase_channel)
        if demodulation_amp != demodulation_phase:
            print('The channels you specified are not from the same demodulation order!\nProceeding anyways...')
            savefile_demod = str(demodulation_amp + ':' + demodulation_phase)
        else:
            savefile_demod = str(demodulation_amp)
        # check if channels are in memory, if not load the data
        if amp_channel not in self.channels:
            amp_data, amp_dict = self._load_data(amp_channel)
        else:
            amp_data = self.all_data[self.channels.index(amp_channel)]
            amp_dict = self.channel_tag_dict[self.channels.index(amp_channel)]
        if phase_channel not in self.channels:
            phase_data, phase_dict = self._load_data(phase_channel)
        else:
            phase_data = self.all_data[self.channels.index(phase_channel)]
            phase_dict = self.channel_tag_dict[self.channels.index(phase_channel)]
        # check if size is identical:
        xres_amp, yres_amp = amp_dict[ChannelTags.PIXELAREA]
        xres_phase, yres_phase = phase_dict[ChannelTags.PIXELAREA]
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
        real_channel = f'O{savefile_demod}' + self.real_indicator
        imag_channel = f'O{savefile_demod}' + self.imag_indicator
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
        gc.collect()

    def create_gif_old(self, amp_channel:str, phase_channel:str, frames:int=20, fps:int=10) -> None:
        """Old gif creation method.

        Args:
            amp_channel (str): _description_
            phase_channel (str): _description_
            frames (int, optional): _description_. Defaults to 20.
            fps (int, optional): _description_. Defaults to 10.
        """
        # Todo
        framenumbers=frames
        Duration=1000/fps # in ms

        realcolorpalette=[]
        # old color palette
        for i in range(0,255):
            realcolorpalette.append(i)
            if (i<127): realcolorpalette.append(i)
            else: realcolorpalette.append(255-i)
            realcolorpalette.append(255-i)

        if self.amp_indicator not in amp_channel or self.phase_indicator not in phase_channel:
            print('The specified channels are not specified as needed!')
            exit()
        demodulation_amp = self._get_demodulation_num(amp_channel)
        demodulation_phase = self._get_demodulation_num(phase_channel)
        if demodulation_amp != demodulation_phase:
            print('The channels you specified are not from the same demodulation order!\nProceeding anyways...')
            savefile_demod = str(demodulation_amp + ':' + demodulation_phase)
        else:
            savefile_demod = str(demodulation_amp)
        # check if channels are in memory, if not load the data
        if amp_channel not in self.channels or phase_channel not in self.channels:
            print('The channels for amplitude or phase were not found in the memory, they will be loaded automatically.\nBe aware that all prior modifications will get deleted.')
            # reload all channels
            self._initialize_data([amp_channel, phase_channel])
        amp_data = self.all_data[self.channels.index(amp_channel)]
        amp_dict = self.channel_tag_dict[self.channels.index(amp_channel)]
        phase_data = self.all_data[self.channels.index(phase_channel)]
        phase_dict = self.channel_tag_dict[self.channels.index(phase_channel)]
        xres_amp, yres_amp = amp_dict[ChannelTags.PIXELAREA]
        xres_phase, yres_phase = phase_dict[ChannelTags.PIXELAREA]
        if xres_amp != xres_phase or yres_amp != yres_phase:
            print('The data of the specified channels has different resolution!')
            exit()
        XRes, YRes = xres_amp, yres_amp
        flattened_amp = amp_data.flatten()
        maxval = max(flattened_amp)

        frames=[]
        for i in range(0,framenumbers):
            phase=i*2*np.pi/framenumbers
            repixels=[]
            colorpixels=[]
            for j in range(0,YRes):
                for k in range(XRes):
                    repixval=amp_data[j][k]*np.cos(phase_data[j][k]-phase)/maxval
                    repixels.append(repixval+1)
            img = Image.new('L', (XRes,YRes))
            # img = Image.fromarray(repixels)
            img.putdata(repixels,256/2,0)
            img.putpalette(realcolorpalette)
            #img=img.rotate(angle)
            #img=img.crop([int(YRes*np.sin(absangle)),int(XRes*np.sin(absangle)),int(XRes-YRes*np.sin(absangle)),int(YRes-XRes*np.sin(absangle))])
            #img.putdata(colorpixels,256,0)
            frames.append(img)
        channel = 'O' + savefile_demod + 'R'
        # self.filename is actually a windows path element not a str filename, to get the string use: self.filename.name
        # print('savefile path: ', self.directory_name / Path(self.filename.name + f'{channel}_gif.gif'))
        frames[0].save(self.directory_name / Path(self.filename.name + f'{channel}_gif_old.gif'), format='GIF', append_images=frames[1:], save_all=True,duration=Duration, loop=0)
        self._display_gif(self.directory_name / Path(self.filename.name + f'{channel}_gif_old.gif'), fps=fps)

    def create_gif(self, amp_channel:str, phase_channel:str, frames:int=20, fps:int=10, dpi=100) -> Path:
        """This function will create a gif from the amplitude and phase channel you specify. The gif will show the animated realpart by repeatedly adding a phase shift.
        The gif will be saved in the same directory as the measurement file and displayed afterwards.

        Args:
            amp_channel (str): Amplitude channel.
            phase_channel (str): Phase channel.
            frames (int, optional): Number of frames the gif should have. Defaults to 20.
            fps (int, optional): Frames per second. Defaults to 10.
            dpi (int, optional): Dots per inch. Defaults to 100.
        
        Returns:
            Path: Path to the saved gif.
        """
        framenumbers=frames
        Duration=1000/fps # in ms

        realcolorpalette=[]
        # old color palette
        for i in range(0,255):
            realcolorpalette.append(i)
            if (i<127): realcolorpalette.append(i)
            else: realcolorpalette.append(255-i)
            realcolorpalette.append(255-i)
        # convert cmap to colorpalette
        # realcolorpalette = SNOM_realpart
        # import matplotlib as mpl
        # norm = mpl.colors.Normalize()
        # from matplotlib import cm

        if self.amp_indicator not in amp_channel or self.phase_indicator not in phase_channel:
            print('The specified channels are not specified as needed!')
            exit()
        demodulation_amp = self._get_demodulation_num(amp_channel)
        demodulation_phase = self._get_demodulation_num(phase_channel)
        if demodulation_amp != demodulation_phase:
            print('The channels you specified are not from the same demodulation order!\nProceeding anyways...')
            savefile_demod = str(demodulation_amp + ':' + demodulation_phase)
        else:
            savefile_demod = str(demodulation_amp)
        # check if channels are in memory, if not load the data
        if amp_channel not in self.channels or phase_channel not in self.channels:
            print('The channels for amplitude or phase were not found in the memory, they will be loaded automatically.\nBe aware that all prior modifications will get deleted.')
            # reload all channels
            self._initialize_data([amp_channel, phase_channel])
        amp_data = self.all_data[self.channels.index(amp_channel)]
        amp_dict = self.channel_tag_dict[self.channels.index(amp_channel)]
        phase_data = self.all_data[self.channels.index(phase_channel)]
        phase_dict = self.channel_tag_dict[self.channels.index(phase_channel)]
        xres_amp, yres_amp = amp_dict[ChannelTags.PIXELAREA]
        xres_phase, yres_phase = phase_dict[ChannelTags.PIXELAREA]
        if xres_amp != xres_phase or yres_amp != yres_phase:
            print('The data of the specified channels has different resolution!')
            exit()
        XRes, YRes = xres_amp, yres_amp
        flattened_amp = amp_data.flatten()
        maxval = max(flattened_amp)

        frames=[]
        for i in range(0,framenumbers):
            phase=i*2*np.pi/framenumbers
            repixels=[]
            for j in range(0,YRes):
                for k in range(XRes):
                    repixval=amp_data[j][k]*np.cos(phase_data[j][k]-phase)/maxval
                    repixels.append(repixval+0.5)
            data = np.array(repixels).reshape(YRes, XRes)
            img = Image.fromarray(SNOM_realpart(data, bytes=True))
            frames.append(img)
        channel = 'O' + savefile_demod + 'R'
        # self.filename is actually a windows path element not a str filename, to get the string use: self.filename.name
        # print('savefile path: ', self.directory_name / Path(self.filename.name + f'{channel}_gif.gif'))
        gif_path = self.directory_name / Path(self.filename.name + f'{channel}_gif.gif')
        frames[0].save(gif_path, format='GIF', append_images=frames[1:], save_all=True,duration=Duration, loop=0, dpi=dpi)
        # plt.show()
        # plt.close(fig)
        if PlotDefinitions.show_plot:
            self._display_gif(gif_path, fps=fps)
        return gif_path

    def _display_gif(self, gif_path, fps=10):
        # Load the gif
        frames = imageio.mimread(gif_path)

        # Create a figure and axis
        fig, ax = plt.subplots()

        # Create a function to update the frame
        def update_image(frame):
            ax.clear()
            ax.imshow(frames[frame])
            # dont show frame around the image
            ax.axis('off')

        # Hide the axes
        ax.axis('off')

        # Create the animation
        ani = FuncAnimation(fig, update_image, frames=len(frames), interval=1000/fps, repeat=True)

        # Display the animation
        plt.show()

    def create_gif_v2(self, amp_channel:str, phase_channel:str, frames:int=20, fps:int=10) -> None:
        # Todo i dont even remember which version is best^^
        """This function will create a gif from the amplitude and phase channel you specify. The gif will show the animated realpart by repeatedly adding a phase shift.
        The gif will be saved in the same directory as the measurement file and displayed afterwards.

        Args:
            amp_channel (str): Amplitude channel.
            phase_channel (str): Phase channel.
            frames (int, optional): Number of frames the gif should have. Defaults to 20.
            fps (int, optional): Frames per second. Defaults to 10.
        """
        frame_numer = frames

        if self.amp_indicator not in amp_channel or self.phase_indicator not in phase_channel:
            print('The specified channels are not specified as needed!')
            exit()
        demodulation_amp = self._get_demodulation_num(amp_channel)
        demodulation_phase = self._get_demodulation_num(phase_channel)
        if demodulation_amp != demodulation_phase:
            print('The channels you specified are not from the same demodulation order!\nProceeding anyways...')
            savefile_demod = str(demodulation_amp + ':' + demodulation_phase)
        else:
            savefile_demod = str(demodulation_amp)
        # check if channels are in memory, if not load the data
        if amp_channel not in self.channels or phase_channel not in self.channels:
            print('The channels for amplitude or phase were not found in the memory, they will be loaded automatically.\nBe aware that all prior modifications will get deleted.')
            # reload all channels
            self._initialize_data([amp_channel, phase_channel])
        amp_data = self.all_data[self.channels.index(amp_channel)]
        amp_dict = self.channel_tag_dict[self.channels.index(amp_channel)]
        phase_data = self.all_data[self.channels.index(phase_channel)]
        phase_dict = self.channel_tag_dict[self.channels.index(phase_channel)]
        xres_amp, yres_amp = amp_dict[ChannelTags.PIXELAREA]
        xres_phase, yres_phase = phase_dict[ChannelTags.PIXELAREA]
        if xres_amp != xres_phase or yres_amp != yres_phase:
            print('The data of the specified channels has different resolution!')
            exit()
        XRes, YRes = xres_amp, yres_amp
        flattened_amp = amp_data.flatten()
        maxval = max(flattened_amp)
        cmap = SNOM_realpart

        # create real data for all frames
        self.all_real_data = []
        for i in range(0, frame_numer):
            phase = i*2*np.pi/frame_numer
            real_data = np.zeros((YRes, XRes))
            for j in range(0, YRes):
                for k in range(XRes):
                    real_data[j][k] = amp_data[j][k]*np.cos(phase_data[j][k]-phase)/maxval
            self.all_real_data.append(real_data)

        # Create figure and axis
        # figsize = 10
        # figsizex = 10
        # figsizey = 10*YRes/XRes
        fig, ax = plt.subplots(tight_layout=True) #, figsize=(figsizex, figsizey)
        
        # Create empty list to store the frames
        frames = []
        # Create the frames
        for i in range(frame_numer):
            ax.clear()
            data = self.all_real_data[i]
            self.cax = ax.pcolormesh(data, cmap=cmap, vmin=-maxval*1.1, vmax=maxval*1.1)
            # self.cax = ax.imshow(data, cmap=cmap, aspect='equal', vmin=-maxval*1.1, vmax=maxval*1.1)
            ax.set_aspect('equal')
            ax.invert_yaxis()
            ax.set_title('Frame {}'.format(i))
            if i == 0: # create colorbar only once
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size=f"{2}%", pad=0.05)
            cbar = plt.colorbar(self.cax, cax=cax)
            cbar.ax.get_yaxis().labelpad = 15
            cbar.ax.set_ylabel('Ez [arb.u.]', rotation=270)
            # remove ticks on x and y axis, they only show pixelnumber anyways, better to add a scalebar
            ax.set_xticks([])
            ax.set_yticks([])
            # disable the black frame around the image
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            # remove the whitespace around the image
            # ax.margins(0)
            # ax.margins(x=0, y=0)
            # ax.spines[['right', 'top']].set_visible(False)
            # disable the black frame around the colorbar
            cbar.outline.set_visible(False)
            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            frames.append(image)


        channel = 'O' + savefile_demod + 'R'
        # Save the frames as a gif
        imageio.mimsave(self.directory_name / Path(self.filename.name + f'{channel}_gif_v2.gif'), frames, fps=fps)
        # alternative:
        # import imageio.v3 as iio
        # iio.imwrite(self.directory_name / Path(self.filename.name + f'{channel}_gif_withimwrite.gif'), frames, fps=fps)
        # try with writer:
        # writer = imageio.get_writer(self.directory_name / Path(self.filename.name + f'{channel}_gif_with_writer.gif'), fps = fps)

        # for im in frames:
        #     writer.append_data(im)
        # writer.close()

        # delete the figure
        plt.close(fig)
        # display the gif
        self._display_gif(self.directory_name / Path(self.filename.name + f'{channel}_gif_v2.gif'), fps=fps)

    def substract_channels(self, channel1:str, channel2:str) -> None:
        """This function will substract the data of channel2 from channel1 and save the result in a new channel.
        The new channel will be named channel1-channel2.

        Args:
            channel1 (str): Channel from which the data will be substracted.
            channel2 (str): Channel which will be substracted from channel1.
        """
        if channel1 not in self.channels or channel2 not in self.channels:
            print('The specified channels are not in memory, they will be loaded automatically.')
            self._initialize_data([channel1, channel2])
        data1 = self.all_data[self.channels.index(channel1)]
        data2 = self.all_data[self.channels.index(channel2)]
        if data1.shape != data2.shape:
            print('The data of the specified channels has different resolution!')
            exit()
        result = data1 - data2
        self.channels.append(channel1 + '-' + channel2)
        self.all_data.append(result)
        self.channel_tag_dict.append(self.channel_tag_dict[self.channels.index(channel1)])
        self.channels_label.append(channel1 + '-' + channel2)

    def _select_data_range(self, channel:str, data:np.array=None, use_memory=True) -> tuple:
        """This function will use the data range selector to select a range of data. If use_memory is True the function will use the data from memory for the specified channel.
        In that case it will ignore the data argument. If use_memory is False the function will use the data argument and ignore the channel argument. The channel argument is only
        used to get the correct colormap. The function will return the selected data.
        Either one or two arrays will be returned depending on the selection.

        Args:
            data (np.array): Data array to select the range from. Defaults to None.
            channel (str): Channel name to get the data from memory or/and colormap from. Defaults to None.
            use_memory (bool, optional): If True the function will use the data from memory for the specified channel. Defaults to True.

        Returns:
            list: List of one or two arrays containing the selected data depending on the selection.
        """
        # identify the data to use for the range selection
        if use_memory:
            data = self.all_data[self.channels.index(channel)]
        elif data is None:
            print('No data was specified!')
            return None
        # get the range selection
        start, end, is_horizontal, inverted = select_data_range(data, channel)
        return start, end, is_horizontal, inverted

    def _get_data_from_selected_range(self, data:np.array, start:int, end:int, is_horizontal:bool, inverted:bool) -> list:
        """This function will return one or two arrays from the data using the coordinates of the range selection.

        Args:
            data (np.ndarray): Data array to create the array/s from.
            start (int): Start coordinate of the range selection.
            end (int): End coordinate of the range selection.
            is_horizontal (bool): Boolean to indicate if the range selection is horizontal.
            inverted (bool): Bollean to indicate if the range selection is inverted.

        Returns:
            list: The list contains one or two arrays depending on the selection. Each array contains the selected data.
        """
        # start, end, is_horizontal, inverted = self._select_data_range(channel, data, use_memory)
        # create one or two arrays from the data using the coordinates
        # print(f'start: <{start}>, end: <{end}>, is_horizontal: <{is_horizontal}>, inverted: <{inverted}>')
        # print(f'start type: <{type(start)}>, end type: <{type(end)}>, is_horizontal type: <{type(is_horizontal)}>, inverted type: <{type(inverted)}>')
        # print(f'data shape: {data.shape}')
        # print(f'data type: {type(data)}')
        reduced_data = []
        if is_horizontal:
            if inverted:
                left_data = data[:,:start]
                right_data = data[:,end:]
                reduced_data.append(left_data)
                reduced_data.append(right_data)
            else:
                selected_data = data[:,start:end]
                reduced_data.append(selected_data)
        else:
            if inverted:
                top_data = data[:start,:]
                bottom_data = data[end:,:]
                reduced_data.append(top_data)
                reduced_data.append(bottom_data)
            else:
                selected_data = data[start:end,:]
        return reduced_data
    
    def level_data_columnwise(self, channel_list:list=None, display_channel:str=None, selection:list=None) -> None:
        """This function will level the data of the specified channels columnwise.
        The function will use the data from the display channel to select the range for leveling.
        If no channels are specified all channels in memory will be leveled. If no display channel is specified the first channel in memory will be used.

        Args:
            channels (list, optional): Channels from memory which should be leveled. Defaults to None.
            display_channel (str, optional): Channel to use to select the range for leveling. Defaults to None.
            selection (list, optional): Selection to use for leveling. Defaults to None.
            You can specify the selection manually as a list with 4 elements like: [bound1 (int), bound2 (int), is_horizontal (bool), inverted (bool)].
        """
        # todo sofar only for the horizontal selection (slow drifts), maybe problematic if the data was rotated...
        # todo does not work yet for phase and amplitude channels
        # almost works but for phase channels phase jumps are an issue
        if channel_list is None:
            print('No channels specified, using all channels in memory.')
            channel_list = self.channels.copy() # make sure to use a copy for the iteration, because the list will be modified
        if display_channel is None:
            # preferably use a height channel:
            for channel in self.channels:
                if self.height_indicator in channel:
                    display_channel = channel
                    break
            if display_channel is None:
                display_channel = self.channels[0]
        # get the selection from the display channel
        if selection is None:
            selection = self._select_data_range(display_channel)
        # now use the selection to level all channels
        for channel in channel_list:
            # get the data from memory
            data = self.all_data[self.channels.index(channel)]
            # get the reduced data
            reduced_data = self._get_data_from_selected_range(data, *selection)
            # level the data
            if len(reduced_data) == 1:
                # print('leveling with one reference area')
                # get the reference data from the mean of the reduced data for each row
                reference_data = np.mean(reduced_data[0], axis=1)
                # create the leveled data
                leveled_data = np.zeros(data.shape)
                for i in range(data.shape[0]):
                    # leveled_data[i] = data[i] - reference_data[i]
                    if i > 0:
                        mean_drift = np.mean(reference_data[i]) - np.mean(reference_data[0])
                        leveled_data[i] = data[i] - mean_drift
                    else:
                        leveled_data[i] = data[i]
            elif len(reduced_data) == 2:
                # print('leveling with two reference areas')
                # get the reference data from the mean of the reduced data for each column and for both sides
                reference_data_left = np.mean(reduced_data[0], axis=1)
                reference_data_right = np.mean(reduced_data[1], axis=1)
                # create the leveled data by interpolating between the two reference data arrays and subtracting them from the data
                leveled_data = np.zeros(data.shape)
                for i in range(data.shape[0]):
                    # if phase is leveled make sure no phase jumps occur otherwise the leveling will not work
                    # first correct the overall drift of the mean per line
                    if i > 0:
                        mean_drift = np.mean([reference_data_left[i], reference_data_right[i]]) - np.mean([reference_data_left[0], reference_data_right[0]])
                        leveled_data[i] = data[i] - mean_drift
                    else:
                        leveled_data[i] = data[i]
                    # then correct the drift within each individual line by interpolating between the two reference data arrays
                    line_drift = np.interp(np.linspace(0, 1, data.shape[1]), [0, 1], [reference_data_left[i], reference_data_right[i]])
                    # shift line_drift such that the mean is zero
                    line_drift = line_drift - np.mean(line_drift)
                    leveled_data[i] = leveled_data[i] - line_drift
            # if phase channel, shift the data to match the leveled data to the original data
            if self.phase_indicator in channel:
                # todo, for now just shift by 0 to make sure the data is within the 0 to 2pi range
                # shift the data such that the mean is pi
                mean_phase = np.mean(leveled_data)
                shift = np.pi - mean_phase
                self._shift_phase_data(leveled_data, shift=shift)
            '''# save the leveled data, add the leveled data to memory and keep old data
            self.channels.append(channel + '_leveled')
            self.all_data.append(leveled_data)
            self.channel_tag_dict.append(self.channel_tag_dict[self.channels.index(channel)])
            self.channels_label.append(channel + '_leveled')'''
            # save the leveled data and replace old data
            # keep original channel name, but change the data and the channels_label
            self.all_data[self.channels.index(channel)] = leveled_data
            self.channels_label[self.channels.index(channel)] = channel + '_leveled'
        self._write_to_logfile('level_data_columnwise_selection', [channel_list, [elem for elem in selection]])

    def create_new_channel(self, data, channel_name:str, channel_tag_dict:dict, channel_label:str=None) -> None:
        """This function will create a new channel from the specified data and add it to memory.

        Args:
            data (np.array): Data array to create the new channel from.
            channel_name (str): Name of the new channel.
            channel_tag_dict (dict): Channel tag dictionary for the new channel.
            channel_label (str, optional): Label for the new channel. Defaults to None.
        """
        if channel_label is None:
            channel_label = channel_name
        self.channels.append(channel_name)
        self.all_data.append(data)
        self.channel_tag_dict.append(channel_tag_dict)
        self.channels_label.append(channel_label)

    # not yet fully implemented, eg. the profile plot function is only ment for full horizontal or vertical profiles only
    def test_profile_selection(self, channel:str=None) -> None:
        """Select a profile from the data. Allows the user to arbitrarily select a profile from the data.

        Args:
            channel (str, optional): channel to get the profile data from. Defaults to None.

        Returns:
            np.array, int, int, int: profile, start, end, width
        """
        if channel is None:
            channel = self.channels[0]
        
        array_2d = self.all_data[self.channels.index(channel)]
        # x, y = np.mgrid[-0:100:1, 0:200:1]
        # z = np.sqrt(x**2 + y**2) + np.sin(x**2 + y**2)
        # z = np.sin(x/2)*np.exp(-x/100)
        # array_2d = z
        # plt.pcolormesh(array_2d)
        # plt.show()
        profile, start, end, width = select_profile(array_2d, channel)
        # plt.plot(profile)
        # plt.show()
        return profile, start, end, width
        '''self.profile_channel = channel
        self.profiles = [profile]
        # find out the orientation of the profile
        if start[0] == end[0]:
            self.profile_orientation = Definitions.horizontal
        elif start[1] == end[1]:
            self.profile_orientation = Definitions.vertical
        else:
            self.profile_orientation = 'unknown'
            print('The profile orientation could not be determined!')'''

class ApproachCurve(FileHandler):
    """This class opens an approach curve measurement and handels all the approach curve related functions.
    
    Args:
        directory_name (str): Directory path of the measurement.
        channels (list, optional): List of channels to load. Defaults to None.
        title (str, optional): Title of the measurement. Defaults to None.
    """
    def __init__(self, directory_name:str, channels:list=None, title:str=None) -> None:
        self.measurement_type = MeasurementTypes.APPROACHCURVE
        if channels == None:
            channels = ['M1A']
        self.channels = channels.copy()
        self.x_channel = 'Z'
        super().__init__(directory_name, title)
        self.header = 27
        self._initialize_measurement_channel_indicators()
        self._load_data()
        # get the plotting style from the mpl style file
        self._load_mpl_style()

    def _initialize_measurement_channel_indicators(self):
        self.height_channel = 'Z'
        self.height_channels = ['Z']
        self.mechanical_channels = ['M1A', 'M1P']
        self.phase_channels = ['O1P','O2P','O3P','O4P','O5P']
        self.amp_channels = ['O1A','O2A','O3A','O4A','O5A']
        self.all_channels_default = self.height_channels + self.mechanical_channels + self.phase_channels + self.amp_channels
        self.height_indicator = self._get_from_config('height_indicator')
        self.amp_indicator = self._get_from_config('amp_indicator')
        self.phase_indicator = self._get_from_config('phase_indicator')
        self.backwards_indicator = self._get_from_config('backwards_indicator')
        self.real_indicator = self._get_from_config('real_indicator')
        self.imag_indicator = self._get_from_config('imag_indicator')
        self.optical_indicator = self._get_from_config('optical_indicator')
        self.mechanical_indicator = self._get_from_config('mechanical_indicator')
        self.channel_prefix_default = self._get_from_config('channel_prefix_default')
        self.channel_prefix_custom = self._get_from_config('channel_prefix_custom')
        self.channel_suffix_default = self._get_from_config('channel_suffix_default')
        self.channel_suffix_custom = self._get_from_config('channel_suffix_custom')
        self.channel_suffix_manipulated = self._get_from_config('channel_suffix_manipulated')
        self.channel_suffix_overlain = self._get_from_config('channel_suffix_overlain')
        self.file_ending = self._get_from_config('file_ending')
        self.phase_offset_default = self._get_from_config('phase_offset_default')
        self.phase_offset_custom = self._get_from_config('phase_offset_custom')
        self.rounding_decimal_amp_default = self._get_from_config('rounding_decimal_amp_default')
        self.rounding_decimal_amp_custom = self._get_from_config('rounding_decimal_amp_custom')
        self.rounding_decimal_phase_default = self._get_from_config('rounding_decimal_phase_default')
        self.rounding_decimal_phase_custom = self._get_from_config('rounding_decimal_phase_custom')
        self.rounding_decimal_complex_default = self._get_from_config('rounding_decimal_complex_default')
        self.rounding_decimal_complex_custom = self._get_from_config('rounding_decimal_complex_custom')
        self.rounding_decimal_height_default = self._get_from_config('rounding_decimal_height_default')
        self.rounding_decimal_height_custom = self._get_from_config('rounding_decimal_height_custom')
        self.height_scaling_default = self._get_from_config('height_scaling_default')
        self.height_scaling_custom = self._get_from_config('height_scaling_custom')

    def _load_data(self):
        self.all_data = {}
        datafile = self.directory_name / Path(self.filename.name + '.txt')
        x_channel_index = self.find_index(datafile, self.x_channel)
        with open(datafile, 'r') as file:
            xdata = np.genfromtxt(file ,skip_header=self.header, usecols=(x_channel_index), delimiter='\t', invalid_raise = False)
        self.all_data[self.x_channel] = xdata
        for channel in self.channels:
            channel_index = self.find_index(datafile, channel)
            with open(datafile, 'r') as file:
                y_data = np.genfromtxt(file ,skip_header=self.header, usecols=(channel_index), delimiter='\t', invalid_raise = False)
                self.all_data[channel] = y_data
        # scale the x data to nm
        x_scaling = 1
        x_unit = self._get_measurement_tag_dict_unit(MeasurementTags.SCANAREA)

        # we want to convert the xaxis to nm
        if x_unit == '[µm]':
            x_scaling = pow(10,3)
        elif x_unit == '[nm]':
            x_scaling = 1
        elif x_unit == '[m]':
                x_scaling = pow(10,9)
        # ok forget about that, the software from neaspec saves the scan area parameters as µm but the actual data is stored in m...
        x_scaling = pow(10,9)
        # scale xdata:
        self.all_data[self.x_channel] = np.multiply(self.all_data[self.x_channel], x_scaling)

    def set_min_to_zero(self) -> None:
        """This function will set the minimum of the xdata array to zero."""
        # set the min of the xdata array to zero
        min_x = np.nanmin(self.all_data[self.x_channel]) # for some reason at least the first value seems to be nan 
        self.all_data[self.x_channel] = self.all_data[self.x_channel] - min_x

    def display_channels(self, y_channels:list=None):
        """This function will display the specified channels in a plot. The x channel is always 'Z'.
        If no y channels are specified all channels in memory will be displayed.
        All channels will share one axis.

        Args:
            y_channels (list, optional): List of channels to display. Defaults to None.
        """
        # get the plotting style from the mpl style file
        self._load_mpl_style()
        if y_channels == None:
            y_channels = self.channels
        x_channel = 'Z'
        
        for channel in y_channels:
            plt.plot(self.all_data[self.x_channel], self.all_data[channel], label=channel)

        # labels for axes:
        plt.xlabel(f'Z [nm]')
        if len(self.channels) == 1:
            plt.ylabel(self.channels[0])
        plt.legend()
        if PlotDefinitions.tight_layout:
            plt.tight_layout()
        
        if PlotDefinitions.show_plot:
            plt.show()
    
    def display_channels_v2(self, y_channels:list=None):
        """This function will display the specified channels in a plot. The x channel is always 'Z'.
        If no y channels are specified all channels in memory will be displayed.
        Each channel will have its own axis. And if you specify two channels it will make use of the left and right axis.
        For more channels only the left axis will be used for the first channel.
        
        Args:
            y_channels (list, optional): List of channels to display. Defaults to None.
        """
        x_channel = 'Z'
        if y_channels == None:
            y_channels = self.channels
        y_data = []
        for channel in y_channels:
            y_data.append(self.all_data[channel])
        self._display_approach_curve(x_data=self.all_data[self.x_channel], y_data=y_data, x_channel=x_channel, y_channels=y_channels)

    def _display_approach_curve(self, x_data, y_data:list, x_channel, y_channels):        
        # import matplotlib.colors as mcolors
        # colors = mcolors.TABLEAU_COLORS
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:olive']
        fig, ax1 = plt.subplots()
        line1, = ax1.plot(x_data, y_data[0], label=y_channels[0], color=colors[0])
        if len(y_channels) == 1:
            ax1.legend()
        elif len(y_channels) == 2:
            ax2 = ax1.twinx()
            line2, = ax2.plot(x_data, y_data[1], label=y_channels[1], color=colors[1])
            ax2.set_ylabel(y_channels[1])
            ax1.legend(handles=[line1, line2])
        else: # deactivate ticks for all except the first or it will get messy
            handles = [line1]
            for channel in y_channels[1:]: # ignore the first as it was plotted already
                # i = self.channels.index(channel)
                i = y_channels.index(channel)
                ax = ax1.twinx()
                ax.tick_params(right=False, labelright=False)
                line, = ax.plot(x_data, y_data[i], label=channel, color=colors[i])
                handles.append(line)
            ax1.legend(handles=handles)

        # labels for axes:
        ax1.set_xlabel(f'Z [nm]')
        ax1.set_ylabel(y_channels[0])
        if PlotDefinitions.tight_layout:
            plt.tight_layout()
        if PlotDefinitions.show_plot:
            plt.show()
        gc.collect()


    def find_index(self, filepath, channel):
        """This function will find the index of the specified channel in the datafile.
        
        Args:
            filepath (str): Path to the datafile.
            channel (str): Channel to find the index for.
        """
        with open(filepath, 'r') as file:
            for i in range(self.header+1):
                line = file.readline()
        split_line = line.split('\t')
        split_line.remove('\n')
        return split_line.index(channel)

class Scan3D(FileHandler):
    """A 3D scan is a measurement where one approach curve is saved per pixel. This class is ment to handle such measurements.

    Args:
        directory_name (str): Directory path of the measurement.
        channels (list, optional): List of channels to load. Defaults to None.
        title (str, optional): Title of the measurement. Defaults to None.
    """
    def __init__(self, directory_name: str, channels:list=None, title: str = None) -> None:
        self.measurement_type = MeasurementTypes.SCAN3D
        # set channelname if none is given
        if channels == None:
            channels = ['Z', 'O2A', 'O2P'] # if you want to plot approach curves 'Z' must be included!
        self.channels = channels.copy()
        self.x_channel = 'Z'
        # call the init constructor of the filehandler class
        super().__init__(directory_name, title)
        # define header, probably same as for approach curve
        self.header = 27
        # initialize the channel indicators
        # print('filetype: ', self.file_type)
        self._initialize_measurement_channel_indicators()
        self._update_measurement_channel_indicators()
        # for some reason the naming convention does not always follow the default for the snom measurements of the same filetype
        try:
            self._create_channel_tag_dict()
        except:
            self.channel_suffix_default = ''
            try:
                self._create_channel_tag_dict()
            except:
                print('The channel tag dict could not be created!')
                exit()
        # load the channels from the datafile
        self._load_data()
        # get the plotting style from the mpl style file
        self._load_mpl_style()

    def _update_measurement_channel_indicators(self):
        self.height_channel = 'Z'
        self.height_channels = ['Z']
        self.mechanical_channels = ['M1A', 'M1P'] # todo
        self.phase_channels = ['O1P','O2P','O3P','O4P','O5P']
        self.amp_channels = ['O1A','O2A','O3A','O4A','O5A']
        self.all_channels_default = self.mechanical_channels + self.phase_channels + self.amp_channels
        self.all_channels_custom = self.height_channels
    
    def _load_data(self):
        datafile = self.directory_name / Path(self.filename.name + '.txt')
        # initialize all data dict
        self.all_data = {} # (key, value) = (channelname, 3d matrix, shape:(xres, yres, zres)) 
        # load the data per channel and add to all_data
        for channel in self.channels:
            index = find_index(self.header, datafile, channel) # find the index of the channels
            file = open(datafile, 'r')
            self.all_data[channel] = np.genfromtxt(file ,skip_header=self.header+1, usecols=(index), delimiter='\t', invalid_raise = False)
            file.close()
            x,y,z = self._get_measurement_tag_dict_value(MeasurementTags.PIXELAREA)
            self.all_data[channel] = np.reshape(self.all_data[channel], (y,x,z))
        # scale the x data to nm
        x_scaling = 1
        # try: x_unit = self.measurement_tag_dict[MeasurementTags.SCANAREA][0]
        x_unit = self._get_measurement_tag_dict_unit(MeasurementTags.SCANAREA)
        # except: x_unit = None
        # else:
        # we want to convert the xaxis to nm
        if x_unit == '[µm]':
            x_scaling = pow(10,3)
        elif x_unit == '[nm]':
            x_scaling = 1
        elif x_unit == '[m]':
                x_scaling = pow(10,9)
        # ok forget about that, the software from neaspec saves the scan area parameters as µm but the actual data is stored in m...
        x_scaling = pow(10,9)
        # scale xdata:
        self.all_data[self.x_channel] = np.multiply(self.all_data[self.x_channel], x_scaling)

    def set_min_to_zero(self) -> None:
        """This function will set the minimum of the xdata array to zero."""
        # set the min of the xdata array to zero
        min_x = np.nanmin(self.all_data[self.x_channel]) # for some reason at least the first value seems to be nan 
        self.all_data[self.x_channel] = self.all_data[self.x_channel] - min_x

    def get_cutplane_data(self, axis:str='x', line:int=0, channel:str=None) -> np.array:
        """This function will return the data of a cutplane of the 3D scan. The cutplane is defined by the axis and the line.
        
        Args:
            axis (str, optional): Axis of the cutplane. Defaults to 'x'.
            line (int, optional): Line of the cutplane. Defaults to 0.
            channel (str, optional): Channel to get the data from. Defaults to None.

        Returns:
            np.array: Data of the cutplane.
        """
        if channel == None:
            channel = self.channels[0]
        x,y,z = self._get_measurement_tag_dict_value(MeasurementTags.PIXELAREA)
        data = self.all_data[channel].copy()
        if axis == 'x':
            cutplane_data = np.zeros((z,x)) 
            for i in range(x):
                for j in range(z):
                    cutplane_data[j][i] = data[line][i][j]
        return cutplane_data

    def generate_all_cutplane_data(self, axis:str='x', line:int=0):
        """This function will generate the data of all cutplanes for all channels and store them in a dictionary.
        
        Args:
            axis (str, optional): Axis of the cutplane. Defaults to 'x'.
            line (int, optional): Line of the cutplane. Defaults to 0.
        """
        self.all_cutplane_data = {}
        for channel in self.channels:
            self.all_cutplane_data[channel] = self.get_cutplane_data(axis=axis, line=line, channel=channel)

    def _create_subplot(self, axis:str='x', line:int=0, channel:str=None, auto_align:bool=False):
        if channel == None:
            channel = self.channels[0]
        cutplane_data = self.all_cutplane_data[channel]
        # sadly the data definitions for this filytype are off, eg. missing 'raw' suffix for 3D scan, also the channel headers are incomplete, z res is false
        # XRes, YRes, ZRes = self._get_channel_tag_dict_value(channel, ChannelTags.PIXELAREA)
        # therefore we use the measurement tag dict
        XRes, YRes, ZRes = self._get_measurement_tag_dict_value(MeasurementTags.PIXELAREA)

        # YRes, XRes = cutplane_data.shape # cutplane data might have been
        XRange, YRange, ZRange = self._get_measurement_tag_dict_value(MeasurementTags.SCANAREA)
        XYZUnit = self._get_measurement_tag_dict_unit(MeasurementTags.SCANAREA)
        # convert Range to nm
        if XYZUnit == '[µm]':
            XRange = XRange*1e3
            YRange = YRange*1e3
            ZRange = ZRange*1e3
        else:
            print('Error! The unit of the scan area is not supported yet!')
        z_pixelsize = ZRange/ZRes

        # now we can try to shift each approach curve by the corresponding z_shift
        # easiest way is to use the z start position of each approach curve
        if auto_align:
            z_shifts = np.zeros(XRes)
            # idea: get all the lowest points of the approach curves and shift them to the same z position, herefore we shift them only upwards relative to the lowest point
            z_data = self.all_cutplane_data[self.x_channel]
            # reshape the data to the correct shape
            for i in range(XRes):
                z_shifts[i] = self._get_z_shift_(z_data[:,i])
            # z_data is in nm
            z_shifts = z_shifts
            z_min = np.min(z_shifts)
            z_shifts = z_shifts - z_min
            # therefore we need to create a new data array which can encorporate the shifted data
            # calculate the new z range
            ZRange_new = ZRange + z_shifts.max()
            ZRes_new = int(ZRange_new/z_pixelsize)
            # print('ZRes_new: ', ZRes_new)
            # create the new data array
            cutplane_data = np.zeros((ZRes_new, XRes))
            data = self.all_cutplane_data[channel].copy()
            for i in range(XRes):
                for j in range(ZRes):
                    cutplane_data[j+int(z_shifts[i]/z_pixelsize)][i] = data[j][i]
            # This shifting is not optimal, since a slow drift or a tilt of the sample would lead to a wrong alignment of the approach curves, although they start at the bottom.
            # Maybe try to use a 2d scan of the same region to align the approach curves.
        
        # import plotting_parameters.json, here the user can tweek some options for the plotting, like automatic titles and colormap choices
        plotting_parameters = self._get_plotting_parameters()

        # update the placeholders in the dictionary
        # the dictionary contains certain placeholders, which are now being replaced with the actual values
        # until now only the channel placeholder is used but more could be added
        # placeholders are indicated by the '<' and '>' characters
        # this step insures, that for example the title contains the correct channel name
        placeholders = {'<channel>': channel}
        plotting_parameters = self._replace_plotting_parameter_placeholders(plotting_parameters, placeholders)

        # set colormap depending on channel
        if self.amp_indicator in channel:
            cmap = plotting_parameters["amplitude_cmap"]
            label = plotting_parameters["amplitude_cbar_label"]
            title = plotting_parameters["amplitude_title"]
        elif self.phase_indicator in channel:
            cmap = plotting_parameters["phase_cmap"]
            label = plotting_parameters["phase_cbar_label"]
            title = plotting_parameters["phase_title"]
        elif self.height_indicator in channel:
            cmap = plotting_parameters["height_cmap"]
            label = plotting_parameters["height_cbar_label"]
            title = plotting_parameters["height_title"]
        else:
            cmap = 'viridis'
            label = 'unknown'
            title = 'unknown'
        return cutplane_data, cmap, label, title
    
    def display_cutplanes(self, axis:str='x', line:int=0, channels:list=None, auto_align:bool=False):
        """This function will display the cutplanes of the specified channels.
        You can also autoalign the data which will apply a shift to align the approach curves, more physically correct but not perfect.
        
        Args:
            axis (str, optional): Axis of the cutplane. Defaults to 'x'.
            line (int, optional): Line of the cutplane. Defaults to 0.
            channels (list, optional): Channel to display, if you don't specify some all channels in memory will be used. Defaults to None.
            align (bool, optional): Alignment of the approach curves.
            If set to True the individual approach curves will be shifted such that they start at the same Z corrdinate. Defaults to False.
        """
        # get the plotting style from the mpl style file
        self._load_mpl_style()
        if channels == None:
            channels = self.channels
        number_of_channels = len(channels)
        if number_of_channels == 1:
            cols = 1
        elif number_of_channels < 5:
            cols = 2
        else:
            cols = 3
        rows = number_of_channels//cols
        if number_of_channels%cols != 0:
            rows += 1
        # print('rows: ', rows)
        # print('cols: ', cols)
        fig, axs = plt.subplots(rows, cols, figsize=(PlotDefinitions.figsizex, PlotDefinitions.figsizey))
        for channel in channels:
            # get column and row index
            if number_of_channels < 5:
                col = channels.index(channel)%2
                row = channels.index(channel)//2
            else:
                col = channels.index(channel)%3
                row = channels.index(channel)//3
            if rows == 1 and cols == 1:
                ax = axs
            elif rows == 1:
                ax = axs[col]
            else:
                ax = axs[row, col]
            if channel not in self.channels:
                print(f'The channel <{channel}> is not in memory!')
                continue
            cutplane_data, cmap, label, title = self._create_subplot(axis=axis, line=line, channel=channel, auto_align=auto_align)
            img = ax.pcolormesh(cutplane_data, cmap=cmap, rasterized=True)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size=f"{calculate_colorbar_size(fig, ax, self.colorbar_width)}%", pad=0.05) # size is the size of colorbar relative to original axis, 100% means same size, 10% means 10% of original
            # f"{calculate_colorbar_size(fig, axis, self.colorbar_width)}%"
            cbar = plt.colorbar(img, aspect=1, cax=cax)
            cbar.ax.get_yaxis().labelpad = 15
            cbar.ax.set_ylabel(label, rotation=270)
            ax.axis('scaled')
            if self.hide_ticks == True:
                # remove ticks on x and y axis, they only show pixelnumber anyways, better to add a scalebar
                ax.set_xticks([])
                ax.set_yticks([])
            if self.show_titles == True:
                ax.set_title(title)
        #turn off all unneeded axes
        counter = 1
        for row in range(rows):
            for col in range(cols):
                if rows == 1 and cols ==1:
                    ax = axs
                elif rows == 1:
                    ax = axs[col]
                else:
                    ax = axs[row, col]
                if counter >= number_of_channels: 
                    ax.axis('off')
                counter += 1
        if self.tight_layout is True:
            plt.tight_layout()
        if PlotDefinitions.show_plot is True:
            plt.show()
        gc.collect()

    def display_cutplane_v2_realpart(self, axis:str='x', line:int=0, demodulation:int=2, align='auto'):
        """This function will display the cutplane of the realpart data of the channels of the specified demodulation order.
        The data will be shifted to align the approach curves.
        
        Args:
            axis (str, optional): Axis of the cutplane. Defaults to 'x'.
            line (int, optional): Line of the cutplane. Defaults to 0.
            demodulation (int, optional): Demodulation order of the data. Defaults to 2.
            align (str, optional): Alignment of the approach curves. Defaults to 'auto'.
        """
        amp_channel = f'O{demodulation}A'
        phase_channel = f'O{demodulation}P'
        x,y,z = self._get_measurement_tag_dict_value(MeasurementTags.PIXELAREA)
        amp_data = self.all_data[amp_channel].copy()
        phase_data = self.all_data[phase_channel].copy()
        if axis == 'x':
            cutplane_amp_data = np.zeros((z,x)) 
            cutplane_phase_data = np.zeros((z,x))
            for i in range(x):
                for j in range(z):
                    cutplane_amp_data[j][i] = amp_data[line][i][j]
                    cutplane_phase_data[j][i] = phase_data[line][i][j]
        # todo: shift each y column by offset value depending on average z position, to correct for varying starting position, due to non flat substrates
        z_shifts = np.zeros(x)
        # idea: get all the lowest points of the approach curves and shift them to the same z position, herefore we shift them only upwards relative to the lowest point
        z_data_raw = self.all_data[self.x_channel]
        # reshape the data to the correct shape
        if axis == 'x':
            z_data = np.zeros((z,x)) 
            for i in range(x):
                for j in range(z):
                    z_data[j][i] = z_data_raw[line][i][j]
        for i in range(x):
            z_shifts[i] = self._get_z_shift_(z_data[:,i])
        z_shifts = z_shifts
        if align == 'auto':
            z_min = np.min(z_shifts)
            z_shifts = z_shifts - z_min
        # now we need to shift each approach curve by the corresponding z_shift
        # therefore we need to create a new data array which can encorporate the shifted data
        XRes, YRes, ZRes = self._get_measurement_tag_dict_value(MeasurementTags.PIXELAREA)
        # print('ZR: ', ZRes)
        XRange, YRange, ZRange = self._get_measurement_tag_dict_value(MeasurementTags.SCANAREA)
        XYZUnit = self._get_measurement_tag_dict_unit(MeasurementTags.SCANAREA)
        # convert Range to nm
        if XYZUnit == '[µm]':
            XRange = XRange*1e3
            YRange = YRange*1e3
            ZRange = ZRange*1e3
        else:
            print('Error! The unit of the scan area is not supported yet!')
        z_pixelsize = ZRange/ZRes
        # print('z_shifts: ', z_shifts)
        # calculate the new z range
        ZRange_new = ZRange + z_shifts.max()
        ZRes_new = int(ZRange_new/z_pixelsize)
        # print('ZRes_new: ', ZRes_new)
        # create the new data array
        cutplane_real_data = np.zeros((ZRes_new, XRes))
        for i in range(XRes):
            for j in range(ZRes):
                cutplane_real_data[j+int(z_shifts[i]/z_pixelsize)][i] = amp_data[line][i][j]*np.cos(phase_data[line][i][j])
        # set the channel 
        channel = f'O{demodulation}Re'
        '''This shifting is not optimal, since a slow drift or a tilt of the sample would lead to a wrong alignment of the approach curves, although they start at the bottom.
        Maybe try to use a 2d scan of the same region to align the approach curves.'''
        
        # import plotting_parameters.json, here the user can tweek some options for the plotting, like automatic titles and colormap choices
        plotting_parameters = self._get_plotting_parameters()

        # update the placeholders in the dictionary
        # the dictionary contains certain placeholders, which are now being replaced with the actual values
        # until now only the channel placeholder is used but more could be added
        # placeholders are indicated by the '<' and '>' characters
        # this step insures, that for example the title contains the correct channel name
        placeholders = {'<channel>': channel}
        plotting_parameters = self._replace_plotting_parameter_placeholders(plotting_parameters, placeholders)

        # set colormap depending on channel
        if self.amp_indicator in channel:
            cmap = plotting_parameters["amplitude_cmap"]
            label = plotting_parameters["amplitude_cbar_label"]
            title = plotting_parameters["amplitude_title"]
        elif self.phase_indicator in channel:
            cmap = plotting_parameters["phase_cmap"]
            label = plotting_parameters["phase_cbar_label"]
            title = plotting_parameters["phase_title"]
        elif self.real_indicator in channel:
            cmap = plotting_parameters["real_cmap"]
            label = plotting_parameters["real_cbar_label"]
            title = plotting_parameters["real_title_real"]
        else:
            cmap = 'viridis'
        fig, ax = plt.subplots()
        max_val = np.max(cutplane_real_data)
        img = plt.pcolormesh(cutplane_real_data, cmap=cmap, vmin=-max_val, vmax=max_val)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=f"{self.colorbar_width}%", pad=0.05) # size is the size of colorbar relative to original axis, 100% means same size, 10% means 10% of original
        cbar = plt.colorbar(img, aspect=1, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(label, rotation=270)
        if self.hide_ticks == True:
            # remove ticks on x and y axis, they only show pixelnumber anyways, better to add a scalebar
            ax.set_xticks([])
            ax.set_yticks([])
        plt.tight_layout()
        # plt.colorbar(img)
        plt.show()
    
    def display_cutplane_realpart(self, axis:str='x', line:int=0, demodulation:int=2, align='auto'):
        """This function will display the cutplane of the realpart data of the channels of the specified demodulation order.
        The data will be shifted to align the approach curves.
        
        Args:
            axis (str, optional): Axis of the cutplane. Defaults to 'x'.
            line (int, optional): Line of the cutplane. Defaults to 0.
            demodulation (int, optional): Demodulation order of the data. Defaults to 2.
            align (str, optional): Alignment of the approach curves. Defaults to 'auto'.
        """
        amp_channel = f'O{demodulation}A'
        phase_channel = f'O{demodulation}P'
        real_channel = f'O{demodulation}Re'
        # set the channel 
        channel = f'O{demodulation}Re'
        if channel == None:
            channel = self.channels[0]
        # create real part cutplane data
        self.all_cutplane_data[real_channel] = np.multiply(self.all_cutplane_data[f'O{demodulation}A'], np.cos(self.all_cutplane_data[f'O{demodulation}P']))
        cutplane_data = self.all_cutplane_data[real_channel]
        XRes, YRes, ZRes = self._get_measurement_tag_dict_value(MeasurementTags.PIXELAREA)
        XRange, YRange, ZRange = self._get_measurement_tag_dict_value(MeasurementTags.SCANAREA)
        XYZUnit = self._get_measurement_tag_dict_unit(MeasurementTags.SCANAREA)
        # convert Range to nm
        if XYZUnit == '[µm]':
            XRange = XRange*1e3
            YRange = YRange*1e3
            ZRange = ZRange*1e3
        else:
            print('Error! The unit of the scan area is not supported yet!')
        z_pixelsize = ZRange/ZRes

        # now we can try to shift each approach curve by the corresponding z_shift
        # easiest way is to use the z start position of each approach curve
        if align == 'auto':
            z_shifts = np.zeros(XRes)
            # idea: get all the lowest points of the approach curves and shift them to the same z position, herefore we shift them only upwards relative to the lowest point
            z_data = self.all_cutplane_data[self.x_channel]
            # reshape the data to the correct shape
            for i in range(XRes):
                z_shifts[i] = self._get_z_shift_(z_data[:,i])
            # z_data is in nm
            z_shifts = z_shifts
            z_min = np.min(z_shifts)
            z_shifts = z_shifts - z_min
            # therefore we need to create a new data array which can encorporate the shifted data
            # calculate the new z range
            ZRange_new = ZRange + z_shifts.max()
            ZRes_new = int(ZRange_new/z_pixelsize)
            # print('ZRes_new: ', ZRes_new)
            # create the new data array
            cutplane_data = np.zeros((ZRes_new, XRes))
            data = self.all_cutplane_data[real_channel].copy()
            for i in range(XRes):
                for j in range(ZRes):
                    cutplane_data[j+int(z_shifts[i]/z_pixelsize)][i] = data[j][i]
            # This shifting is not optimal, since a slow drift or a tilt of the sample would lead to a wrong alignment of the approach curves, although they start at the bottom.
            # Maybe try to use a 2d scan of the same region to align the approach curves.
        
        '''This shifting is not optimal, since a slow drift or a tilt of the sample would lead to a wrong alignment of the approach curves, although they start at the bottom.
        Maybe try to use a 2d scan of the same region to align the approach curves.'''
        
        # import plotting_parameters.json, here the user can tweek some options for the plotting, like automatic titles and colormap choices
        plotting_parameters = self._get_plotting_parameters()

        # update the placeholders in the dictionary
        # the dictionary contains certain placeholders, which are now being replaced with the actual values
        # until now only the channel placeholder is used but more could be added
        # placeholders are indicated by the '<' and '>' characters
        # this step insures, that for example the title contains the correct channel name
        placeholders = {'<channel>': channel}
        plotting_parameters = self._replace_plotting_parameter_placeholders(plotting_parameters, placeholders)

        # set colormap depending on channel
        if self.amp_indicator in channel:
            cmap = plotting_parameters["amplitude_cmap"]
            label = plotting_parameters["amplitude_cbar_label"]
            title = plotting_parameters["amplitude_title"]
        elif self.phase_indicator in channel:
            cmap = plotting_parameters["phase_cmap"]
            label = plotting_parameters["phase_cbar_label"]
            title = plotting_parameters["phase_title"]
        elif self.real_indicator in channel:
            cmap = plotting_parameters["real_cmap"]
            label = plotting_parameters["real_cbar_label"]
            title = plotting_parameters["real_title_real"]
        else:
            cmap = 'viridis'
        fig, ax = plt.subplots()
        max_val = np.max(cutplane_data)
        img = plt.pcolormesh(cutplane_data, cmap=cmap, vmin=-max_val, vmax=max_val)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=f"{self.colorbar_width}%", pad=0.05) # size is the size of colorbar relative to original axis, 100% means same size, 10% means 10% of original
        cbar = plt.colorbar(img, aspect=1, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(label, rotation=270)
        if self.hide_ticks == True:
            # remove ticks on x and y axis, they only show pixelnumber anyways, better to add a scalebar
            ax.set_xticks([])
            ax.set_yticks([])
        if self.tight_layout is True:
            plt.tight_layout()
        if PlotDefinitions.show_plot is True:
            plt.show()
        gc.collect()

    def _get_z_shift_(self, z_data):
        # get the average z position for each approach curve
        # might change in the future to a more sophisticated method
        # return np.mean(z_data)

        # return the shift of the starting point of the approach curve
        return z_data[0]

    def display_approach_curve(self, x_pixel, y_pixel, x_channel:str=None, y_channels:list=None):
        if x_channel == None:
            x_channel = 'Z'
        if x_channel not in self.channels:
            print('The specified x channel is not in the channels of the measurement! Can not display approach curve.')
            return None
        if y_channels == None:
            y_channels = self.channels
        x_data = self.all_data[x_channel][y_pixel][x_pixel]
        y_data = []
        for channel in y_channels:
            y_data.append(self.all_data[channel][y_pixel][x_pixel])
        self._display_approach_curve(x_data, y_data, x_channel, y_channels)

    def _display_approach_curve(self, x_data, y_data:list, x_channel, y_channels):
        
        # x_channel = 'Depth'
        
        # import matplotlib.colors as mcolors
        # colors = mcolors.TABLEAU_COLORS
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:olive']
        fig, ax1 = plt.subplots()
        line1, = ax1.plot(x_data, y_data[0], label=y_channels[0], color=colors[0])
        if len(y_channels) == 1:
            ax1.legend()
        elif len(y_channels) == 2:
            ax2 = ax1.twinx()
            line2, = ax2.plot(x_data, y_data[1], label=y_channels[1], color=colors[1])
            ax2.set_ylabel(y_channels[1])
            ax1.legend(handles=[line1, line2])
        else: # deactivate ticks for all except the first or it will get messy
            handles = [line1]
            for channel in y_channels[1:]: # ignore the first as it was plotted already
                # i = self.channels.index(channel)
                i = y_channels.index(channel)
                # plt.plot(x_data, self.all_data[channel], label=channel)
                ax = ax1.twinx()
                ax.tick_params(right=False, labelright=False)
                line, = ax.plot(x_data, y_data[i], label=channel, color=colors[i])
                handles.append(line)
            ax1.legend(handles=handles)
            
        # print(x_data)
        # print(self.all_data[y_channels[0]])
        # print(self.channels)

        # labels for axes:
        ax1.set_xlabel(f'Z [nm]')
        ax1.set_ylabel(y_channels[0])
        # plt.xlabel(f'Depth [px]')
        # if len(self.channels) == 1:
        #     plt.ylabel(self.channels[0])
        # plt.legend()
        if PlotDefinitions.tight_layout:
            plt.tight_layout()
        
        if PlotDefinitions.show_plot:
            plt.show()
        gc.collect()

    def match_phase_offset(self, channels:list=None, reference_channel=None, reference_area=None, manual_width=5, axis='x', line=0) -> None:
        """This function matches the phase offset of all phase channels in memory to the reference channel.
        The reference channel is the first phase channel in memory if not specified.

        Args:
            channels (list, optional): list of channels, will override the already existing channels
            reference_channel ([type], optional): The reference channel to which all other phase channels will be matched.
                If not specified the first phase channel in memory will be used. Defaults to None.
            reference_area ([type], optional): The area in the reference channel which will be used to calculate the phase offset. If not specified the whole image will be used.
                You can also specify 'manual' then you will be asked to click on a point in the image. The area around that pixel will then be used as reference. Defaults to None.
            manual_width (int, optional): The width of the manual reference area. Only applies if reference_area='manual'. Defaults to 5.
        """
        # if a list of channels is specified those will be loaded and the old ones will be overwritten
        # self._initialize_data(channels)
        # define local list of channels to use for leveling
        channels = self.channels
        if reference_channel == None:
            for channel in channels:
                if self.phase_indicator in channel:
                    reference_channel = channel
                    break
        cutplane_data = self.get_cutplane_data(axis=axis, line=line, channel=reference_channel)
        if reference_area is None:
            # reference_area = [[xmin, xmax][ymin, ymax]]
            reference_area = [[0, len(cutplane_data[0])],[0, len(cutplane_data)]]
        elif reference_area == 'manual':
            # use pointcklicker to get the reference area
            fig, ax = plt.subplots()
            ax.pcolormesh(cutplane_data, cmap=SNOM_phase)
            klicker = clicker(ax, ["event"], markers=["x"])
            ax.legend()
            ax.axis('scaled')
            # ax.invert_yaxis()
            plt.title('Please click in the area to use as reference.')
            plt.show()
            klicker_coords = klicker.get_positions()['event']
            klick_coordinates = [[round(element[0]), round(element[1])] for element in klicker_coords]
            # make sure only one point is selected
            if len(klick_coordinates) != 1 and type(klick_coordinates[0]) != list:
                print('You must specify one point which should define the reference area!')
                print('Do you want to try again?')
                user_input = self._user_input_bool()
                if user_input == True:
                    self.match_phase_offset(channels, reference_channel, 'manual', manual_width, axis, line)
                else:
                    exit()
            reference_area = [[klick_coordinates[0][0] - manual_width,klick_coordinates[0][0] + manual_width],[klick_coordinates[0][1] - manual_width, klick_coordinates[0][1] + manual_width]]
        
        reference_data = cutplane_data
        reference_phase = np.mean([cutplane_data[reference_area[0][0]:reference_area[0][1]] for i in range(reference_area[1][0], reference_area[1][1])])
        
        # display the reference area
        fig, ax = plt.subplots()
        img = ax.pcolormesh(reference_data, cmap=SNOM_phase)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('phase', rotation=270)
        ax.legend()
        ax.axis('scaled')  
        rect = patches.Rectangle((reference_area[0][0], reference_area[1][0]), reference_area[0][1]-reference_area[0][0], reference_area[1][1]-reference_area[1][0], linewidth=1, edgecolor='g', facecolor='none')
        ax.add_patch(rect)
        ax.invert_yaxis()
        plt.title('Reference Area: ' + reference_channel)
        plt.show()

        for channel in channels:
            if self.phase_indicator in channel:
                # phase_data = self.get_cutplane_data(axis=axis, line=line, channel=channel)
                phase_data = self.all_cutplane_data[channel]
                # phase_offset = np.mean(phase_data) - reference_phase
                phase_offset = np.mean([phase_data[i][reference_area[0][0]:reference_area[0][1]] for i in range(reference_area[1][0], reference_area[1][1])]) - reference_phase
                self.all_cutplane_data[channel] = self._shift_phase_data(phase_data, -phase_offset)
        self._write_to_logfile('match_phase_offset_reference_area', reference_area)
        gc.collect()

    def _shift_phase_data(self, data, shift) -> np.array:
        """This function adds a phaseshift to the specified phase data. The phase data is automatically kept in the 0 to 2 pi range.
        Could in future be extended to show a live view of the phase data while it can be modified by a slider...
        e.g. by shifting the colorscale in the preview rather than the actual data..."""
        yres = len(data)
        xres = len(data[0])
        for y in range(yres):
            for x in range(xres):
                data[y][x] = (data[y][x] + shift) % (2*np.pi)
        return data

    def shift_phase(self, shift:float=None, channels:list=None) -> None:
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
        # self._initialize_data(channels)
        if shift == None:
            shift_known = False
        else:
            shift_known = True
        if shift_known is False:
            if self.preview_phasechannel in channels:
                    # phase_data = np.copy(self.all_data[self.channels.index(self.preview_phasechannel)])
                    phase_data = np.copy(self.all_cutplane_data[self.preview_phasechannel])
            else:
                # check if corrected phase channel is present
                # just take the first phase channel in memory
                for channel in channels:
                    if self.phase_indicator in channel:
                        # phase_data = np.copy(self.all_data[self.channels.index(channel)])
                        phase_data = np.copy(self.all_cutplane_data[channel])
                        # print(len(phase_data))
                        # print(len(phase_data[0]))
                        break
            shift = get_phase_offset(phase_data)
            # print('The phase shift you chose is:', shift)
            shift_known = True

        # export shift value to logfile
        self._write_to_logfile('phase_shift', shift)
        # shift all phase channels in memory
        # could also be implemented to shift each channel individually...
        
        for channel in channels:
            # print(channel)
            if self.phase_indicator in channel:
                # print('Before phase shift: ', channel)
                # print('Min phase value:', np.min(self.all_cutplane_data[channel]))
                # print('Max phase value:', np.max(self.all_cutplane_data[channel]))
                # self.all_data[self.channels.index(channel)] = self._shift_phase_data(self.all_data[self.channels.index(channel)], shift)
                self.all_cutplane_data[channel] = self._shift_phase_data(self.all_cutplane_data[channel], shift)
                # print('After phase shift: ', channel)
                # print('Min phase value:', np.min(self.all_cutplane_data[channel]))
                # print('Max phase value:', np.max(self.all_cutplane_data[channel]))
        gc.collect()

    def cut_data(self):
        pass

    def average_data(self, channels:list=None):
        if channels == None:
            channels = self.channels
        # create a cutplane of the data by averaging over the y axis
        # create a new data array with the averaged data
        self.all_cutplane_data = {}
        for channel in channels:
            if self.amp_indicator in channel:
                amp_data = self.all_data[channel]
                averaged_amp_data = np.mean(amp_data, axis=0)
                self.all_cutplane_data[channel] = np.transpose(averaged_amp_data, axes=(1,0))
            elif self.phase_indicator in channel:
                phase_data = self.all_data[channel]
                averaged_phase_data = np.mean(phase_data, axis=0)
                self.all_cutplane_data[channel] = np.transpose(averaged_phase_data, axes=(1,0))
            elif self.real_indicator in channel:
                real_data = self.all_data[channel]
                averaged_real_data = np.mean(real_data, axis=0)
                self.all_cutplane_data[channel] = np.transpose(averaged_real_data, axes=(1,0))
            elif self.height_indicator in channel:
                height_data = self.all_data[channel]
                averaged_height_data = np.mean(height_data, axis=0)
                self.all_cutplane_data[channel] = np.transpose(averaged_height_data, axes=(1,0))


        
        # averaged_height_data = np.mean(new_data, axis=2)
        # # plot the averaged height data
        # fig, ax = plt.subplots()
        # ax.pcolormesh(averaged_height_data)
        # ax.invert_yaxis()
        # plt.show()
        
    def align_lines(self):
        # idea: take the height channel and average each approach curve, then compare the averaged lines to each other and aplly a shift to align them
        height_data = self.all_data[self.height_channel]
        averaged_height_data = np.mean(height_data, axis=2)
        # plot the averaged height data
        fig, ax = plt.subplots()
        ax.pcolormesh(averaged_height_data)
        ax.invert_yaxis()
        plt.show()

        # get the index which minimized the deviation of the height channels
        indices = []
        for line in averaged_height_data:
            # calculate the index which minimizes the deviation of the height data
            index = realign.minimize_deviation_1d(averaged_height_data[0], line, 5, False)
            indices.append(index)
        # make a new data array with the shifted data
        # apply the shift to all channels
        XRes, YRes, ZRes = self._get_measurement_tag_dict_value(MeasurementTags.PIXELAREA)
        # ac_zeros = np.zeros(ZRes)
        # idea: create a new data array where each approach curve is shifted by the corresponding index
        # get the biggest differnce in indices
        max_shift = np.max(indices) - np.min(indices)
        # apply the shift to each channel
        for channel in self.channels:
            new_data = np.zeros((YRes, XRes+max_shift, ZRes))
            for y in range(YRes):
                shift = indices[y] - np.min(indices)
                for x in range(XRes):
                    new_data[y][x+shift] = self.all_data[channel][y][x]
            self.all_data[channel] = new_data
        # self.measurement_tag_dict[MeasurementTags.PIXELAREA] = (XRes+max_shift, YRes, ZRes)
        self._set_measurement_tag_dict_value(MeasurementTags.PIXELAREA, [XRes+max_shift, YRes, ZRes])


