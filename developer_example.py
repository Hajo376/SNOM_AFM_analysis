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

import tkinter as tk
from tkinter import filedialog
import pathlib
this_files_path = pathlib.Path(__file__).parent.absolute()

# from snom_analysis.main import*
from src.snom_analysis.main import SnomMeasurement, ApproachCurve, Scan3D, FileHandler
from src.snom_analysis.lib.definitions import Definitions, MeasurementTags, ChannelTags, PlotDefinitions

'''
This is an example script of how to use this package to display and manipulate example data.
The package supports SNOM and AFM data. Sofar the main functionality is for the 2D Scans.
Ther is also limited support for approach/deproach curves and 3D scans. Spectra are not supported yet.
However, the package is designed to be easily extendable, and spectra support could be added in the future.
The measurement class you decide to use will automatically try to figure out the file type and the parameters type.
If your filetype is not covered by the package defaults (based on my experience with the software) you can set the file type yourself.
The package will create a folder in the users home directory with a config.ini. This file will store the default file types and definitions.
You can add more file types and definitions to this file. The package will then use these definitions. But keep the same format.
Also these are not permanent changes, if the package does not find the config it will recreate it, so keep a copy of your definitions.

To load a measurement you can the use the corresponding class from the main file. This will only require a directory path.
I tend to use the tkinter filedialog to select the directory. This will open a file dialog and you can select the directory.
The measurement instance of the measurement class should be stored in a variable.
The class will automatically load data on initialisation. The class will try to find out the default channels to use if you don't specify any.
You can reload the data anytime by calling the >>initialize_channels<< function on the instance.
Once the data is loaded you can display the data by calling the >>display_channels<< function on the instance.
Similarly you can apply manipulations to the data by calling the corresponding functions on the instance.
If you want to review changes you made we got you covered. You can display all subplots from memory by calling >>display_all_subplots<< on the instance.
Everytime you call the >>display_channels<< function the data will be added to the memory. You can then display all subplots from memory.
You can delete the plot memory by calling >>remove_all_subplots<< on the instance.
Alternatively everytime you create a new measurement instance the plot memory will be cleared if the PlotDefinitions.autodelete_all_subplots is set to True, which is the default.
If you want to compare multiple measurements you might want to set this to False by calling PlotDefinitions.autodelete_all_subplots=False before the first measurement.
You can give a measurement title, which will then be displayed in the plots.





You can then access functions on the instance to manipulate data or to plot the data.
Common functions are:
    display_channels() # Display all specified channels or the ones in memory, plot data will be added to all subplots and can be replotted via display_all_subplots()
    scale_channels() # Scale channels by a factor, each pixel will then consist of factor*factor identical subpixels, should be applied before the blurring
    gauss_filter_channels_complex() # will blurr the complex values of the specified channels, if optical and height channels
    gauss_filter_channels() # not ideal, will just blurr all channels independently, good for amplitude channels or height channels
    correct_phase_drift() # Correct linear phase drift along slow axis and apply correction to all phase channels
    shift_phase() # Shift the phase with a slider and apply the shift to all phase channels
    heigth_mask_channels() # Create a height mask fromt the height channel and applies it to all channels in memory
    cut_channels() # Cut all specified channels, either manually or auto. The height channel should be in the Measurement istance for this to make sense.
    set_min_to_zero() # Sets the minimum value of the specified channels to zero
    heigth_level_3point() # Levels the data for all height channels
    level_data_columnwise() # new function appliable to all channels, levels the data columnwise and nonlinearly, also takes dift on the x-axis into account
    scalebar() # creates a scalebar in the subplot for the specified channels
    remove_subplots() # Remove subplots from memory, takes a list of the indices of the subplots to delete
    remove_last_subplots() # Remove the last subplot from memory n-times
    switch_supplots() # Switch positions of two subplots
    display_all_subplots() # Display all subplots from memory
    create_new_channel() # Create your own channel with arbitrary data, it can then be handeled like the other channels, but it will require the channel name, the data, the channel tag dict and a chnnel label used for plotting
    save_to_gsf() # save specified channels to .gsf files, with or without manipulations

Underscore functions are only for class internal use, but if you know what you are doing feel free to also use the other functions.
Before you gauss blurr your data you should always scale the data first.
Further functions can be added as long as they don't interfere with the other functions.
...

Have fun.
Hajo Schill

'''

def test_realign():
    directory_name = 'tests/testdata/2022-04-29 1613 PH topol_FB_horizontal_interf_synchronize_nanoFTIR_mixedres_long'
    # Example to realign horizontal waveguides
    # channels = ['O3P', 'O3A', 'Z C']
    channels = ['O2P', 'O2A', 'Z C']
    measurement = SnomMeasurement(directory_name, channels)
    # measurement.rotate_90_deg()
    measurement.display_channels()
    measurement.set_min_to_zero(['Z C'])
    measurement.scale_channels()
    measurement.gauss_filter_channels_complex()
    # measurement.realign()
    # measurement.realign(bounds=[182, 238], axis=0)
    measurement.realign(bounds=[147, 238])
    # measurement.realign(bounds=[38, 58])
    # measurement.display_channels()
    # measurement.level_height_channels_3point()
    measurement.level_height_channels_3point(coords=[[514, 175], [1763, 375], [2635, 169]])
    # measurement.cut_channels()
    measurement.cut_channels(coords=[[64, 25], [3093, 388]])
    # measurement.heigth_mask_channels()
    measurement.display_channels()
    # measurement.display_all_subplots()

def test_cut():
    directory_name = 'tests/testdata/2022-04-25 1212 PH pentamer_840nm_s50_1'
    channels = ['O2P', 'O2A', 'Z C']
    measurement = SnomMeasurement(directory_name, channels)
    # measurement.cut_channels(autocut=False) # autocut will remove all empty lines and columns
    measurement.cut_channels(coords=[[6, 8], [43, 41]]) # autocut will remove all empty lines and columns
    measurement.display_channels()

def test_cut_comsol():
    directory_name = 'tests/testdata/DLSPPW_bragg_8slits_Ex_ongold'
    channels = ['abs', 'arg', 'Z']
    # PlotDefinitions.tight_layout = False
    measurement = SnomMeasurement(directory_name, channels)
    measurement.display_channels()
    measurement.cut_channels()
    measurement.display_channels(ncols=3)
    measurement.display_all_subplots()

def test_height_masking():
    directory_name = 'tests/testdata/2022-04-25 1212 PH pentamer_840nm_s50_1'
    channels = ['O2P', 'O2A', 'Z C']
    measurement = SnomMeasurement(directory_name, channels)
    # measurement.scale_channels()
    # measurement.gauss_filter_channels_complex()
    # measurement.level_height_channels_3point(coords=[[42, 42], [84, 186], [184, 107]])
    # measurement.set_min_to_zero(['Z C'])
    # measurement.heigth_mask_channels(threshold=0.58)

    measurement.set_min_to_zero(['Z C'])
    measurement.level_height_channels_3point(coords=[[12, 9], [13, 43], [42, 34]])
    # measurement.heigth_mask_channels()
    # measurement.heigth_mask_channels(channels=['O2A'], threshold=0.45)
    measurement.heigth_mask_channels(threshold=0.45)
    measurement.display_channels()
    # measurement.cut_channels(autocut=False) # autocut will remove all empty lines and columns
    # measurement.cut_channels(coords=[[6, 8], [43, 41]]) # autocut will remove all empty lines and columns
    # measurement.display_channels()
    # measurement.display_all_subplots()

def test_scalebar():
    # directory_name = 'tests/testdata/2022-04-29 1613 PH topol_FB_horizontal_interf_synchronize_nanoFTIR_mixedres_long'
    # directory_name = 'tests/testdata/2020-01-08 1337 PH denmark_skurve_02_synchronize'
    directory_name = 'tests/testdata/2022-04-25 1212 PH pentamer_840nm_s50_1'
    channels = ['O2P', 'O2A', 'Z C']
    measurement = SnomMeasurement(directory_name, channels, autoscale=True)
    PlotDefinitions.colorbar_width = 2
    measurement.scalebar(['Z C'], length_fraction=0.5)
    measurement.display_channels()
    # measurement.scale_channels()
    # measurement.gauss_filter_channels_complex()
    # measurement.scalebar(['Z C'], length_fraction=0.5)
    # measurement.cut_channels(coords=[[20, 29], [174, 167]])
    # measurement.display_channels()    
    # measurement.display_all_subplots()
    
def test_phaseshift():
    directory_name = 'tests/testdata/2022-04-25 1212 PH pentamer_840nm_s50_1'
    channels = ['O2P', 'O2A', 'Z C']
    # PlotDefinitions.full_phase_range = False
    # PlotDefinitions.shared_phase_range = False
    measurement = SnomMeasurement(directory_name, channels)
    # measurement.scale_channels()
    # measurement.gauss_filter_channels_complex()
    measurement.display_channels()
    # measurement.shift_phase()
    measurement.shift_phase(shift=4.82)
    measurement.display_channels()
    measurement.display_all_subplots()

def compare_measurements():
    channels = ['O2A', 'O2P', 'Z C']
    measurement_titles = ['measurment1: ', 'measurement2: '] # the measurment title just precedes the generic subplot titles
    directories = ['tests/testdata/2022-04-25 1212 PH pentamer_840nm_s50_1', 'tests/testdata/2022-04-25 1227 PH pentamer_840nm_s50_2']
    PlotDefinitions.autodelete_all_subplots = False # keep subplots from previous measurement in memory!
    N = 2
    for i in range(N):
        # directory_name = filedialog.askdirectory(initialdir='tests/testdata')
        directory_name = directories[i]
        measurement_title = measurement_titles[i]
        measurement = SnomMeasurement(directory_name, channels, measurement_title)
        if i == 0:
            # just to make shure plot memory is cleared before the first measurement is displayed
            measurement._delete_all_subplots()
        measurement.display_channels()
    # display measurements
    measurement.display_all_subplots()

def correct_phase_drift():
    directory_name = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    # the linear slope is not sufficient in this case, we need a nonlinear correction to do better
    channels = ['O2P', 'O3P', 'Z C']
    measurement = SnomMeasurement(directory_name, channels)
    # measurement.display_channels()
    # measurement.shift_phase()
    measurement.display_channels()
    # measurement.correct_phase_drift(zone=2)
    # measurement.correct_phase_drift(zone=2, point_based=False)
    measurement.correct_phase_drift(phase_slope=0.016)
    measurement.display_channels()
    measurement.display_all_subplots()

def correct_phase_drift_nonlinear():
    directory_name = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    # the linear slope is not sufficient in this case, we need a nonlinear correction to do better
    channels = ['O2P', 'O3P', 'Z C']
    measurement = SnomMeasurement(directory_name, channels)
    # measurement.display_channels()
    # measurement.shift_phase()
    measurement.display_channels()
    # measurement.correct_phase_drift(zone=2)
    measurement.correct_phase_drift_nonlinear()
    measurement.display_channels()
    measurement.display_all_subplots()
    
def synccorrection():
    directory_name = 'tests/testdata/2020-01-08 1337 PH denmark_skurve_02_synchronize'
    channels = ['O2P', 'O2A', 'Z C']
    measurement = SnomMeasurement(directory_name, channels, autoscale=False)
    # measurement.synccorrection(wavelength=1.6)
    measurement.synccorrection(wavelength=1.6, phasedir=1)
    # measurement.display_channels(['O2Re_corrected', 'O2A', 'Z C'])
    measurement.initialize_channels(['O2P_corrected', 'O2A', 'Z C'])
    measurement.display_channels()

def test_aachen_files():
    # limited support for aachen files, only ascii files are supported, could not load .dump files
    channels = ['O2-F-abs', 'O2-F-arg', 'MT-F-abs']
    directory_name = 'tests/testdata/2018-09-10_16-44-27_scan'

    measurement = SnomMeasurement(directory_name, channels)
    PlotDefinitions.full_phase_range = False
    # measurement.display_channels()
    # measurement.set_min_to_zero()
    # measurement.scale_channels()
    # measurement.gauss_filter_channels_complex() # will blurr the complex values of the specified channels, if optical
    measurement.scalebar(['MT-F-abs'])
    measurement.display_channels()
    # measurement.gauss_filter_channels() # not ideal, will just blurr all channels independently
    # measurement.correct_phase_drift()
    # measurement.shift_phase()
    # measurement.heigth_mask_channels()
    # measurement.fourier_filter_channels()
    measurement.cut_channels()
    # measurement.set_min_to_zero(['MT-F-abs'])
    # measurement.level_height_channels()
    measurement.display_channels()
    # measurement.remove_subplots([1])
    # measurement.remove_last_subplots(2)
    # measurement.display_all_subplots()
    # measurement.switch_supplots(2, 3)
    # measurement.display_all_subplots()

def test_export_to_gsf():
    directory_name = 'tests/testdata/2022-04-25 1227 PH pentamer_840nm_s50_2'
    channels = ['O2P', 'O2A', 'Z C']
    measurement = SnomMeasurement(directory_name, channels)
    measurement.scale_channels()
    measurement.gauss_filter_channels_complex()
    measurement.set_min_to_zero(['Z C'])
    measurement.heigth_mask_channels(threshold=0.45)
    measurement.cut_channels(autocut=True) # autocut will remove all empty lines and columns
    measurement.display_channels()
    measurement.save_to_gsf()

def test_gif():
    directory_name =  'tests/testdata/2020-01-08 1337 PH denmark_skurve_02_synchronize'
    # directory_name = 'tests/testdata/2022-04-25 1212 PH pentamer_840nm_s50_1'
    channels = ['O2A', 'O2P_corrected', 'Z C']
    measurement = SnomMeasurement(directory_name, channels, autoscale=False)
    # measurement.synccorrection(wavelength=1.6)
    # measurement.display_channels()
    # measurement.rotate_90_deg(orientation='left')
    # measurement.gauss_filter_channels_complex()
    measurement.create_gif('O2A', 'O2P_corrected', frames=20, fps=10, dpi=100)
    # measurement.create_gif_V2('O2A', 'O2P_corrected', 20, 10)
    # measurement.create_gif_Old('O2A', 'O2P_corrected', 20, 10)

def test_3d_scan():
    directory = 'tests/testdata/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    # directory = 'tests/testdata/2025-06-02 152505 PH 3D Gold_refl_1600nm' # newer version with one more header line
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    measurement = Scan3D(directory, channels)
    measurement.set_min_to_zero()
    PlotDefinitions.colorbar_width = 3
    measurement.generate_all_cutplane_data()
    measurement.match_phase_offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)
    measurement.display_cutplanes(axis='x', line=0, channels=['O2P'])
    measurement.display_cutplanes(axis='x', line=0, channels=['O3P'])
    measurement.shift_phase()
    measurement.display_cutplanes(axis='x', line=0, channels=['O2P'])
    measurement.display_cutplanes(axis='x', line=0, channels=['O3P'])
    measurement.display_cutplane_realpart(axis='x', line=0, demodulation=2)
    measurement.display_cutplane_realpart(axis='x', line=0, demodulation=3)

def average_3d_scan():
    directory = 'tests/testdata/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    measurement = Scan3D(directory, channels)
    PlotDefinitions.colorbar_width = 3
    measurement.average_data()
    measurement.set_min_to_zero()
    measurement.display_cutplanes(axis='x', line=0, channels=['O2A', 'O2P'], auto_align=False)
    measurement.display_cutplanes(axis='x', line=0, channels=['O3A', 'O3P'])

def test_phase_drift_correction():
    directory = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    channels = ['O2P', 'O3P', 'O4P']
    measurement = SnomMeasurement(directory, channels)
    measurement.display_channels()
    measurement.correct_phase_drift_nonlinear(channels=['O2P', 'O3P', 'O4P']) # , reference_area=[0, 50]
    # measurement.level_data_columnwise(selection=[32, 169, True, True]) # optional to get rid of dirft in individual lines
    measurement.display_channels()
    # measurement.match_phase_offset(reference_channel='O2P', reference_area=[[142, 182], [587, 627]], manual_width=20)
    # measurement.shift_phase(shift=3.11)
    # measurement.display_channels()

def test_amplitude_drift_correction():
    directory = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    channels = ['O2A', 'O3A', 'O4A']
    measurement = SnomMeasurement(directory, channels)
    PlotDefinitions.amp_cbar_range = False
    measurement.display_channels()
    measurement.correct_amplitude_drift_nonlinear(channels=['O2A', 'O3A', 'O4A'], reference_area=[0, 30])
    measurement.display_channels()

def test_height_drift_correction():
    directory = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    channels = ['Z C']
    measurement = SnomMeasurement(directory, channels)
    measurement.display_channels()
    measurement.correct_height_drift_nonlinear(channels=['Z C'], reference_area=[0, 30])
    measurement.display_channels()

def test_channel_substraction():
    directory = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    channels = ['O3P', 'O4P']
    PlotDefinitions.amp_cbar_range = False
    PlotDefinitions.colorbar_width = 3
    measurement = SnomMeasurement(directory, channels)
    measurement.display_channels()
    measurement.substract_channels('O3P', 'O4P')
    measurement.shift_phase(channels=['O3P-O4P'], shift=3.22)
    measurement.display_channels()

def use_data_external_example():
    import matplotlib.pyplot as plt
    directory = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    channels = ['Z C']
    measurement = SnomMeasurement(directory, channels, autoscale=False)
    height_data = measurement.all_data[0]
    # do your own stuff using the data as a numpy array, e.g. plot the data
    plt.pcolormesh(height_data)
    plt.show()
    # you can also use the channels and measurement tag dicts to get additional information about the data
    measurement_tag_dict = measurement.measurement_tag_dict
    channel_tag_dict = measurement.channel_tag_dict[0]
    # e.g. print the content of the measurement tag dict
    print('Measurement tag dict:')
    for key, value in measurement_tag_dict.items():
        print(f'{key} = {value}')
    print('Channel tag dict:')
    # e.g. print the content of the channel tag dict
    for key, value in channel_tag_dict.items():
        print(f'{key} = {value}')
    # but you could also use the build in functions for that
    measurement.print_measurement_tag_dict()
    measurement.print_channel_tag_dict()
    
def test_level_columnwise():
    directory = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    channels = ['O2P', 'O2A', 'Z C']
    measurement = SnomMeasurement(directory, channels)
    measurement.level_data_columnwise(channel_list=channels, selection=[32, 169, True, True])
    # its expected that phase jumps lead to problems, is not yet fully implemented
    measurement.display_channels() 

def test_get_pixel_value():
    directory = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    channels = ['Z C']
    measurement = SnomMeasurement(directory, channels)
    coords = measurement.get_pixel_coordinates(channels[0])
    val = measurement.get_pixel_value('Z C', coords)
    # val = measurement.get_pixel_value('Z C')
    print(val)

def test_profile_selector():
    import matplotlib.pyplot as plt
    directory = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    channels = ['Z C']
    measurement = SnomMeasurement(directory, channels)
    profile, start, end, width = measurement.test_profile_selection()
    plt.plot(profile)
    plt.show()
    # this tests a new funtion to select arbitrary profiles from the data, it is however not yet fully implemented
    # measurement.Display_Profiles()

def test_gauss_filter_v2():
    directory = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'
    channels = ['Z C', 'O2P', 'O2A']
    # channels = ['Z C', 'O2P_corrected_aligned', 'O2A'] # should work with all kinds of phase data as long as the shape is identical
    measurement = SnomMeasurement(directory, channels)
    measurement.scale_channels()
    measurement.gauss_filter_channels_complex()
    measurement.display_channels()

def test_comsol_data():
    directory = 'tests/testdata/DLSPPW_bragg_8slits_Ex_ongold'
    # channels = ['Z C', 'O2P_corrected_aligned', 'O2A']
    measurement = SnomMeasurement(directory)
    measurement.display_channels()
    # print('All channels: ', measurement.channels)

def test_comsol_height_data():
    import numpy as np
    def create_comsol_height_data():
        height = 0.14 # in um
        wg_width = 0.3
        wg_length = 5
        grating_width = 0.3
        grating_length = 1.5
        grating_period = 0.455
        num_of_gratings = 8
        simulation_area = [501, 301] # x, y
        physical_area = [10, 6] # x, y in um
        pixel_resolution = [physical_area[0]/(simulation_area[0]-1), physical_area[1]/(simulation_area[1]-1)] # x, y in um/pixel
        waveguide_area = [[0,wg_length], [3-wg_width/2,3+wg_width/2]] # [[x_left,x_right], [y_low, y_top]] in um
        grating_positions = [wg_length - grating_width/2 + i*grating_period for i in range(num_of_gratings)]
        grating_areas = [[[pos-grating_width/2, pos+grating_width/2], [3-grating_length/2, 3+grating_length/2]] for pos in grating_positions]
        # create the height data
        height_data = np.zeros((simulation_area[1], simulation_area[0]))
        for i in range(simulation_area[0]):
            for j in range(simulation_area[1]):
                if waveguide_area[0][0] <= i*pixel_resolution[0] <= waveguide_area[0][1] and waveguide_area[1][0] <= j*pixel_resolution[1] <= waveguide_area[1][1]:
                    height_data[j,i] = height*1000 # in nm
                for grating_area in grating_areas:
                    if grating_area[0][0] <= i*pixel_resolution[0] <= grating_area[0][1] and grating_area[1][0] <= j*pixel_resolution[1] <= grating_area[1][1]:
                        height_data[j,i] = height*1000 # in nm
        # return the height data
        return height_data

    directory = 'tests/testdata/DLSPPW_bragg_8slits_Ex_ongold'
    channels = ['abs', 'arg']
    measurement = SnomMeasurement(directory, channels)
    measurement.display_channels()
    # create the height data
    height_data = create_comsol_height_data()
    height_channel_tag_dict = measurement.channel_tag_dict[0] # just copy the amp channel tag dict
    measurement.create_new_channel(height_data, 'Z', height_channel_tag_dict, 'Height')
    measurement.rotate_90_deg(orientation='right')
    measurement.display_channels()
    # amp = measurement.all_data[0]
    # height_data = measurement.all_data[2]
    measurement.display_overlay('abs', 'Z', alpha=0.2)
    # measurement.save_to_gsf(['Z'], appendix='')

def test_approach_curve():
    # directory_name = 'tests/testdata/2024-04-03 133202 PH AC topol_20mufromcoupler_right_interf_peak'
    directory_name = 'tests/testdata/2025-06-02 142837 PH AC TGQ1_test_refl_1600nm' # newer version with one more header line
    # channels = ['M1A', 'O2P', 'O2A']
    channels = ['M1A', 'O2A', 'O3A', 'O2P', 'O3P']
    measurement = ApproachCurve(directory_name, channels)
    measurement.set_min_to_zero()
    measurement.display_channels_v2()

def test_find_measurement_type():
    directory_name = 'tests/testdata/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    handler = FileHandler(directory_name)
    print(handler.measurement_type)
    print(handler.file_type)

def test_delete_data():
    # the copy has to be remade everytime you want to delete data, otherwise the data will not be deleted
    # directory_name = 'tests/testdata/test_delete/2024-03-28 164507 PH just_incoupler_square'
    directory_name = 'tests/testdata/test_delete/2025-07-18 130120 AFM gold_disc_test_1600nm'
    measurement = SnomMeasurement(directory_name)
    measurement.display_channels()
    measurement.delete_unwanted_files(mechanical_channels=True, optical_channels=True, images_folder=True, gwy_file=True)
    # measurement.cut_channels()
    # measurement.set_min_to_zero(['Z C'])
    # measurement.display_channels()

def test_regex_recognition():
    directory_name = 'tests/testdata/2022-04-25 1212 PH pentamer_840nm_s50_1'
    channels = ['O2P', 'O2A', 'Z C']
    measurement = SnomMeasurement(directory_name, channels)
    print('all filenames: ', measurement._get_all_filenames_in_directory())
    # print('prefix and suffix for channel Z C: ', measurement._get_prefix_and_suffix('Z C'))
    measurement.text_regex_file_recognition('Z C_manipulated')

#########################################
#### Examples used in documentation: ####
#########################################

# example data links:
# pentamer:
data_path_1 = 'tests/testdata/2022-04-25 1212 PH pentamer_840nm_s50_1'

# pmma wedge with amp drift:
data_path_2 = 'tests/testdata/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm'

# transmission measurement with lower parabolo in syncronized mode:
# alternative : data_path_3 = 'tests/testdata/2020-01-08 1337 PH denmark_skurve_02_synchronize'
data_path_3 = 'tests/testdata/2024-03-28 164507 PH just_incoupler_square'

# approach curve:
data_path_4 = 'tests/testdata/2024-04-03 133202 PH AC topol_20mufromcoupler_right_interf_peak'

# 3d scan:
data_path_5 = 'tests/testdata/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'


def example_snommeasurement_1():
    """This is an example of how to use the SnomMeasurement class to load, modify and save data.
    """
    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import SnomMeasurement

    # Create an instance of the SnomMeasurement class by providing the path to the measurement folder.
    # measurement = SnomMeasurement('path/to/your/measurement/folder') # for documentation purposes we will use a placeholder
    measurement = SnomMeasurement(data_path_1)

    # Now you can access the data and the functions of the measurement instance.
    # For example you can plot the data by calling the plot function.
    measurement.display_channels() # if you don't specify any arguments all channels will be plotted.

    # You can also apply some modifications to the data, for example you can crop the data.
    # If you don't specify any arguments you will be asked to provide the cropping range using
    # some default channel and the crop will be applied to all channels.
    measurement.cut_channels()

    # You can also blurr the data using a gaussian filter.
    # If you want to blurr phase data you will need to have
    # the phase and amplitude channels of the same demodulation in memory.
    measurement.gauss_filter_channels_complex()

    # You can use a simple 3 point correction to level the height data.
    # For better leveling you can always use some other software like Gwyddion and import already leveled data.
    measurement.level_height_channels_3point()

    # I always like to set the minimum of the height channel to zero.
    # This should be applied after the leveling. Otherwise the leveling will also change the minimum.
    measurement.set_min_to_zero(['Z C'])

    # Now you can plot the data again to see the modifications.
    measurement.display_channels()

    # You can also compare the data before and after the modifications.
    measurement.display_all_subplots()

    # If you are happy with the modifications you can save the data.
    # This will save the data to a new .gsf file in the measurement folder.
    measurement.save_to_gsf()
    # Per default the '_manipulated' suffix will be added to the filename,
    # to ensure you don't overwrite the original data.

def example_snommeasurement_2():
    """This is an example of how to correct a nonlinear drift in the measurements.
    """
    # Import filedialog to open a file dialog to select the measurement folder.
    # from tkinter import filedialog # for documentation purposes

    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import SnomMeasurement, PlotDefinitions
    PlotDefinitions.colorbar_width = 3


    # Open a file dialog to select the measurement folder.
    # directory = filedialog.askdirectory() # for documentation purposes
    directory = data_path_2

    # It is always a good idea to select the channels you want to use before loading the data.
    channels = ['O2P', 'O2A', 'Z C']

    # Create an instance of the SnomMeasurement class by providing the path to the measurement folder.
    measurement = SnomMeasurement(directory, channels)

    # Plot the data without any modifications.
    measurement.display_channels()

    # You can also add a scalebar, which will be saved to the plot memory.
    measurement.scalebar(['Z C'])

    # for phase data the level_data_columnwise function is not yet working properly if the phase data is not in the range of 0 to 2pi (drifts more than 2pi)
    # but we can correct a linear drift in the phase data before.
    measurement.correct_phase_drift_nonlinear()

    # You can also aplly corrections such as a linear 2 point based y-phase drift correction.
    # measurement.correct_phase_drift()
    # or a nonlinear dift correction, wich takes into account each line and also applies
    # a linear correction in x-direction if both sides are specified
    measurement.level_data_columnwise()

    # I always like to set the minimum of the height channel to zero.
    # This should be applied after the drift correction. Otherwise the drift correction will also change the minimum.
    measurement.set_min_to_zero(['Z C'])

    # You can shift the phase-offset of the phase channel.
    # Some functions will take additional parameters, if given they will apply directly.
    # If not they will promt the user with a popup window to select the parameters manually.
    measurement.shift_phase()
    # or you can specify the shift directly.
    # measurement.shift_phase(shift=0.3)

    # Now you can plot the data.
    measurement.display_channels()

    # You can also compare the data before and after the modifications.
    measurement.display_all_subplots()

def example_snommeasurement_3():
    """This is an example of how to use the synccorrection for transmission measurements with the lower parabolo in syncronized mode.
    """
    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import SnomMeasurement, PlotDefinitions
    # PlotDefinitions.colorbar_width = 4 # colorbar width for long thin measurements looks too big

    # Open a file dialog to select the measurement folder.
    # directory = filedialog.askdirectory() # for documentation purposes
    directory = data_path_3

    # It is always a good idea to select the channels you want to use before loading the data.
    channels = ['O3P']

    # Create an instance of the SnomMeasurement class by providing the path to the measurement folder.
    # If we want to apply the synccorrection we need to set autoscale to False.
    measurement = SnomMeasurement(directory, channels, autoscale=False)

    # Plot the data without any modifications.
    measurement.display_channels()

    # Apply the synccorrection to the data. But we don't know the direction yet. The interferometer sometimes goes in the wron direction.
    # This will create corrected channels and save them as .gsf with the appendix '_corrected'.
    # measurement.synccorrection(1.6) # for chiral coupler
    measurement.synccorrection(0.97)

    # We want to display the corrected channels, so we reinitialize the channels with the corrected channels.
    measurement.autoscale = True
    measurement.initialize_channels(['Z C', 'O3A', 'O3P', 'O3P_corrected'])
    
    # adjust the height channel somewhat
    measurement.level_data_columnwise(['Z C'], 'Z C')
    measurement.set_min_to_zero(['Z C'])

    # Plot the corrected data.
    measurement.display_channels()

    # You can also compare the data before and after the modifications.
    # measurement.display_all_subplots()

def example_snommeasurement_4():
    """This is an example of how to create a realpart gif.
    """
    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import SnomMeasurement, PlotDefinitions
    PlotDefinitions.colorbar_width = 4 # colorbar width for long thin measurements looks too big

    # Open a file dialog to select the measurement folder.
    # directory = filedialog.askdirectory() # for documentation purposes
    directory = data_path_3

    # Create an instance of the SnomMeasurement class by providing the path to the measurement folder.
    # If we want to apply the synccorrection we need to set autoscale to False.
    measurement = SnomMeasurement(directory, autoscale=False)

    # Plot the data without any modifications.
    measurement.display_channels()

    # Apply the synccorrection to the data. But we don't know the direction yet. The interferometer sometimes goes in the wron direction.
    # This will create corrected channels and save them as .gsf with the appendix '_corrected'.
    # measurement.synccorrection(1.6) # for chiral coupler
    measurement.synccorrection(0.97)

    # We want to use the corrected channels, so we reinitialize the channels with the corrected channels.
    measurement.autoscale = True
    measurement.initialize_channels(['O3A', 'O3P_corrected'])
    # measurement.scale_channels()
    # measurement.gauss_filter_channels_complex()

    # create a gif of the realpart of the O3A channel and the O3P_corrected channel.
    measurement.create_gif('O3A', 'O3P_corrected', frames=20, fps=10, dpi=100)
    
def example_approachcurve_1():
    # Import filedialog to open a file dialog to select the measurement folder.
    from tkinter import filedialog

    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import ApproachCurve

    # Open a file dialog to select the measurement folder.
    # directory_name = filedialog.askdirectory()

    # It is always a good idea to select the channels you want to use before loading the data.
    channels = ['M1A', 'O2A']

    # Create an instance of the ApproachCurve class by providing the path to the measurement folder.
    # measurement = ApproachCurve(directory_name, channels)
    measurement = ApproachCurve(data_path_4, channels)

    # Set the minimum value of the data to zero.
    measurement.set_min_to_zero()

    # Display the channels in a plot. And scale each data set to the whole image size.
    measurement.display_channels_v2()

    # Alternatively you can use the display_channels() function, which will display the channels in a plot without rescaling.
    # measurement.display_channels()

def example_scan3d_1():
    # Import filedialog to open a file dialog to select the measurement folder.
    from tkinter import filedialog

    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import Scan3D

    # Open a file dialog to select the measurement folder.
    # directory = filedialog.askdirectory()

    # It is always a good idea to select the channels you want to use before loading the data.
    channels = ['O2A', 'O2P', 'Z', 'O3A', 'O3P']

    # Create an instance of the Scan3D class by providing the path to the measurement folder.
    # measurement = Scan3D(directory, channels)
    measurement = Scan3D(data_path_5, channels)

    # Set the minimum value of the data to zero.
    measurement.set_min_to_zero()

    # Change the colorbar width to 3 on the fly will apply to all following plots.
    # PlotDefinitions.colorbar_width = 3

    # Generate all cutplane data for the channels.
    measurement.generate_all_cutplane_data()

    # For example display just the first cutplane of the O2P channel. Then display the first cutplane of the O3P channel.
    measurement.display_cutplanes(axis='x', line=0, channels=['O2A', 'O2P'])
    measurement.display_cutplanes(axis='x', line=0, channels=['O2A', 'O2P'], auto_align=True)

def example_scan3d_2():
    # Import filedialog to open a file dialog to select the measurement folder.
    from tkinter import filedialog

    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import Scan3D

    # Open a file dialog to select the measurement folder.
    # directory = filedialog.askdirectory()

    # It is always a good idea to select the channels you want to use before loading the data.
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']

    # Create an instance of the Scan3D class by providing the path to the measurement folder.
    # measurement = Scan3D(directory, channels)
    measurement = Scan3D(data_path_5, channels)

    # Set the minimum value of the data to zero.
    measurement.set_min_to_zero()

    # Change the colorbar width to 3 on the fly will apply to all following plots.
    PlotDefinitions.colorbar_width = 3

    # Generate all cutplane data for the channels.
    measurement.generate_all_cutplane_data()

    # Match the phase offset of the channels O2P and O3P to the reference channel O2P.
    measurement.match_phase_offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)

    # For example display just the first cutplane of the O2P channel. Then display the first cutplane of the O3P channel.
    measurement.display_cutplanes(axis='x', line=0, channels=['O2P'])
    measurement.display_cutplanes(axis='x', line=0, channels=['O3P'])

    # Shift the phase of the channels O2P and O3P by an arbitrary amount to make it visually clearer what you want to see.
    measurement.shift_phase()

    # Display the first cutplane again with all the changes applied.
    measurement.display_cutplanes(axis='x', line=0, channels=['O2P'])
    measurement.display_cutplanes(axis='x', line=0, channels=['O3P'])

    # Display the real part of the first cutplane of the O2 channel and the O3 channel.
    measurement.display_cutplane_realpart(axis='x', line=0, demodulation=2)
    measurement.display_cutplane_realpart(axis='x', line=0, demodulation=3)


def main():
    #######################################################################################
    #### Testes functions, which can be used to test the functionality of the package. ####
    #######################################################################################
    # test_realign()
    # test_cut()
    # test_cut_comsol()
    # test_height_masking()
    # test_scalebar()
    # test_phaseshift()
    # compare_measurements()
    # correct_phase_drift()
    # correct_phase_drift_nonlinear()
    # synccorrection()
    # test_aachen_files()
    # test_export_to_gsf()
    # test_gif()
    # test_3d_scan()
    # average_3d_scan()
    # test_phase_drift_correction()
    # test_amplitude_drift_correction()
    # test_height_drift_correction()
    # test_channel_substraction()
    # use_data_external_example()
    # test_level_columnwise()
    # test_get_pixel_value()
    # test_profile_selector()
    # test_gauss_filter_v2()
    # test_comsol_data()
    # test_comsol_height_data()
    # test_approach_curve()
    # test_find_measurement_type()
    test_delete_data()
    # test_regex_recognition()

    ################################
    #### Documentation examples ####
    ################################
    # example_snommeasurement_1()
    # example_snommeasurement_2()
    # example_snommeasurement_3()
    # example_snommeasurement_4()
    # example_approachcurve_1()
    # example_scan3d_1()
    # example_scan3d_2()

    ############################
    #### Untested functions ####
    ############################


if __name__ == '__main__':
    main()


