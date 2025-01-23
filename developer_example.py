##########################################################################
# This code was created by Hans-Joachim Schill, University of Bonn, 2022 #
##########################################################################

import tkinter as tk
from tkinter import filedialog
import pathlib
this_files_path = pathlib.Path(__file__).parent.absolute()

# from snom_analysis.main import*
from snom_analysis.main import SnomMeasurement, ApproachCurve, Scan3D
from snom_analysis.lib.definitions import Definitions, MeasurementTags, ChannelTags, PlotDefinitions

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
    measurement.display_channels()
    measurement.set_min_to_zero(['Z C'])
    measurement.scale_channels()
    # measurement.realign()
    measurement.realign(bounds=[131, 255])
    # measurement.display_channels()
    measurement.gauss_filter_channels_complex()
    # measurement.level_height_channels_3point()
    measurement.level_height_channels_3point(coords=[[514, 175], [1763, 375], [2635, 169]])
    # measurement.cut_channels()
    measurement.cut_channels(coords=[[64, 25], [3093, 388]])
    # measurement.heigth_mask_channels()
    measurement.display_channels()
    measurement.display_all_subplots()

def test_cut():
    directory_name = 'tests/testdata/2022-04-25 1212 PH pentamer_840nm_s50_1'
    channels = ['O2P', 'O2A', 'Z C']
    measurement = SnomMeasurement(directory_name, channels)
    # measurement.cut_channels(autocut=False) # autocut will remove all empty lines and columns
    measurement.cut_channels(coords=[[6, 8], [43, 41]]) # autocut will remove all empty lines and columns
    measurement.display_channels()

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
    directory_name = 'tests/testdata/2022-04-25 1212 PH pentamer_840nm_s50_1'
    channels = ['O2P', 'O2A', 'Z C']
    measurement = SnomMeasurement(directory_name, channels)
    measurement.scalebar(['Z C'], length_fraction=0.5)
    measurement.display_channels()
    measurement.scale_channels()
    measurement.gauss_filter_channels_complex()
    measurement.scalebar(['Z C'], length_fraction=0.5)
    measurement.cut_channels(coords=[[20, 29], [174, 167]])
    measurement.display_channels()    
    measurement.display_all_subplots()
    
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
    measurement = SnomMeasurement(directory_name, channels)
    # measurement.display_channels()
    # measurement.rotate_90_deg(orientation='left')
    # measurement.gauss_filter_channels_complex()
    measurement.create_gif('O2A', 'O2P_corrected', frames=20, fps=10, dpi=100)
    # measurement.create_gif_V2('O2A', 'O2P_corrected', 20, 10)
    # measurement.create_gif_Old('O2A', 'O2P_corrected', 20, 10)

def test_3d_scan():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 144100 PH 3D single_wg_20mu_3d'
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    # directory = 'C:/Users/hajos/git_projects/SNOM_AFM_analysis/tests/testdata/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    # directory = 'C:/Users/hajos/git_projects/SNOM_AFM_analysis/tests/testdata/2024-05-08 163342 PH 3D single_wg_20mu_3d_5ypx'
    directory = filedialog.askdirectory(initialdir='C:/Users/hajos/git_projects/SNOM_AFM_analysis/tests/testdata')
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    # channels = ['Z']
    measurement = Scan3D(directory, channels)
    measurement.set_min_to_zero()
    # measurement.display_approach_curve(20, 0, 'Z', ['Z', 'O2A', 'O3A']) 
    # measurement.display_cutplane(axis='x', line=0, channel='O3A')
    PlotDefinitions.colorbar_width = 3
    # measurement.display_cutplane_V2(axis='x', line=0, channel='O2A')
    # measurement.display_cutplane_V2(axis='x', line=0, channel='O2P')
    # measurement.display_cutplane_V2_Realpart(axis='x', line=0, demodulation=2)
    # measurement.display_cutplane_V2(axis='x', line=0, channel='O3A')
    # measurement.display_cutplane_V2(axis='x', line=0, channel='O3P')
    # measurement.display_cutplane_V2_Realpart(axis='x', line=0, demodulation=3)
    measurement.generate_all_cutplane_data()
    measurement.match_phase_offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)
    measurement.display_cutplanes(axis='x', line=0, channels=['O2P'])
    measurement.display_cutplanes(axis='x', line=0, channels=['O3P'])
    measurement.shift_phase()
    measurement.display_cutplanes(axis='x', line=0, channels=['O2P'])
    measurement.display_cutplanes(axis='x', line=0, channels=['O3P'])
    measurement.display_cutplane_realpart(axis='x', line=0, demodulation=2)
    measurement.display_cutplane_realpart(axis='x', line=0, demodulation=3)

def phase_correction_3d_scan():
    directory = 'C:/Users/hajos/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    # channels = ['Z']
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
    # directory = 'C:/Users/hajos/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    directory = filedialog.askdirectory()
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    # channels = ['Z']
    measurement = Scan3D(directory, channels)
    measurement.set_min_to_zero()
    PlotDefinitions.colorbar_width = 3
    # measurement.generate_all_cutplane_data()
    # measurement.display_cutplanes(axis='x', line=0, channels=['O2A'])
    # measurement.display_cutplanes(axis='x', line=0, channels=['O2P'])
    # measurement.display_cutplanes(axis='x', line=0, channels=['O3A'])
    # measurement.display_cutplanes(axis='x', line=0, channels=['O3P'])
    measurement.average_data()
    # measurement.display_cutplanes(axis='x', line=0, channels=['O2A'])
    # measurement.display_cutplanes(axis='x', line=0, channels=['O2P'])
    # measurement.display_cutplanes(axis='x', line=0, channels=['O3A'])
    # measurement.display_cutplanes(axis='x', line=0, channels=['O3P'])

    measurement.match_phase_offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)
    measurement.display_cutplanes(axis='x', line=0, channels=['O2P'])
    measurement.display_cutplanes(axis='x', line=0, channels=['O3P'])
    measurement.shift_phase()
    measurement.display_cutplanes(axis='x', line=0, channels=['O2P'])
    measurement.display_cutplanes(axis='x', line=0, channels=['O3P'])
    measurement.display_cutplane_realpart(axis='x', line=0, demodulation=2)
    measurement.display_cutplane_realpart(axis='x', line=0, demodulation=3)

def test_phase_drift_correction():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-23 113254 PH single_wv-on-wg_long'
    directory = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-05-23-ssh-reflection')
    channels = ['O2P', 'O3P', 'O4P']
    measurement = SnomMeasurement(directory, channels)
    measurement.display_channels()
    measurement.correct_phase_drift_nonlinear(channels=['O2P', 'O3P', 'O4P'], reference_area=[0, 50])
    measurement.display_channels()
    # measurement.match_phase_offset(channels=['O2P', 'O3P', 'O4P'], reference_channel='O2P', reference_area=[[0,50],[0,50]])
    # print('channels: ', measurement.channels)
    measurement.match_phase_offset(reference_channel='O2P', reference_area='manual', manual_width=20)
    measurement.display_channels()

def test_amplitude_drift_correction():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-23 113254 PH single_wv-on-wg_long'
    directory = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-05-23-ssh-reflection')
    channels = ['O2A', 'O3A', 'O4A']
    measurement = SnomMeasurement(directory, channels)
    PlotDefinitions.amp_cbar_range = False
    measurement.display_channels()
    # measurement.correct_amplitude_drift_nonlinear(channels=['O2A'], reference_area=[0, 50])
    measurement.correct_amplitude_drift_nonlinear(channels=['O2A', 'O3A', 'O4A'], reference_area=[140, 160])
    measurement.display_channels()

def test_height_drift_correction():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-23 113254 PH single_wv-on-wg_long'
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-10-23 113133 PH wg_wv_long_refl_15slits_No3_fine'
    directory = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-05-23-ssh-reflection')
    channels = ['Z C']
    # channels = ['O2A']
    measurement = SnomMeasurement(directory, channels)
    PlotDefinitions.amp_cbar_range = False
    measurement.display_channels()
    # measurement.correct_amplitude_drift_nonlinear(channels=['O2A'], reference_area=[0, 50])
    # measurement.correct_height_drift_nonlinear(channels=['Z C'], reference_area=[20, 40])
    measurement.level_data_columnwise(channel_list=channels)
    measurement.display_channels()

def test_channel_substraction():
    initialdir = 'C:/Users/Hajo/sciebo/Phd/Paper/Dielectric_Waveguides/raw_data/reflection_mode/24-07-10/970nm'
    directory = filedialog.askdirectory(initialdir=initialdir)
    channels = ['O3P', 'O4P']
    PlotDefinitions.amp_cbar_range = False
    PlotDefinitions.colorbar_width = 3
    measurement = SnomMeasurement(directory, channels)
    measurement.display_channels()
    measurement.substract_channels('O3P', 'O4P')
    measurement.shift_phase(channels=['O3P-O4P'])
    measurement.display_channels()

def simple_afm_example():
    directory = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-05-23-ssh-reflection')
    channels = ['Z C']
    measurement = SnomMeasurement(directory, channels, autoscale=False)
    height_data = measurement.all_data[0]
    plt.pcolormesh(height_data)
    plt.show()

    # measurement.display_channels()

def test_data_range_selector():
    '''array = [[1, 2, 3, 4, 5], 
             [6, 7, 8, 9, 10], 
             [11, 12, 13, 14, 15], 
             [16, 17, 18, 19, 20], 
             [21, 22, 23, 24, 25]]
    array = np.array(array)
    # slice array:
    # only display center three columns
    array_1 = array[:,1:4]
    print(array_1)'''
    directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2022-04-25 1227 PH pentamer_840nm_s50_2'
    channels = ['Z C']
    # channels = ['O2A']
    measurement = SnomMeasurement(directory, channels)
    measurement.test_data_range_selector()

def test_level_columnwise():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-10-24 104052 PH wg_wv_No3_15slits_pol45deg_anal90deg_fine' # for height
    directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm' # for height, amplitude and phase
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-08-21 132732 PH pmma_wedge_on_gold_thick_1600nm' # for height, amplitude and phase
    channels = ['Z C']
    # channels = ['O2A']
    # channels = ['O2P']
    # PlotDefinitions.height_cbar_range = False
    measurement = SnomMeasurement(directory, channels)
    measurement.level_data_columnwise(channel_list=channels)
    measurement.level_data_columnwise(channel_list=channels)
    measurement.display_channels() 
    # print(measurement.all_data[1])

def test_get_pixel_value():
    directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm' # for height, amplitude and phase
    channels = ['Z C']
    measurement = SnomMeasurement(directory, channels)
    coords = measurement.get_pixel_coordinates(channels[0])
    val = measurement.get_pixel_value('Z C', coords)
    # val = measurement.get_pixel_value('Z C')
    print(val)

def test_profile_selector():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm' # for height, amplitude and phase
    directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-10-24 104052 PH wg_wv_No3_15slits_pol45deg_anal90deg_fine'
    channels = ['Z C']
    measurement = SnomMeasurement(directory, channels)
    measurement.test_profile_selection()
    # measurement.Display_Profiles()

def test_gauss_filter_v2():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-07-25 114001 PH pmma_wedge_on_gold_thin_970nm' # for height, amplitude and phase
    directory = 'C:/Users/Hajo/sciebo/Phd/Paper/Dielectric_Waveguides/raw_data/transmission_mode/reflection_grating/2024-10-24 100547 PH wg_wv_No3_15slits_pol45deg_anal0deg_fine'
    # channels = ['Z C', 'O2P', 'O2A']
    channels = ['Z C', 'O2P_corrected_aligned', 'O2A']

    measurement = SnomMeasurement(directory, channels)
    measurement.gauss_filter_channels_complex()
    measurement.display_channels()

def test_comsol_data():
    directory = 'example_measurements/DLSPPW_bragg_8slits_Ex_ongold'
    # channels = ['Z C', 'O2P_corrected_aligned', 'O2A']
    measurement = SnomMeasurement(directory)
    measurement.display_channels()
    # print('All channels: ', measurement.channels)

def test_comsol_height_data():
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

    channels = ['abs', 'arg']
    directory = 'example_measurements/DLSPPW_bragg_8slits_Ex_ongold'
    measurement = SnomMeasurement(directory, channels)
    measurement.display_channels()
    # create the height data
    height_data = create_comsol_height_data()
    height_channel_tag_dict = measurement.channel_tag_dict[0] # just copy the amp channel tag dict
    measurement.create_new_channel(height_data, 'Z', height_channel_tag_dict, 'Height')
    measurement.display_channels()
    amp = measurement.all_data[0]
    height_data = measurement.all_data[2]
    measurement.display_overlay('abs', 'Z', alpha=0.2)
    # measurement.save_to_gsf(['Z'], appendix='')

def test_config():
    # directory_name = 'example_measurements/2018-02-09 1506 PH dt20nmwindow1' # version 1.6.3359.1
    # directory_name = 'example_measurements/2022-04-25 1227 PH pentamer_840nm_s50_2' # version 1.8.5017.0
    # directory_name = 'example_measurements/2022-08-30 1454 PH cc_BV_No3_interf_sync_CP1R' # version 1.8.5017.0
    # directory_name = 'example_measurements/2024-10-24 104052 PH wg_wv_No3_15slits_pol45deg_anal90deg_fine' # version 1.10.9592.0
    # channels = ['O2P', 'O2A', 'Z C']
    channels = ['abs', 'arg']
    directory_name = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo')
    # directory_name = 'example_measurements/DLSPPW_bragg_8slits_Ex_ongold' # comsol

    measurement = SnomMeasurement(directory_name, channels)
    measurement.display_channels()
    # measurement._create_default_config()
    # measurement._print_measurement_tags()
    # measurement._print_config()

def test_approach_curve():
    directory_name = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo')
    channels = ['M1A', 'O2P', 'O2A']
    measurement = ApproachCurve(directory_name, channels)
    measurement.set_min_to_zero()
    measurement.display_channels_v2()

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
    measurement.level_height_channels()

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
    from snom_analysis.main import SnomMeasurement

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
    PlotDefinitions.colorbar_width = 4 # colorbar width for long thin measurements looks too big

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
    measurement.initialize_channels(['O2A', 'O2P_corrected'])

    # create a gif of the realpart of the O2A channel and the O2P_corrected channel.
    measurement.create_gif('O2A', 'O2P_corrected', frames=20, fps=10, dpi=100)
    
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
    PlotDefinitions.colorbar_width = 3

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
    channels = ['O2A', 'O2P', 'Z']

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
    measurement.display_cutplane_v3(axis='x', line=0, channel='O2P')
    measurement.display_cutplane_v3(axis='x', line=0, channel='O3P')

    # Shift the phase of the channels O2P and O3P by an arbitrary amount to make it visually clearer what you want to see.
    measurement.shift_phase()

    # Display the first cutplane again with all the changes applied.
    measurement.display_cutplane_v3(axis='x', line=0, channel='O2P')
    measurement.display_cutplane_v3(axis='x', line=0, channel='O3P')

    # Display the real part of the first cutplane of the O2 channel and the O3 channel.
    measurement.display_cutplane_v3_realpart(axis='x', line=0, demodulation=2)
    measurement.display_cutplane_v3_realpart(axis='x', line=0, demodulation=3)


def main():
     
    # test_realign()
    # test_cut()
    # test_height_masking()
    test_scalebar()
    # test_phaseshift()
    # compare_measurements()
    # correct_phase_drift()
    # correct_phase_drift_nonlinear()
    # synccorrection()
    # test_aachen_files()
    # test_export_to_gsf()
    # test_gif()
    # test_3d_scan()
    # phase_correction_3d_scan()
    # average_3d_scan()
    # test_phase_drift_correction()
    # test_amplitude_drift_correction()
    # test_height_drift_correction()
    # test_channel_substraction()
    # simple_afm_example()
    # test_data_range_selector()
    # test_level_columnwise()
    # test_get_pixel_value()
    # test_profile_selector()
    # test_gauss_filter_v2()
    # test_comsol_data()
    # test_comsol_height_data()
    # test_config()
    # test_approach_curve()

    # examples for documentation
    # example_snommeasurement_1()
    # example_snommeasurement_2()
    # example_snommeasurement_3()
    # example_snommeasurement_4()
    # example_approachcurve_1()
    # example_scan3d_1()
    # example_scan3d_2()


if __name__ == '__main__':
    main()


