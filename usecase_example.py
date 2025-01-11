##########################################################################
# This code was created by Hans-Joachim Schill, University of Bonn, 2022 #
##########################################################################
import tkinter as tk
from tkinter import filedialog
import pathlib
this_files_path = pathlib.Path(__file__).parent.absolute()

from src.snom_analysis.main import*

'''
This is an example script of how to access an use the snom_pyhton_classes.
Before you open a measurement you should specify the filetype and parameter type you are using.
The standard is that all data are saved in .gsf files and the parameters are extracted from a .html file.
If your data is stored in an .ascii file you can try setting the file type to that:
    File_Definitions.file_type = File_Type.aachen_ascii

If that does not work you must specify the file type yourself.

If your parameters are stored in an .parameters.txt file you can try setting the file type to that:
    File_Definitions.parmeters_type = File_Type.txt

If that does not work you must specify the parameters type yourself.

After that you can open a measurement. This will create a measurement instance.
For the measurement you must specify the measurement folder. You can additionally specify a list of channels you want to display and manipulate.
If you don't specify channels a set of standard channels will be used instead.
Furthermore you can give a measurement title, which will then be displayed in the plots.

You can then access functions on the instance to manipulate data or to plot the data.
The usable functions are:
    display_channels() # Display all specified channels or the ones in memory, plot data will be added to all subplots and can be replotted via display_all_subplots()
    scale_channels() # Scale channels by a factor, each pixel will then consist of factor*factor identical subpixels, should be applied before the blurring
    gauss_filter_channels_complex() # will blurr the complex values of the specified channels, if optical and height channels
    gauss_filter_channels() # not ideal, will just blurr all channels independently
    correct_phase_drift() # Correct linear phase drift along slow axis and apply correction to all phase channels
    shift_phase() # Shift the phase with a slider and apply the shift to all phase channels
    heigth_mask_channels() # Create a height mask fromt the height channel and applies it to all channels in memory
    fourier_filter_channels() # Applies simple fourier filter to all channels, need amplitude and phase
    cut_channels() # Cut all specified channels, either manually or auto. The height channel should be in the Measurement istance for this to make sense.
    set_min_to_zero(['MT-F-abs']) # Sets the minimum value of the specified channels to zero
    level_height_channels() # Levels the data for all height channels
    scalebar(['MT-F-abs']) # creates a scalebar in the subplot for the specified channels
    remove_subplots([1]) # Remove subplots from memory, takes a list of the indices of the subplots to delete
    remove_last_subplots(2) # Remove the last subplot from memory n-times
    switch_supplots(2, 3) # Switch positions of two subplots
    display_all_subplots() # Display all subplots from memory
    save_to_gsf() # save specified channels to .gsf files, with or without manipulations

Double underscore functions are only for class internal use.
The functions are mostly self explanatory but here are some useful tips:
Most functions will accept a channels list. You don't need to specify that, but if you do, all data in the instance memory will be overwritten with the specified channels.
Before you gauss blurr your data you should always scale the data first.
Further functions can be added as long as they don't interfere with the other functions.

Futhter improvement could include a gui and the creation of a .txt file to save all the values from the manipulations,
such that a manipulated dataset could be reproduced exactly. Also a function to open or convert .dump datafiles could be handy.
...

Have fun.
Hajo Schill

'''


def realign():
    directory_name = 'example_measurements/2022-04-29 1613 PH topol_FB_horizontal_interf_synchronize_nanoFTIR_mixedres_long'
    # Example to realign horizontal waveguides
    # channels = ['O3P', 'O3A', 'Z C']
    channels = ['O2P', 'O2A', 'Z C']
    Measurement = SnomMeasurement(directory_name, channels)
    Measurement.display_channels()
    Measurement.set_min_to_zero(['Z C'])
    Measurement.scale_channels()
    Measurement.realign()
    # Measurement.display_channels()
    Measurement.gauss_filter_channels_complex()
    Measurement.level_height_channels()
    Measurement.cut_channels()
    # Measurement.heigth_mask_channels()
    Measurement.display_channels()
    Measurement.display_all_subplots()
 
def cut_masked():
    directory_name = 'example_measurements/2022-04-25 1227 PH pentamer_840nm_s50_2'
    channels = ['O2P', 'O2A', 'Z C']
    Measurement = SnomMeasurement(directory_name, channels)
    Measurement.set_min_to_zero(['Z C'])
    # Measurement.scale_channels()
    # Measurement.gauss_filter_channels_complex()
    Measurement.heigth_mask_channels()
    # Measurement.display_channels()
    Measurement.cut_channels(autocut=True) # autocut will remove all empty lines and columns
    Measurement.display_channels()
    # Measurement.display_all_subplots()

def test_scalebar():
    directory_name = 'example_measurements/2022-04-25 1227 PH pentamer_840nm_s50_2'
    channels = ['O2P', 'O2A', 'Z C']
    # Measurement = SnomMeasurement(directory_name, channels)
    Measurement = SnomMeasurement(directory_name, channels)
    Measurement.scalebar(['Z C'], length_fraction=0.5)
    Measurement.display_channels()
    Measurement.scale_channels()
    Measurement.gauss_filter_channels_complex()
    Measurement.scalebar(['Z C'], length_fraction=0.5)
    # Measurement.cut_channels()
    Measurement.display_channels()    
    Measurement.display_all_subplots()
    
def test_phaseshift():
    # directory_name = 'example_measurements/2022-04-25 1227 PH pentamer_840nm_s50_2'
    directory_name = 'C:/Users/Hajo/sciebo/Phd/Paper/Dielectric_Waveguides/raw_data/reflection_mode/2024-05-23 161457 PH single_wv-on-wg_-45deg_thicc'
    
    # channels = ['O2P', 'O3P', 'O2A', 'Z C']
    # Plot_Definitions.full_phase_range = False
    Plot_Definitions.shared_phase_range = False
    channels = ['O2P_manipulated', 'O3P_manipulated', 'O2A', 'Z C']
    # Measurement = SnomMeasurement(directory_name, channels)
    Measurement = SnomMeasurement(directory_name, channels)
    # Measurement.scale_channels()
    # Measurement.gauss_filter_channels_complex()
    Measurement.display_channels()
    Measurement.shift_phase()
    # Measurement.shift_phase(1.96)
    Measurement.display_channels()
    # Measurement.remove_subplots([2,3,6,7])
    # Measurement.display_all_subplots()

def compare_measurements():
    channels = ['O2A', 'O2P', 'Z C']
    measurement_titles = ['measurment1: ', 'measurement2: '] # the measurment title just precedes the generic subplot titles
    Plot_Definitions.autodelete_all_subplots = False # keep subplots from previous measurement in memory!
    N = 2
    for i in range(N):
        # directory_name = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2022_07_27')
        directory_name = filedialog.askdirectory(initialdir='testdata')
        measurement_title = measurement_titles[i]
        Measurement = SnomMeasurement(directory_name, channels, measurement_title)
        Measurement.display_channels()

    Measurement.display_all_subplots()

def test_rectangle_selector():
    directory_name = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2022-04-25 1227 PH pentamer_840nm_s50_2'
    channels = ['O2P', 'O2A', 'Z C']
    Measurement = SnomMeasurement(directory_name, channels)
    Measurement.display_channels()
    # print('Select a rectangle in the plot')
    Measurement.cut_channels(reset_mask=True)
    Measurement.display_channels()
    Measurement.display_all_subplots()

def correct_phase_drift():
    directory_name = 'example_measurements/2020-03-04 1810 PH Arrow EFM Tip 9K_815nm_hires'
    channels = ['O2P', 'O3P', 'Z C']
    Measurement = SnomMeasurement(directory_name, channels)
    Measurement.display_channels()
    Measurement.shift_phase()
    Measurement.display_channels()
    Measurement.correct_phase_drift()
    Measurement.display_channels()
    Measurement.shift_phase()
    Measurement.display_channels()
    Measurement.display_all_subplots()
    
def synccorrection():
    # directory_name = 'example_measurements/2020-01-08 1337 PH denmark_skurve_02_synchronize'
    directory_name = 'example_measurements/2022-08-30 1454 PH cc_BV_No3_interf_sync_CP1R'
    channels = ['O2P', 'O2A', 'Z C']
    # channels = ['O2P_corrected', 'O2A', 'Z C']
    # channels = ['O2Re_corrected', 'O2A', 'Z C']
    Measurement = SnomMeasurement(directory_name, channels, autoscale=False)
    Measurement.synccorrection(1.6)
    # Measurement.display_channels(['O2Re_corrected', 'O2A', 'Z C'])
    Measurement.display_channels(['O2A', 'Z C'])

def complete_example_1():
    # examples from aachen group Taubner
    # directory_name = 'example_measurements/2022-07-07_16-10-33_scan_AaronLukas_2D_4x4_array' # version unknown, .dump files
    # directory_name = 'example_measurements/2018-09-10_16-44-27_scan' # version unknown, .dump files and .ascii files
    # channels = ['O2-F-abs', 'O2-F-arg', 'MT-F-abs']
    # own examples from Bonn for different software versions
    # directory_name = 'example_measurements/2018-02-09 1506 PH dt20nmwindow1' # version 1.6.3359.1
    # directory_name = 'example_measurements/2022-04-25 1227 PH pentamer_840nm_s50_2' # version 1.8.5017.0
    # directory_name = 'example_measurements/2022-08-30 1454 PH cc_BV_No3_interf_sync_CP1R' # version 1.8.5017.0
    directory_name = 'example_measurements/2024-10-24 104052 PH wg_wv_No3_15slits_pol45deg_anal90deg_fine' # version 1.10.9592.0
    channels = ['O2P', 'O2A', 'Z C']

    Measurement = SnomMeasurement(directory_name, channels)
    # Measurement.set_min_to_zero(['Z C'])
    # Measurement.display_channels()
    Measurement.scale_channels()
    Measurement.gauss_filter_channels_complex()
    # Measurement.shift_phase()
    # Measurement.heigth_mask_channels()
    Measurement.scalebar([channels[-1]])
    Measurement.display_channels()
    # Measurement.cut_channels(autocut=True) # autocut will remove all empty lines and columns
    # Measurement.display_channels()
    # Measurement.display_all_subplots()
    # Measurement._export_all_subplots()

def test_aachen_files():
    # File_Definitions.file_type = File_Type.aachen_gsf
    # File_Definitions.file_type = File_Type.aachen_ascii
    # File_Definitions.parmeters_type = File_Type.txt
    channels = ['O2-F-abs', 'O2-F-arg', 'MT-F-abs']
    # directory_name = '2018-09-10_16-44-27_scan'
    # directory_name = 'example_measurements/2022-07-07_16-10-33_scan_AaronLukas_2D_4x4_array'
    directory_name = 'example_measurements/2018-09-10_16-44-27_scan'

    Measurement = SnomMeasurement(directory_name, channels)
    Plot_Definitions.full_phase_range = False
    # Measurement.display_channels()
    # Measurement.set_min_to_zero()
    # Measurement.scale_channels()
    # Measurement.gauss_filter_channels_complex() # will blurr the complex values of the specified channels, if optical
    Measurement.scalebar(['MT-F-abs'])
    Measurement.display_channels()
    # Measurement.gauss_filter_channels() # not ideal, will just blurr all channels independently
    # Measurement.correct_phase_drift()
    # Measurement.shift_phase()
    # Measurement.heigth_mask_channels()
    # Measurement.fourier_filter_channels()
    # Measurement.cut_channels()
    # Measurement.set_min_to_zero(['MT-F-abs'])
    # Measurement.level_height_channels()
    # Measurement.display_channels()
    # Measurement.remove_subplots([1])
    # Measurement.remove_last_subplots(2)
    # Measurement.display_all_subplots()
    # Measurement.switch_supplots(2, 3)
    # Measurement.display_all_subplots()

def test_export_to_gsf():
    channels = ['O2-F-abs', 'O2-F-arg', 'MT-F-abs']
    # directory_name = 'example_measurements/2022-04-25 1227 PH pentamer_840nm_s50_2'
    # channels = ['O2P', 'O2A', 'Z C']
    directory_name = 'example_measurements/2018-09-10_16-44-27_scan'
    Measurement = SnomMeasurement(directory_name, channels)
    # Measurement.set_min_to_zero(['Z C'])
    Measurement.scale_channels()
    Measurement.gauss_filter_channels_complex()
    # Measurement.heigth_mask_channels()
    Measurement.display_channels()
    # Measurement.cut_channels(autocut=True) # autocut will remove all empty lines and columns
    # Measurement.display_channels()
    Measurement.display_all_subplots()
    Measurement.save_to_gsf()

def test_export_and_load_all_subplots():
    directory_name = 'example_measurements/2022-08-30 1454 PH cc_BV_No3_interf_sync_CP1R'
    channels = ['O2P', 'O2A', 'Z C']
    Measurement = SnomMeasurement(directory_name, channels)
    Measurement.scale_channels()
    Measurement.gauss_filter_channels_complex()
    Measurement.scalebar(['Z C'])
    Measurement.display_channels()
    Measurement.display_channels()
    Measurement.display_channels()
    Measurement.display_channels()
    Measurement.display_channels()
    # Measurement._export_all_subplots()
    # Measurement._load_all_subplots()
    # Measurement.display_all_subplots()
    # Measurement._delete_all_subplots()

def gif():
    directory_name = 'C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-03-14-ssh-snom/2024-03-14 112623 PH single_wg_lowest_long_interf_sync'
    # directory_name =  'C:/Users/Hajo/sciebo/Phd/python/SNOM/example_measurements/2020-01-08 1337 PH denmark_skurve_02_synchronize'
    # directory_name = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2022-04-25 1212 PH pentamer_840nm_s50_1'
    channels = ['O2A', 'O2P_corrected', 'Z C']
    Measurement = SnomMeasurement(directory_name, channels)
    # Measurement.display_channels()
    Measurement.rotate_90_deg(orientation='left')
    # Measurement.gauss_filter_channels_complex()
    Measurement.create_gif('O2A', 'O2P_corrected', frames=20, fps=10, dpi=100)
    # Measurement.create_gif_V2('O2A', 'O2P_corrected', 20, 10)
    # Measurement.create_gif_Old('O2A', 'O2P_corrected', 20, 10)

def test_3d_scan():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 144100 PH 3D single_wg_20mu_3d'
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    directory = 'C:/Users/hajos/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    # channels = ['Z']
    measurement = Scan3D(directory, channels)
    measurement.set_min_to_zero()
    # measurement.display_approach_curve(20, 0, 'Z', ['Z', 'O2A', 'O3A']) 
    # measurement.display_cutplane(axis='x', line=0, channel='O3A')
    Plot_Definitions.colorbar_width = 3
    # measurement.display_cutplane_V2(axis='x', line=0, channel='O2A')
    # measurement.display_cutplane_V2(axis='x', line=0, channel='O2P')
    # measurement.display_cutplane_V2_Realpart(axis='x', line=0, demodulation=2)
    # measurement.display_cutplane_V2(axis='x', line=0, channel='O3A')
    # measurement.display_cutplane_V2(axis='x', line=0, channel='O3P')
    # measurement.display_cutplane_V2_Realpart(axis='x', line=0, demodulation=3)
    measurement.generate_all_cutplane_data()
    measurement.match_phase_offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)
    measurement.display_cutplane_v3(axis='x', line=0, channel='O2P')
    measurement.display_cutplane_v3(axis='x', line=0, channel='O3P')
    measurement.shift_phase()
    measurement.display_cutplane_v3(axis='x', line=0, channel='O2P')
    measurement.display_cutplane_v3(axis='x', line=0, channel='O3P')
    measurement.display_cutplane_v3_realpart(axis='x', line=0, demodulation=2)
    measurement.display_cutplane_v3_realpart(axis='x', line=0, demodulation=3)

def phase_correction_3d_scan():
    directory = 'C:/Users/hajos/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    # channels = ['Z']
    measurement = Scan3D(directory, channels)
    measurement.set_min_to_zero()
    Plot_Definitions.colorbar_width = 3
    measurement.generate_all_cutplane_data()
    measurement.match_phase_offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)
    measurement.display_cutplane_v3(axis='x', line=0, channel='O2P')
    measurement.display_cutplane_v3(axis='x', line=0, channel='O3P')
    measurement.shift_phase()
    measurement.display_cutplane_v3(axis='x', line=0, channel='O2P')
    measurement.display_cutplane_v3(axis='x', line=0, channel='O3P')
    measurement.display_cutplane_v3_realpart(axis='x', line=0, demodulation=2)
    measurement.display_cutplane_v3_realpart(axis='x', line=0, demodulation=3)

def average_3d_scan():
    directory = 'C:/Users/hajos/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    # channels = ['Z']
    measurement = Scan3D(directory, channels)
    measurement.set_min_to_zero()
    Plot_Definitions.colorbar_width = 3
    # measurement.generate_all_cutplane_data()
    # measurement.display_cutplane_V3(axis='x', line=0, channel='O2A')
    # measurement.display_cutplane_V3(axis='x', line=0, channel='O2P')
    # measurement.display_cutplane_V3(axis='x', line=0, channel='O3A')
    # measurement.display_cutplane_V3(axis='x', line=0, channel='O3P')
    measurement.average_data()
    # measurement.display_cutplane_V3(axis='x', line=0, channel='O2A')
    # measurement.display_cutplane_V3(axis='x', line=0, channel='O2P')
    # measurement.display_cutplane_V3(axis='x', line=0, channel='O3A')
    # measurement.display_cutplane_V3(axis='x', line=0, channel='O3P')

    measurement.match_phase_offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)
    measurement.display_cutplane_v3(axis='x', line=0, channel='O2P')
    measurement.display_cutplane_v3(axis='x', line=0, channel='O3P')
    measurement.shift_phase()
    measurement.display_cutplane_v3(axis='x', line=0, channel='O2P')
    measurement.display_cutplane_v3(axis='x', line=0, channel='O3P')
    measurement.display_cutplane_v3_realpart(axis='x', line=0, demodulation=2)
    measurement.display_cutplane_v3_realpart(axis='x', line=0, demodulation=3)

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
    Plot_Definitions.amp_cbar_range = False
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
    Plot_Definitions.amp_cbar_range = False
    measurement.display_channels()
    # measurement.correct_amplitude_drift_nonlinear(channels=['O2A'], reference_area=[0, 50])
    # measurement.correct_height_drift_nonlinear(channels=['Z C'], reference_area=[20, 40])
    measurement.level_data_columnwise(channel_list=channels)
    measurement.display_channels()

def test_channel_substraction():
    initialdir = 'C:/Users/Hajo/sciebo/Phd/Paper/Dielectric_Waveguides/raw_data/reflection_mode/24-07-10/970nm'
    directory = filedialog.askdirectory(initialdir=initialdir)
    channels = ['O3P', 'O4P']
    Plot_Definitions.amp_cbar_range = False
    Plot_Definitions.colorbar_width = 3
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
    # Plot_Definitions.height_cbar_range = False
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

    Measurement = SnomMeasurement(directory_name, channels)
    Measurement.display_channels()
    # Measurement._create_default_config()
    # Measurement._print_measurement_tags()
    # Measurement._print_config()


def main():
     
    # realign()
    # cut_masked()
    # test_scalebar()
    # test_phaseshift()
    # compare_measurements()
    # test_rectangle_selector()
    # correct_phase_drift()
    # synccorrection()
    # complete_example_1()
    # test_aachen_files()
    # test_export_to_gsf()

    # test_export_and_load_all_subplots()
    # gif()
    # test_3d_scan()
    # phase_correction_3d_scan()
    # average_3d_scan()
    # test_phase_drift_correction()
    # test_amplitude_drift_correction()
    test_height_drift_correction()
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


if __name__ == '__main__':
    main()


