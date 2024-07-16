##########################################################################
# This code was created by Hans-Joachim Schill, University of Bonn, 2022 #
##########################################################################
import tkinter as tk
from tkinter import filedialog
import pathlib
this_files_path = pathlib.Path(__file__).parent.absolute()

from SNOM_AFM_analysis.python_classes_snom import*

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
    Display_Channels() # Display all specified channels or the ones in memory, plot data will be added to all subplots and can be replotted via Display_All_Subplots()
    Scale_Channels() # Scale channels by a factor, each pixel will then consist of factor*factor identical subpixels, should be applied before the blurring
    Gauss_Filter_Channels_complex() # will blurr the complex values of the specified channels, if optical and height channels
    Gauss_Filter_Channels() # not ideal, will just blurr all channels independently
    Correct_Phase_Drift() # Correct linear phase drift along slow axis and apply correction to all phase channels
    Shift_Phase() # Shift the phase with a slider and apply the shift to all phase channels
    Heigth_Mask_Channels() # Create a height mask fromt the height channel and applies it to all channels in memory
    Fourier_Filter_Channels() # Applies simple fourier filter to all channels, need amplitude and phase
    Cut_Channels() # Cut all specified channels, either manually or auto. The height channel should be in the Measurement istance for this to make sense.
    Set_Min_to_Zero(['MT-F-abs']) # Sets the minimum value of the specified channels to zero
    Level_Height_Channels() # Levels the data for all height channels
    Scalebar(['MT-F-abs']) # creates a scalebar in the subplot for the specified channels
    Remove_Subplots([1]) # Remove subplots from memory, takes a list of the indices of the subplots to delete
    Remove_Last_Subplots(2) # Remove the last subplot from memory n-times
    Switch_Supplots(2, 3) # Switch positions of two subplots
    Display_All_Subplots() # Display all subplots from memory
    Save_to_gsf() # save specified channels to .gsf files, with or without manipulations

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


def Realign():
    directory_name = 'testdata/2022-04-29 1613 PH topol_FB_horizontal_interf_synchronize_nanoFTIR_mixedres_long'
    # Example to realign horizontal waveguides
    # channels = ['O3P', 'O3A', 'Z C']
    channels = ['O2P', 'O2A', 'Z C']
    Measurement = Open_Measurement(directory_name, channels)
    Measurement.Display_Channels()
    Measurement.Set_Min_to_Zero(['Z C'])
    Measurement.Scale_Channels()
    Measurement.Realign()
    # Measurement.Display_Channels()
    Measurement.Gauss_Filter_Channels_complex()
    Measurement.Level_Height_Channels()
    Measurement.Cut_Channels()
    # Measurement.Heigth_Mask_Channels()
    Measurement.Display_Channels()
    Measurement.Display_All_Subplots()
 
def Cut_Masked():
    directory_name = 'testdata/2022-04-25 1227 PH pentamer_840nm_s50_2'
    channels = ['O2P', 'O2A', 'Z C']
    Measurement = Open_Measurement(directory_name, channels)
    Measurement.Set_Min_to_Zero(['Z C'])
    Measurement.Scale_Channels()
    Measurement.Gauss_Filter_Channels_complex()
    Measurement.Heigth_Mask_Channels()
    Measurement.Display_Channels()
    Measurement.Cut_Channels(autocut=True) # autocut will remove all empty lines and columns
    Measurement.Display_Channels()
    Measurement.Display_All_Subplots()

def Test_Scalebar():
    directory_name = 'testdata/2022-04-25 1227 PH pentamer_840nm_s50_2'
    channels = ['O2P', 'O2A', 'Z C']
    Measurement = Open_Measurement(directory_name, channels)
    Measurement.Scalebar(['Z C'], length_fraction=0.5)
    Measurement.Display_Channels()
    Measurement.Scale_Channels()
    Measurement.Gauss_Filter_Channels_complex()
    Measurement.Scalebar(['Z C'], length_fraction=0.5)
    Measurement.Cut_Channels()
    Measurement.Display_Channels()    
    Measurement.Display_All_Subplots()
    
def Test_Phaseshift():
    # directory_name = 'testdata/2022-04-25 1227 PH pentamer_840nm_s50_2'
    directory_name = 'C:/Users/Hajo/sciebo/Phd/Paper/Dielectric_Waveguides/raw_data/reflection_mode/2024-05-23 161457 PH single_wv-on-wg_-45deg_thicc'
    
    # channels = ['O2P', 'O3P', 'O2A', 'Z C']
    # Plot_Definitions.full_phase_range = False
    Plot_Definitions.shared_phase_range = False
    channels = ['O2P_manipulated', 'O3P_manipulated', 'O2A', 'Z C']
    # Measurement = Open_Measurement(directory_name, channels)
    Measurement = SnomMeasurement(directory_name, channels)
    # Measurement.Scale_Channels()
    # Measurement.Gauss_Filter_Channels_complex()
    Measurement.Display_Channels()
    Measurement.Shift_Phase()
    # Measurement.Shift_Phase(1.96)
    Measurement.Display_Channels()
    # Measurement.Remove_Subplots([2,3,6,7])
    # Measurement.Display_All_Subplots()

def Compare_Measurements():
    channels = ['O2A', 'O2P', 'Z C']
    measurement_titles = ['measurment1: ', 'measurement2: '] # the measurment title just precedes the generic subplot titles
    File_Definitions.autodelete_all_subplots = False # keep subplots from previous measurement in memory!
    N = 2
    for i in range(N):
        # directory_name = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2022_07_27')
        directory_name = filedialog.askdirectory(initialdir='testdata')
        measurement_title = measurement_titles[i]
        Measurement = Open_Measurement(directory_name, channels, measurement_title)
        Measurement.Display_Channels()

    Measurement.Display_All_Subplots()

def Test_Rectangle_Selector():
    directory_name = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2022-04-25 1227 PH pentamer_840nm_s50_2'
    channels = ['O2P', 'O2A', 'Z C']
    Measurement = Open_Measurement(directory_name, channels)
    Measurement.Display_Channels()
    Measurement.Cut_Channels(reset_mask=True)
    Measurement.Display_Channels()
    Measurement.Display_All_Subplots()

def Correct_Phase_Drift():
    directory_name = 'testdata/2020-03-04 1810 PH Arrow EFM Tip 9K_815nm_hires'
    channels = ['O2P', 'O3P', 'Z C']
    Measurement = Open_Measurement(directory_name, channels)
    Measurement.Display_Channels()
    Measurement.Shift_Phase()
    Measurement.Display_Channels()
    Measurement.Correct_Phase_Drift()
    Measurement.Display_Channels()
    Measurement.Shift_Phase()
    Measurement.Display_Channels()
    Measurement.Display_All_Subplots()
    
def Synccorrection():
    # directory_name = 'testdata/2020-01-08 1337 PH denmark_skurve_02_synchronize'
    directory_name = 'example_measurements/2022-08-30 1454 PH cc_BV_No3_interf_sync_CP1R'
    channels = ['O2P', 'O2A', 'Z C']
    # channels = ['O2P_corrected', 'O2A', 'Z C']
    channels = ['O2Re_corrected', 'O2A', 'Z C']
    Measurement = Open_Measurement(directory_name, channels, autoscale=False)
    # Measurement.Synccorrection(1.6)
    Measurement.Display_Channels(['O2Re_corrected', 'O2A', 'Z C'])

def Complete_Example_1():
    # directory_name = 'example_measurements/2022-04-25 1227 PH pentamer_840nm_s50_2'
    directory_name = 'example_measurements/2022-08-30 1454 PH cc_BV_No3_interf_sync_CP1R'
    channels = ['O2P', 'O2A', 'Z C']
    Measurement = Open_Measurement(directory_name, channels)
    # Measurement.Set_Min_to_Zero(['Z C'])
    # Measurement.Display_Channels()
    Measurement.Scale_Channels()
    Measurement.Gauss_Filter_Channels_complex()
    # Measurement.Shift_Phase()
    # Measurement.Heigth_Mask_Channels()
    Measurement.Scalebar(['Z C'])
    Measurement.Display_Channels()
    # Measurement.Cut_Channels(autocut=True) # autocut will remove all empty lines and columns
    # Measurement.Display_Channels()
    # Measurement.Display_All_Subplots()
    Measurement._Export_All_Subplots()

def Test_Aachen_files():
    # File_Definitions.file_type = File_Type.aachen_gsf
    # File_Definitions.file_type = File_Type.aachen_ascii
    # File_Definitions.parmeters_type = File_Type.txt
    channels = ['O2-F-abs', 'O2-F-arg', 'MT-F-abs']
    # directory_name = '2018-09-10_16-44-27_scan'
    # directory_name = 'testdata/2022-07-07_16-10-33_scan_AaronLukas_2D_4x4_array'
    directory_name = 'example_measurements/2018-09-10_16-44-27_scan'

    Measurement = Open_Measurement(directory_name, channels)
    Plot_Definitions.full_phase_range = False
    # Measurement.Display_Channels()
    # Measurement.Set_Min_to_Zero()
    # Measurement.Scale_Channels()
    # Measurement.Gauss_Filter_Channels_complex() # will blurr the complex values of the specified channels, if optical
    Measurement.Scalebar(['MT-F-abs'])
    Measurement.Display_Channels()
    # Measurement.Gauss_Filter_Channels() # not ideal, will just blurr all channels independently
    # Measurement.Correct_Phase_Drift()
    # Measurement.Shift_Phase()
    # Measurement.Heigth_Mask_Channels()
    # Measurement.Fourier_Filter_Channels()
    # Measurement.Cut_Channels()
    # Measurement.Set_Min_to_Zero(['MT-F-abs'])
    # Measurement.Level_Height_Channels()
    # Measurement.Display_Channels()
    # Measurement.Remove_Subplots([1])
    # Measurement.Remove_Last_Subplots(2)
    # Measurement.Display_All_Subplots()
    # Measurement.Switch_Supplots(2, 3)
    # Measurement.Display_All_Subplots()

def Test_Export_to_gsf():
    # File_Definitions.file_type = File_Type.aachen_gsf
    File_Definitions.file_type = File_Type.aachen_ascii
    File_Definitions.parmeters_type = File_Type.txt
    channels = ['O2-F-abs', 'O2-F-arg', 'MT-F-abs']
    # directory_name = 'testdata/2022-04-25 1227 PH pentamer_840nm_s50_2'
    # channels = ['O2P', 'O2A', 'Z C']
    directory_name = 'testdata/2018-09-10_16-44-27_scan'
    Measurement = Open_Measurement(directory_name, channels)
    # Measurement.Set_Min_to_Zero(['Z C'])
    Measurement.Scale_Channels()
    Measurement.Gauss_Filter_Channels_complex()
    # Measurement.Heigth_Mask_Channels()
    Measurement.Display_Channels()
    # Measurement.Cut_Channels(autocut=True) # autocut will remove all empty lines and columns
    # Measurement.Display_Channels()
    Measurement.Display_All_Subplots()
    Measurement.Save_to_gsf()

def Test_Export_and_Load_all_subplots():
    directory_name = 'example_measurements/2022-08-30 1454 PH cc_BV_No3_interf_sync_CP1R'
    channels = ['O2P', 'O2A', 'Z C']
    Measurement = Open_Measurement(directory_name, channels)
    Measurement.Scale_Channels()
    Measurement.Gauss_Filter_Channels_complex()
    Measurement.Scalebar(['Z C'])
    Measurement.Display_Channels()
    Measurement.Display_Channels()
    Measurement.Display_Channels()
    Measurement.Display_Channels()
    Measurement.Display_Channels()
    # Measurement._Export_All_Subplots()
    # Measurement._Load_All_Subplots()
    # Measurement.Display_All_Subplots()
    # Measurement._Delete_All_Subplots()

def Gif():
    directory_name = 'C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-03-14-ssh-snom/2024-03-14 112623 PH single_wg_lowest_long_interf_sync'
    # directory_name =  'C:/Users/Hajo/sciebo/Phd/python/SNOM/testdata/2020-01-08 1337 PH denmark_skurve_02_synchronize'
    # directory_name = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2022-04-25 1212 PH pentamer_840nm_s50_1'
    channels = ['O2A', 'O2P_corrected', 'Z C']
    Measurement = SnomMeasurement(directory_name, channels)
    # Measurement.Display_Channels()
    Measurement.Rotate_90_deg(orientation='left')
    # Measurement.Gauss_Filter_Channels_complex()
    Measurement.Create_Gif('O2A', 'O2P_corrected', frames=20, fps=10, dpi=100)
    # Measurement.Create_Gif_V2('O2A', 'O2P_corrected', 20, 10)
    # Measurement.Create_Gif_Old('O2A', 'O2P_corrected', 20, 10)

def Test_3D_Scan():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 144100 PH 3D single_wg_20mu_3d'
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-08 151547 PH 3D single_wg_20mu_3d_10ypx'
    directory = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-05-14-ssh-snom')
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    # channels = ['Z']
    measurement = Scan_3D(directory, channels)
    measurement.Set_Min_to_Zero()
    # measurement.Display_Approach_Curve(20, 0, 'Z', ['Z', 'O2A', 'O3A']) 
    # measurement.Display_Cutplane(axis='x', line=0, channel='O3A')
    Plot_Definitions.colorbar_width = 3
    # measurement.Display_Cutplane_V2(axis='x', line=0, channel='O2A')
    # measurement.Display_Cutplane_V2(axis='x', line=0, channel='O2P')
    # measurement.Display_Cutplane_V2_Realpart(axis='x', line=0, demodulation=2)
    # measurement.Display_Cutplane_V2(axis='x', line=0, channel='O3A')
    # measurement.Display_Cutplane_V2(axis='x', line=0, channel='O3P')
    # measurement.Display_Cutplane_V2_Realpart(axis='x', line=0, demodulation=3)
    measurement.Generate_All_Cutplane_Data()
    measurement.Match_Phase_Offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O2P')
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O3P')
    measurement.Shift_Phase()
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O2P')
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O3P')
    measurement.Display_Cutplane_V3_Realpart(axis='x', line=0, demodulation=2)
    measurement.Display_Cutplane_V3_Realpart(axis='x', line=0, demodulation=3)

def Phase_Correction_3D_Scan():
    directory = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-05-14-ssh-snom')
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    # channels = ['Z']
    measurement = Scan_3D(directory, channels)
    measurement.Set_Min_to_Zero()
    Plot_Definitions.colorbar_width = 3
    measurement.Generate_All_Cutplane_Data()
    measurement.Match_Phase_Offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O2P')
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O3P')
    measurement.Shift_Phase()
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O2P')
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O3P')
    measurement.Display_Cutplane_V3_Realpart(axis='x', line=0, demodulation=2)
    measurement.Display_Cutplane_V3_Realpart(axis='x', line=0, demodulation=3)

def Average_3D_Scan():
    directory = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-05-14-ssh-snom')
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    # channels = ['Z']
    measurement = Scan_3D(directory, channels)
    measurement.Set_Min_to_Zero()
    Plot_Definitions.colorbar_width = 3
    # measurement.Generate_All_Cutplane_Data()
    # measurement.Display_Cutplane_V3(axis='x', line=0, channel='O2A')
    # measurement.Display_Cutplane_V3(axis='x', line=0, channel='O2P')
    # measurement.Display_Cutplane_V3(axis='x', line=0, channel='O3A')
    # measurement.Display_Cutplane_V3(axis='x', line=0, channel='O3P')
    measurement.Average_Data()
    # measurement.Display_Cutplane_V3(axis='x', line=0, channel='O2A')
    # measurement.Display_Cutplane_V3(axis='x', line=0, channel='O2P')
    # measurement.Display_Cutplane_V3(axis='x', line=0, channel='O3A')
    # measurement.Display_Cutplane_V3(axis='x', line=0, channel='O3P')

    measurement.Match_Phase_Offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O2P')
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O3P')
    measurement.Shift_Phase()
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O2P')
    measurement.Display_Cutplane_V3(axis='x', line=0, channel='O3P')
    measurement.Display_Cutplane_V3_Realpart(axis='x', line=0, demodulation=2)
    measurement.Display_Cutplane_V3_Realpart(axis='x', line=0, demodulation=3)

    

def Test_Phase_Drift_Correction():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-23 113254 PH single_wv-on-wg_long'
    directory = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-05-23-ssh-reflection')
    channels = ['O2P', 'O3P', 'O4P']
    measurement = SnomMeasurement(directory, channels)
    measurement.Display_Channels()
    measurement.Correct_Phase_Drift_Nonlinear(channels=['O2P', 'O3P', 'O4P'], reference_area=[0, 50])
    measurement.Display_Channels()
    # measurement.Match_Phase_Offset(channels=['O2P', 'O3P', 'O4P'], reference_channel='O2P', reference_area=[[0,50],[0,50]])
    measurement.Match_Phase_Offset(channels=['O2P', 'O3P', 'O4P'], reference_channel='O2P', reference_area='manual', manual_width=20)
    measurement.Display_Channels()

def Test_Amplitude_Drift_Correction():
    # directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-23 113254 PH single_wv-on-wg_long'
    directory = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-05-23-ssh-reflection')
    channels = ['O2A', 'O3A', 'O4A']
    measurement = SnomMeasurement(directory, channels)
    Plot_Definitions.amp_cbar_range = False
    measurement.Display_Channels()
    # measurement.Correct_Amplitude_Drift_Nonlinear(channels=['O2A'], reference_area=[0, 50])
    measurement.Correct_Amplitude_Drift_Nonlinear(channels=['O2A', 'O3A', 'O4A'], reference_area=[140, 160])
    measurement.Display_Channels()

def Test_Height_Drift_Correction():
    directory = 'C:/Users/Hajo/git_projects/SNOM_AFM_analysis/example_measurements/2024-05-23 113254 PH single_wv-on-wg_long'
    # directory = filedialog.askdirectory(initialdir='C:/Users/Hajo/sciebo/Exchange/s-SNOM Measurements/Hajo/PhD/ssh/2024-05-23-ssh-reflection')
    channels = ['Z C']
    measurement = SnomMeasurement(directory, channels)
    Plot_Definitions.amp_cbar_range = False
    measurement.Display_Channels()
    # measurement.Correct_Amplitude_Drift_Nonlinear(channels=['O2A'], reference_area=[0, 50])
    measurement.Correct_Height_Drift_Nonlinear(channels=['Z C'], reference_area=[20, 40])
    measurement.Display_Channels()

def Test_Channel_Substraction():
    initialdir = 'C:/Users/Hajo/sciebo/Phd/Paper/Dielectric_Waveguides/raw_data/reflection_mode/24-07-10/970nm'
    directory = filedialog.askdirectory(initialdir=initialdir)
    channels = ['O3P', 'O4P']
    Plot_Definitions.amp_cbar_range = False
    Plot_Definitions.colorbar_width = 3
    measurement = SnomMeasurement(directory, channels)
    measurement.Display_Channels()
    measurement.Substract_Channels('O3P', 'O4P')
    measurement.Shift_Phase(channels=['O3P-O4P'])
    measurement.Display_Channels()

def main():
     
    # Realign()
    # Cut_Masked()
    # Test_Scalebar()
    # Test_Phaseshift()
    # Compare_Measurements()
    # Test_Rectangle_Selector()
    # Correct_Phase_Drift()
    # Synccorrection()
    # Complete_Example_1()
    # Test_Aachen_files()
    # Test_Export_to_gsf()

    # Test_Export_and_Load_all_subplots()
    # Gif()
    # Test_3D_Scan()
    # Phase_Correction_3D_Scan()
    # Average_3D_Scan()
    # Test_Phase_Drift_Correction()
    # Test_Amplitude_Drift_Correction()
    # Test_Height_Drift_Correction()
    Test_Channel_Substraction()


if __name__ == '__main__':
    main()


