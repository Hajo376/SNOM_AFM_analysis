
.. code-block:: python
    :linenos:
    
    # Import filedialog to open a file dialog to select the measurement folder.
    from tkinter import filedialog
    
    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import Scan3D 

    # Open a file dialog to select the measurement folder.
    directory = filedialog.askdirectory()

    # It is always a good idea to select the channels you want to use before loading the data.
    channels = ['O2A', 'O2P', 'O3A', 'O3P', 'Z']
    
    # Create an instance of the Scan3D class by providing the path to the measurement folder.
    measurement = Scan3D(directory, channels)

    # Set the minimum value of the data to zero.
    measurement.set_min_to_zero()

    # Change the colorbar width to 3 on the fly will apply to all following plots.
    PlotDefinitions.colorbar_width = 3

    # Average all individual cutplanes of the channels O2P and O3P.
    # This is an alternative approach to generating individual cutplanes.
    # All other functions will work with the averaged data, just make shure to use the 'line=0' parameter.
    measurement.average_data()

    # Match the phase offset of the channels O2P and O3P to the reference channel O2P.
    measurement.match_phase_offset(channels=['O2P', 'O3P'], reference_channel='O2P', reference_area='manual', manual_width=3)

    # For example display just the first cutplane of the O2P channel. Then display the first cutplane of the O3P channel.
    measurement.display_cutplane_v3(axis='x', line=0, channel='O2P')
    measurement.display_cutplane_v3(axis='x', line=0, channel='O3P')


    