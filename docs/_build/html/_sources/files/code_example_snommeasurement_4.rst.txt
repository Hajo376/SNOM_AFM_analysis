
.. code-block:: python
    :linenos:
    
    # Import filedialog to open a file dialog to select the measurement folder.
    from tkinter import filedialog

    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import SnomMeasurement, PlotDefinitions
    PlotDefinitions.colorbar_width = 4 # colorbar width for long thin measurements looks too big

    # Open a file dialog to select the measurement folder.
    directory_name = filedialog.askdirectory()

    # Create an instance of the SnomMeasurement class by providing the path to the measurement folder.
    # If we want to apply the synccorrection we need to set autoscale to False.
    measurement = SnomMeasurement(directory_name, autoscale=False)

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
    