
.. code-block:: python
    :linenos:
    
    # Import filedialog to open a file dialog to select the measurement folder.
    from tkinter import filedialog
    
    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import SnomMeasurement, PlotDefinitions
    PlotDefinitions.colorbar_width = 2 # colorbar width for long thin measurements looks too big

    # Open a file dialog to select the measurement folder.
    directory = filedialog.askdirectory()

    # It is always a good idea to select the channels you want to use before loading the data.
    channels = ['O2P', 'O2A', 'Z C']

    # Create an instance of the SnomMeasurement class by providing the path to the measurement folder.
    # If we want to apply the synccorrection we need to set autoscale to False.
    measurement = SnomMeasurement(directory, channels, autoscale=False)

    # Plot the data without any modifications.
    measurement.display_channels()

    # Apply the synccorrection to the data. But we don't know the direction yet.
    # The interferometer sometimes goes in the wron direction.
    # This will create corrected channels and save them as .gsf with the appendix '_corrected'.
    measurement.synccorrection(1.6)

    # We want to display the corrected channels, so we reinitialize the channels with the corrected channels.
    measurement.initialize_channels(['O2P_corrected', 'O2A', 'Z C'])

    # Plot the corrected data.
    measurement.display_channels()

    # You can also compare the data before and after the modifications.
    measurement.display_all_subplots()

    