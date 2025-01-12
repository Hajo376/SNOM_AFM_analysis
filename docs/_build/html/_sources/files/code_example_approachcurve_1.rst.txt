
.. code-block:: python
    :linenos:
    
    # Import filedialog to open a file dialog to select the measurement folder.
    from tkinter import filedialog
    
    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import ApproachCurve 

    # Open a file dialog to select the measurement folder.
    directory = filedialog.askdirectory()

    # It is always a good idea to select the channels you want to use before loading the data.
    channels = ['M1A', 'O2P', 'O2A']

    # Create an instance of the ApproachCurve class by providing the path to the measurement folder.
    measurement = ApproachCurve(directory_name, channels)

    # Set the minimum value of the data to zero.
    measurement.set_min_to_zero()

    # Display the channels in a plot. And scale each data set to the whole image size.
    measurement.display_channels_v2()
    
    # Alternatively you can use the display_channels() function, which will display the channels in a plot without rescaling.
    measurement.display_channels()

    