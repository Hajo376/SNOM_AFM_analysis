
.. code-block:: python
    :linenos:
    
    # Import filedialog to open a file dialog to select the measurement folder.
    from tkinter import filedialog
    
    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import SnomMeasurement 

    # Open a file dialog to select the measurement folder.
    directory = filedialog.askdirectory()

    # It is always a good idea to select the channels you want to use before loading the data.
    channels = ['O2P', 'O2A', 'Z C']

    # Create an instance of the SnomMeasurement class by providing the path to the measurement folder.
    measurement = SnomMeasurement(directory, channels)

    # Plot the data without any modifications.
    measurement.display_channels()

    # You can also add a scalebar, which will be saved to the plot memory.
    measurement.scalebar(['Z C'])

    # for phase data the level_data_columnwise function is not yet working properly,
    # if the phase data is not in the range of 0 to 2pi (drifts more than 2pi)
    # but we can correct a linear drift in the phase data before.
    measurement.correct_phase_drift_nonlinear()

    # You can also aplly corrections such as a linear 2 point based y-phase drift correction.
    # measurement.correct_phase_drift()
    # or a nonlinear dift correction, wich takes into account each line and also applies
    # a linear correction in x-direction if both sides are specified
    measurement.level_data_columnwise()

    # I always like to set the minimum of the height channel to zero.
    # This should be applied after the drift correction.
    # Otherwise the drift correction will also change the minimum.
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

    