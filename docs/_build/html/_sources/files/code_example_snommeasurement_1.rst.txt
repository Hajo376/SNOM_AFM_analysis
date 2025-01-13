
.. code-block:: python
    :linenos:
    
    # Load the main functionality from the package, in this case the SnomMeasurement class.
    from snom_analysis.main import SnomMeasurement 

    # Create an instance of the SnomMeasurement class by providing the path to the measurement folder.
    measurement = SnomMeasurement('path/to/your/measurement/folder')

    # Now you can access the data and the functions of the measurement instance.
    # For example you can plot the data by calling the plot function.
    # If you don't specify any arguments all channels will be plotted.
    measurement.display_channels() 

    # You can also apply some modifications to the data, for example you can crop the data.
    # If you don't specify any arguments you will be asked to provide the cropping range using
    # some default channel and the crop will be applied to all channels.
    measurement.cut_channels()

    # You can also blurr the data using a gaussian filter.
    # If you want to blurr phase data you will need to have
    # the phase and amplitude channels of the same demodulation in memory.
    measurement.gauss_filter_channels_complex()

    # You can use a simple 3 point correction to level the height data.
    # For better leveling you can always use some other software
    # like Gwyddion and import already leveled data.
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