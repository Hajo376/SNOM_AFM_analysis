This package is a toolset to load, view and manipulate SNOM or AFM data.

Currently SNOM/AFM data and approach curves are supported, although the focus so far was on the standard SNOM measurements created by a standard Neaspec (Attocube) s-SNOM.
The specific class will load one measurement and allows to handle multiple channels at once.
Such a measurement instance will contain the raw measurement data as numpy arrays, the channel names, and also a dictionary containing information about the measurement and one for each channel.
The dictionaries are created automatically from the measurement parameter '.txt' file and the individual header section in each '.gsf' file.
File format support is somewhat limited, but by extending the 'config.ini' more file formats can be added without the need to edit the source code.
The config file has one section for each filetype including parameters such as the channel names and identifiers.
Once the data has been loaded it can be modified using functions, it can also be displayed and saved with the modifications.
When applying modifications to the data a log file will be created in the measurement folder and extended to save which modifications have been made. That makes it easier to track back what modifications you have made to your data once it was saved.
The package also creates a folder in the users 'home' directory to store the 'config.ini' and also a file containing the plot memory.
Each plotting command generates a plot and that plot is automatically appended to the plot memory file, allowing to view multiple measurements besides each other.
This memory can of course be deleted at will, which is by default every time a new measurement is opened.

For ease of use i also created a GUI version of this package which allows to use many functions more easily.
I typically use it to quickly generate images of my measurements with minimal modifications.
It's not yet available but i will try to make it available in the future.

I created this package during my PhD for various data analysis purposes. I want to make it available for the general public and have no commercial interests. 
I am not a professional programmer, please keep that in mind if you encounter problems or bugs, which you most probably will.
I encourage you to use and extend the functionality of this script to turn it into something more useful, which can be used by the SNOM/AFM community.

How is the package structured:
------------------------------

The package is structured in a way that the main functionality is in the 'main.py' file.
This file contains the FileHandler class which does all the file handling and contains some basic functions used by the measurement classes.
The main measurement class is the 'SnomMeasurement' class, which is the main class to load and manipulate the standard 2d data.
There are also the 'ApproachCurve' and 'Scan3D' class which are for the other types of measurements.
Spectra are not yet implemented.
There is also a PlotDefinitions class in the lib.definitions file, which contains the plot definitions for the plotting functions.
This can also be changed by the user to adapt the plotting to his needs. This will also be saved as a json file in the users home directory.
Changes made here will become the default settings for the plotting functions.
The other classes contain enums and are only for internal referencing.
Addtional functionality is in the 'lib' folder to keep the main file as small as possible.

When loading a measurement for the first time a 'config.ini' file will be created in the users home directory.
This file contains a predefined set of measurement filetypes and definitions such as the channel names and identifiers.
This file can be extended by the user to add more filetypes or change the channel names.
So far only .gsf and .txt files are supported, this could be improved in the future.
This is relatively complicated, because the measurement instance will create a dictionary from the parameters.txt file in the 
measurement folder and also from the header section in the .gsf files. One is saved per measurement and one per channel.
For example you might want to load channels with different resolutions if you modified the data.
Most functions also add a print statement to the log file in the measurement folder to keep track of what has been done to the data.
You can then use the parameters in the log file to exactly reproduce the modifications you made to the data, as most functions allow a direct parameter input.

When saving the data make shure to use an appendix to the filename, as the package will overwrite the original file if you don't.
However, the function will use a default appendix if you don't specify one.
Nevertheles, i would recommend to work on a copy of the data and not on the original data.

Also for plotting the package makes use of a matplotlib style file which is also located in the users home directory. Feel free to change this file to your needs.
You can also experiment with latex mode if you want to create plots for your thesis, then however the scalebar will not work properly since symbols like a 'µ' aren't rendered properly...
If you delete this file it will be recreated with the default settings.

Installation:
-------------

The package can be installed via ``pip``::

    pip install snom_analysis

or via ``conda``::

    conda install -c conda-forge snom_analysis

Documentation and examples:
---------------------------

The documentation can be found at https://snom-analysis.readthedocs.io/en/latest/

Key areas to improve:
---------------------

Improve the documentation. 
Some sort of Unittesting or Integration testing. 
More universal data loading. (Similar to gwyddion) 
Support for spectra and interferogramms. Maybe incorporate pysnom by Quasars. 
Better support for approach curves and 3D scans. 
Better and adjustable plotting, (adaptable layout, custom fonts, sizes...) 
Allow for custom colormaps. 
Speed improvements. 
Maybe switch to imshow instead of pcolormesh. 
And of course, more functionality.


