# This package is a toolset to load, view and manipulate SNOM or AFM data. #
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

For ease of use i also created a GUI version of this package which allows to use many functions more easily. I typically use it to quickly generate images of my measurements with minimal modifications.

I created this package during my PhD for various data analysis purposes. I want to make it available for the general public and have no commercial interests. 
I am not a professional programmer, please keep that in mind if you encounter problems or bugs, which you most probably will.
I encourage you to use and extend the functionality of this script to turn it into something useful which can be used by the SNOM/AFM community.

## Key areas to improve: ##
Documentation is missing.
Some sort of Unittesting or Integration testing.
More universal data loading. (Similar to gwyddion)
Support for spectra and interferogramms. Maybe incorporate pysnom by Quasars.
Better support for approach curves and 3D scans.
Better and adjustable plotting, (adaptable layout, custom fonts, sizes...)
Allow for custom colormaps.
Speed improvements.
Maybe switch to imshow instead of pcolormesh.

## Documentation will follow soon. ##

## Still need to implement testing. ##
