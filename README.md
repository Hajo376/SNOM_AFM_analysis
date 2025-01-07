# This package is a toolset to load, view and manipulate SNOM or AFM data. #
Currently SNOM/AFM data and approach curves are supported, although the focus so far was on the standard SNOM measurements created by a standard Neaspec (Attocube) s-SNOM.
The specific class will load one measurement and allows to handle multiple channels at once.
Such a measurement instance will contain the raw measurement data as numpy arrays, the channel names, and also a dictionary containing information about the measurement and one for each channel.
The dictionaries are created automatically from the measurement parameter '.txt' file and the individual header section in each '.gsf' file.
File format support is somewhat limited, but by extending the 'config.ini' more file formats can be added without the need to edit the source code.
The config file has one section for each filetype including parameters such as the channel names and identifiers.
Once the data has been loaded it can be modified using functions, it can also be displayed and saved with the made modifications.
The package also creates a folder in the users 'home' directory to store the 'config.ini' and also a file containing the plot memory.
Each plotting command generates a plot and that plot is automatically appended to the plot memory file, allowing to view multiple measurements besides each other.
This memory can of course be deleted at will, which is by default every time a new measurement is opened.

## Documentation will follow soon. ##
