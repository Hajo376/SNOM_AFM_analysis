A package to load, manipulate and visualize SNOM and AFM data.

Overview
--------

The package contains several classes, one for each implemented measurement type (so far: SNOM/AFM images, approach/deproach curves and 3D scans (2D approach curves)). 
These classes need the path to the measurement folder as well as the channels as input. 
The classes will then load the data of all specified channels as well as the measurement parameters and the header information of the measurement files. 
The data can then be manipulated and plotted. Each manipulation changes the data in the memory and also the parameter dictionaries if necessary. 
The data can then also be saved with the changes.

The package will also create a folder in the users home directory to store several files like a config file, a plot memory, a matplotlib style file and a general 
plotting parameters file. Making it easier to adjust the package to your needs.

Installation
------------

The package can be installed via ``pip``::

    pip install snom_analysis

or via ``conda`` (soon)::

    conda install -c conda-forge snom_analysis

If you install via pip all dependencies will be installed automatically. I recommend to use a virtual environment to install the package.

Documentation
-------------

The documentation can be found at https://snom-analysis.readthedocs.io/en/latest/index.html