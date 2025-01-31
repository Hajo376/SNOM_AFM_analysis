.. _code example:
How to use:
-----------

First of all install the package either from PyPi or by using a wheel from the dist folder.
You can also clone the repository and use it as is or create your own wheel by running pyhton setup.py bdist_wheel.
I would always recommend to use a virtual environment to install the package and use the wheel.

Then try out the example script which should be somewhere in the repository. (#todo)
You can also just use the script as a loader and do whatever you want in between.
The data is just stored as a list of np.array in the instance.all_data variable, the channel names are in a correlated list instance.channels.
Additional information is in the two dictionaries instance.measurement_tag_dict and instance.channels_tag_dict.
These are just dictionaries containing the parameters from the parameters.txt file (measurement_tag_dict)
and the individual headers from the .gsf files (channels_tag_dict). Note that the channels_tag_dict is also a list correlated with instance.channels.
This gives you a lot of freedome to implement your own functionality.

Anyways, here is a short example of how to use the package.

.. _code example 1:
Simple usage example using the SnomMeasurement class:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a very simple example using just some base functions without any extra parameters.


.. include:: ../docs/files/code_example_snommeasurement_1.rst

.. _code example 2:
More advanced usage example using the SnomMeasurement class:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is an example using the tkinter filedialog and some more advanced fuctions.

.. include:: ../docs/files/code_example_snommeasurement_2.rst

.. _code example 3:
Example showing how to do a synccorrection of transmission mode measurements:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is an example using the tkinter filedialog and some more advanced fuctions.

.. include:: ../docs/files/code_example_snommeasurement_3.rst

.. _code example 4:
Example showing how to create a gif from the realpart data, useful to visualize traveling waves:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a very simple example using just some base functions without any extra parameters.

.. include:: ../docs/files/code_example_snommeasurement_4.rst

.. _code example 5:
Simple usage example using the ApproachCurve class:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a very simple example using just some base functions without any extra parameters.

.. include:: ../docs/files/code_example_approachcurve_1.rst

.. _code example 6:
Simple usage example using the Scan3D class:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a very simple example using just some base functions to display x-z-cutplanes from the 3D scans.

.. include:: ../docs/files/code_example_scan3d_1.rst

.. _code example 7:
Simple usage example using the Scan3D class and averaging the data:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is an example using just some base functions and using the multiple y lines to average each x-z-cutplane.

.. include:: ../docs/files/code_example_scan3d_2.rst




:ref:`code example`