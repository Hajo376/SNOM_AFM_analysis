How to use:
-----------

First of all install the package either from PyPi or by using a wheel from the dist folder.
You can also clone the repository and use it as is or create your own wheel by running pyhton setup.py bdist_wheel.
I would always recommend to use a virtual environment to install the package and use the wheel.

Then try out the example script which should be somewhere in the repository. (#todo)
Anyways, here is a short example of how to use the package:

Simple usage example using the SnomMeasurement class:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a very simple example using just some base functions without any extra parameters.

.. include:: ../docs/files/code_example_snommeasurement_1.rst

More advanced usage example using the SnomMeasurement class:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is an example using the tkinter filedialog and some more advanced fuctions.

.. include:: ../docs/files/code_example_snommeasurement_2.rst

Simple usage example using the ApproachCurve class:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a very simple example using just some base functions without any extra parameters.

.. include:: ../docs/files/code_example_approachcurve_1.rst

Simple usage example using the Scan3D class:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a very simple example using just some base functions to display x-z-cutplanes from the 3D scans.

.. include:: ../docs/files/code_example_scan3d_1.rst

Simple usage example using the Scan3D class and averaging the data:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is an example using just some base functions and using the multiple y lines to average each x-z-cutplane.

.. include:: ../docs/files/code_example_scan3d_2.rst