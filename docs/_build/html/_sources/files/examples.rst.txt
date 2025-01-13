Here are some examples what you can do with the library:
--------------------------------------------------------

These examples are just using built in functionality.

Example 1:
~~~~~~~~~~

.. figure:: ../docs/images/snom_example_1_comparison.png
   :scale: 70 %
   :alt: Image of the first example code discussed in the how to section.

   This is a comparison of the data before and after modifications, this is possible because each created image is saved in a file.
   In the second image modifications such as simple 3 point height leveling, cutting the data, scaling and adding a gaussian filter.

Example 2:
~~~~~~~~~~

.. figure:: ../docs/images/snom_example_2_comparison.png
   :scale: 100 %
   :alt: Image of the second example code discussed in the how to section.

   This is the data before and after some modifications such as adding a scalebar, a nonlinear correction to height, amplitude and phase,
   and also a phase shift to make slight phase changes more obious.



Example 3:
~~~~~~~~~~

.. figure:: ../docs/images/snom_example_3_comparison.png
   :scale: 80 %
   :alt: Image of the third example code discussed in the how to section.

   This is the data before and after the synccorrection. The synccorrection gets rid of the linear phase drift caused by the movement of the lower parabola.

Example 4:
~~~~~~~~~~

.. figure:: ../docs/images/approachcurve_example_1.png
   :scale: 80 %
   :alt: Image of the fouth example code discussed in the how to section.

   This shows basic plotting of approach curves. The data is loaded and the minimum is set to zero. The data is then displayed in a plot.

Example 5:
~~~~~~~~~~

.. figure:: ../docs/images/3dscan_example_1.png
   :scale: 70 %
   :alt: Image of the fouth example code discussed in the how to section.

   This shows basic plotting of 3D scans. The data is loaded, cutplanes are created and the minimum is set to zero. A single cutplane is then displayed in a plot.

Example 6:
~~~~~~~~~~

.. figure:: ../docs/images/3dscan_example_1_shifted.png
   :scale: 70 %
   :alt: Image of the fouth example code discussed in the how to section.

   This shows basic plotting of 3D scans. The data is loaded, cutplanes are created and the minimum is set to zero. The cutplanes are then autoaligned,
   such that the start point of each individual approach curve is identical in z, and then displayed in a plot. This leads
   to a much better physical representation of the data. In this case a waveguide was in the center of the scan.


