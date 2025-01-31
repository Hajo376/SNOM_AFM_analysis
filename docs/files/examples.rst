Here are some examples what you can do with the library:
--------------------------------------------------------

These examples are just using built in functionality.

.. _plot example 1:
Example 1:
~~~~~~~~~~

.. figure:: images/snom_example_1_comparison.png
   :width: 1000 px
   :height: 500 px
   :scale: 90 %
   :alt: Image of the first example code discussed in the how to section.

   This is a comparison of the data before and after modifications, this is possible because each created image is saved in a file.
   In the second image modifications such as simple 3 point height leveling, cutting the data, scaling and adding a gaussian filter.
   This measurement was performed in transmission mode on a thin gold pentamer structur on a glass substrate.

.. dropdown:: Code

   .. include:: ../docs/files/code_example_snommeasurement_1.rst

.. _plot example 2:
Example 2:
~~~~~~~~~~

.. figure:: images/snom_example_2_comparison.png
   :width: 601 px
   :height: 903 px
   :scale: 100 %
   :alt: Image of the second example code discussed in the how to section.

   This is the data before and after some modifications such as adding a scalebar, a nonlinear correction to height, amplitude and phase,
   and also a phase shift to make slight phase changes more obvious. This measurement was performed in reflection mode on a thin PMMA wedge on top of a gold film.
   The fringes you see are an interference of the surface plasmon polaritons excited by the tip and the reflection from the PMMA edges.
   Due to poor AFM stability at that time the amplitude and phase data drifted quite a bit over time and need correction. Since the drifts are not linear we cannot simply
   use a linear correction.

.. dropdown:: Code

   .. include:: ../docs/files/code_example_snommeasurement_2.rst

.. _plot example 3:
Example 3:
~~~~~~~~~~

.. figure:: images/snom_example_3_comparison.png
   :width: 648 px
   :height: 490 px
   :scale: 100 %
   :alt: Image of the third example code discussed in the how to section.

   This shows the phase data before and after the synccorrection. The synccorrection gets rid of the linear phase drift caused by the movement of the lower parabola.
   This measurement was performed in transmission mode on a grating milled inside of a gold film. What you see is the excitation of surface plasmon polaritons propagating
   to the left and right of the grating.

.. dropdown:: Code

   .. include:: ../docs/files/code_example_snommeasurement_3.rst

.. _plot example 4:
Example 4:
~~~~~~~~~~

.. figure:: images/snom_example_4.gif
   :width: 200 px
   :height: 200 px
   :scale: 100 %
   :alt: Image of the fouth example code discussed in the how to section.

   This shows a gif created from the realpart of the O3A channel and the O3P_corrected channel.

.. dropdown:: Code

   .. include:: ../docs/files/code_example_snommeasurement_4.rst

.. _plot example 5:
Example 5:
~~~~~~~~~~

.. figure:: images/approachcurve_example_1.png
   :width: 640 px
   :height: 480 px
   :scale: 80 %
   :alt: Image of the fouth example code discussed in the how to section.

   This shows basic plotting of approach curves. The data is loaded and the minimum is set to zero. The data is then displayed in a plot.

.. dropdown:: Code

   .. include:: ../docs/files/code_example_approachcurve_1.rst

.. _plot example 6:
Example 6:
~~~~~~~~~~

.. figure:: images/3dscan_example_1.png
   :width: 1000 px
   :height: 500 px
   :scale: 80 %
   :alt: Image of the fouth example code discussed in the how to section.

   This shows basic plotting of 3D scans. The data is loaded, cutplanes are created and the minimum is set to zero. A single cutplane is then displayed in a plot.
   The measurement was performed on a dielectric loaded surface plasmon polariton waveguide on top of a gold film. The measurement is a cut perpendicular to the waveguide.

.. dropdown:: Code

   .. include:: ../docs/files/code_example_scan3d_1.rst

.. _plot example 7:
Example 7:
~~~~~~~~~~

.. figure:: images/3dscan_example_1_shifted.png
   :width: 1000 px
   :height: 500 px
   :scale: 80 %
   :alt: Image of the fouth example code discussed in the how to section.

   This shows basic plotting of 3D scans. The data is loaded, cutplanes are created and the minimum is set to zero. The cutplanes are then autoaligned,
   such that the start point of each individual approach curve is identical in z, and then displayed in a plot. This leads
   to a much better physical representation of the data. In this case a waveguide was in the center of the scan. This image is equivalent to the previous one
   but the data is shifted to align the waveguide in the center of the image.

.. dropdown:: Code

   .. include:: ../docs/files/code_example_scan3d_2.rst
