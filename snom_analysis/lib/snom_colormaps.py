##############################################################################
# Copyright (C) 2020-2025 Hans-Joachim Schill

# This file is part of snom_analysis.

# snom_analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# snom_analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with snom_analysis.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
# import matplotlib as mpl
# to use perceptual colormaps:
import colorcet as cc
"""
SNOM_height = 'gray'

cmap_snom_amplitude = {'red':   ((0.0, 0.0, 0.0),
           (0.33, 1.0, 1.0),
           (0.66, 1.0, 1.0),
           (1.0, 1.0, 1.0)),

  'green': ((0.0, 0.0, 0.0),
           (0.33, 0.0, 0.0),
           (0.66, 1.0, 1.0),
           (1.0, 1.0, 1.0)),

  'blue':  ((0.0, 0.0, 0.0),
           (0.33, 0.0, 0.1),
           (0.66, 0.0, 0.0),
           (1.0, 0.0, 1.0)),

  'alpha': ((0.0, 1.0, 1.0),
            (0.33, 1.0, 1.0),
            (0.66, 1.0, 1.0),
            (1.0, 0.0, 0.0)),
 }
SNOM_amplitude = LinearSegmentedColormap('SNOM_amplitude', cmap_snom_amplitude)
plt.register_cmap(cmap=SNOM_amplitude)

cmap_snom_phase = {'red':   ((0.0, 0.0, 0.0),
           (0.33, 0.0, 0.0),
           (0.66, 1.0, 1.0),
           (1.0, 1.0, 0.0)),

  'green': ((0.0, 0.0, 0.0),
           (0.33, 0.0, 0.0),
           (0.66, 1.0, 1.0),
           (1.0, 0.0, 0.0)),

  'blue':  ((0.0, 0.0, 0.0),
           (0.33, 1.0, 1.0),
           (0.66, 1.0, 1.0),
           (1.0, 0.0, 0.0)),
 }
SNOM_phase = LinearSegmentedColormap('SNOM_phase', cmap_snom_phase)
plt.register_cmap(cmap=SNOM_phase)

cmap_snom_realpart = {
  'red':   ((0.0, 0.0, 0.0),
           (0.5, 1.0, 1.0),
           (1.0, 1.0, 0.0)),

  'green': ((0.0, 0.0, 0.0),
           (0.5, 1.0, 1.0),
           (1.0, 0.0, 0.0)),

  'blue':  ((0.0, 1.0, 1.0),
           (0.5, 1.0, 1.0),
           (1.0, 0.0, 0.0)),
}
SNOM_realpart = LinearSegmentedColormap('SNOM_realpart', cmap_snom_realpart)
plt.register_cmap(cmap=SNOM_realpart)
"""

# replace old maps with perceptual colormaps
SNOM_amplitude = cc.cm.CET_L3
# SNOM_phase = cc.cm.CET_C3s
SNOM_phase = cc.cm.CET_C3
# SNOM_realpart = cc.cm.CET_D1A 
SNOM_realpart = cc.cm.CET_D1
# SNOM_realpart = cc.cm.CET_D9
SNOM_height = cc.cm.CET_L2

# all implemented colormaps
all_colormaps = {
    "<SNOM_amplitude>": SNOM_amplitude,
    "<SNOM_height>": SNOM_height,
    "<SNOM_phase>": SNOM_phase,
    "<SNOM_realpart>": SNOM_realpart
}
