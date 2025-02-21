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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from .snom_colormaps import *



def get_phase_offset(preview_data):
    phase_shift = PhaseSlider(preview_data)
    return phase_shift.phase_shift

class PhaseSlider():
    def __init__(self, preview_data):
        self.data = preview_data
        initial_shift = 0 # initial phase shift
        self.previous_shift = 0
        self.phase_shift = 0
        # create the plot
        self.fig, axis = plt.subplots()
        self.plot = plt.pcolormesh(self.data, cmap=SNOM_phase, vmin=0, vmax=np.pi*2)
        axis.invert_yaxis()
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(self.plot, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        label = 'Phase'
        title = 'Shift the phase'
        cbar.ax.set_ylabel(label, rotation=270)
        axis.set_title(title)
        axis.axis('scaled')


        # adjust the main plot to make room for the sliders
        plt.subplots_adjust(left=0.25, bottom=0.25)
        # Make a horizontal slider to control the frequency.
        ax_phaseshift = plt.axes([0.25, 0.1, 0.65, 0.03])
        self.phase_slider = Slider(
            ax=ax_phaseshift,
            label='Frequency [Hz]',
            valmin=0.0,
            valmax=np.pi*2,
            valinit=initial_shift,
        )
        # register the update function with each slider
        self.phase_slider.on_changed(self.update)

        # Create a `matplotlib.widgets.Button` to accept the current value and close the window.
        accept = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(accept, 'Accept', hovercolor='0.975')
        button.on_clicked(self.accept)
        plt.show()
    
    def accept(self, event):
        #close the plot, the recent phase shift is saved in self.phase_shift and can be accessed from outside
        plt.close()        
    
    def update(self, val):
        self.phase_shift = val 
        self.shift_phase()
        self.plot.set_array(self.data)
        self.fig.canvas.draw_idle()
    
    def shift_phase(self):
        yres = len(self.data)
        xres = len(self.data[0])
        for y in range(yres):
            for x in range(xres):
                self.data[y][x] = (self.data[y][x] + (self.phase_shift - self.previous_shift)) % (2*np.pi)
        self.previous_shift = self.phase_shift
