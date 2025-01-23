import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from .snom_colormaps import *



def get_height_treshold(preview_data):
    popup = Phase_Slider(preview_data)
    return popup.threshold

class Phase_Slider():
    def __init__(self, preview_data):
        self.data = preview_data
        self.threshold = None
        
        # create the plot
        self.fig, axis = plt.subplots()
        self.plot = plt.pcolormesh(self.data, cmap=SNOM_height)
        axis.invert_yaxis()
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(self.plot, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        label = 'Height [nm]'
        title = 'Select the height threshold'
        cbar.ax.set_ylabel(label, rotation=270)
        # disable ticks
        axis.yaxis.set_ticks([])
        axis.xaxis.set_ticks([])
        axis.set_title(title)
        axis.axis('scaled')


        # adjust the main plot to make room for the sliders
        plt.subplots_adjust(left=0.25, bottom=0.25)
        # Make a horizontal slider to control the frequency.
        ax_threshold = plt.axes([0.25, 0.1, 0.65, 0.03])
        self.threshold_slider = Slider(
            ax=ax_threshold,
            label='Threshold',
            valmin=0.0,
            valmax=1,
            valinit=0,
        )
        # register the update function with each slider
        self.threshold_slider.on_changed(self.update)

        # Create a `matplotlib.widgets.Button` to accept the current value and close the window.
        accept = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(accept, 'Accept', hovercolor='0.975')
        button.on_clicked(self.accept)
        plt.show()
    
    def accept(self, event):
        #close the plot, the recent phase shift is saved in self.phase_shift and can be accessed from outside
        plt.close()        
    
    def update(self, val):
        self.threshold = val 
        masked_data = self.apply_threshold()
        self.plot.set_array(masked_data)
        self.fig.canvas.draw_idle()
    
    def apply_threshold(self):
        height_flattened = self.data.flatten()
        height_threshold = self.threshold*(max(height_flattened)-min(height_flattened))+min(height_flattened)

        # create an array containing 0 and 1 depending on wether the height value is below or above threshold
        # mask = np.copy(self.data)
        mask = np.zeros_like(self.data)
        yres = len(mask)
        xres = len(mask[0])
        for y in range(yres):
            for x in range(xres):
                if self.data[y][x] >= height_threshold:
                    mask[y][x] = 1
        return np.multiply(mask, self.data)
