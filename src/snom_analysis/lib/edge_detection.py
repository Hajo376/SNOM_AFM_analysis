import matplotlib.pyplot as plt
# import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from skimage import feature
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable # for colorbar
import scipy.ndimage as ndimage


# initial values for the sliders
initial_threshold_low = 0.2
initial_threshold_high = 0.8
initial_sigma = 4

class EdgeDetection():
    """This class creates a GUI to adjust the parameters for the edge detection algorithm.
    The user can adjust the sigma value, the low and high threshold for the edge detection.
    The user can accept the values or exit the GUI.
    Uppon accepting the values the GUI is closed and the values are saved to the class instance variables.
    """
    def __init__(self, data):
        self.data = data
        self.iterations = None
        # create the gui
        self.create_gui()
    
    def create_gui(self):
        self.fig, axis = plt.subplots()
        axis.imshow(self.data, cmap='gray')
        self.edges = self.calculate_edges(initial_sigma, initial_threshold_low, initial_threshold_high)
        masked_data = np.ma.masked_where(self.edges != 1, self.edges)
        self.plot = axis.imshow(masked_data, interpolation='none', cmap='viridis', vmin=0, vmax=1)
        axis.invert_yaxis()
        # divider = make_axes_locatable(axis)
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # cbar = plt.colorbar(self.plot, cax=cax)
        # cbar.ax.get_yaxis().labelpad = 15
        # label = 'Height'
        title = 'Adjust the sliders to optimize the edge detection. Press "Accept" to save the settings.'
        # cbar.ax.set_ylabel(label, rotation=270)
        axis.set_title(title)
        axis.axis('scaled')

        # adjust the main plot to make room for the sliders
        plt.subplots_adjust(left=0.25, bottom=0.25)
        # Make a horizontal slider to control the frequency.
        ax_threshold_low = plt.axes([0.25, 0.1, 0.65, 0.03])
        self.slider_threshold_low = Slider(
            ax=ax_threshold_low,
            label='Threshold low',
            valmin=0.0,
            valmax=1,
            valinit=initial_threshold_low,
        )
        # register the update function with each slider
        self.slider_threshold_low.on_changed(self.update)

        ax_threshold_high = plt.axes([0.25, 0.15, 0.65, 0.03])
        self.slider_threshold_high = Slider(
            ax=ax_threshold_high,
            label='Threshold high',
            valmin=0.0,
            valmax=1,
            valinit=initial_threshold_high,
        )
        # register the update function with each slider
        self.slider_threshold_high.on_changed(self.update)

        ax_sigma = plt.axes([0.25, 0.05, 0.4, 0.03])
        self.slider_sigma = Slider(
            ax=ax_sigma,
            label='Sigma',
            valmin=0.0,
            valmax=10,
            valinit=initial_sigma,
        )
        # register the update function with each slider
        self.slider_sigma.on_changed(self.update)

        # Create a `matplotlib.widgets.Button` to accept the current value and close the window.
        accept = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(accept, 'Accept', hovercolor='0.975')
        button.on_clicked(self.accept)

        exit = plt.axes([0.7, 0.025, 0.1, 0.04])
        button_exit = Button(exit, 'Exit', hovercolor='0.975')
        button_exit.on_clicked(self.exit)

        # show the plot in full screen, important for high resolution data
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()

        plt.show()
        
    def calculate_edges(self, sigma, threshold_low, threshold_high):
        edges = feature.canny(self.data, sigma=sigma, use_quantiles=True, low_threshold=threshold_low, high_threshold=threshold_high)
        return edges

    def accept(self, event):
        #close the plot, the recent values are saved to class instance variables and can be accessed from outside
        plt.close()   
        self.threshold_high = self.slider_threshold_high.val
        self.threshold_low = self.slider_threshold_low.val
        self.sigma = self.slider_sigma.val   

    def exit(self, event):
        #close the plot, the recent values are saved to class instance variables and can be accessed from outside
        plt.close()
        # set the values to the initial values
        self.threshold_high = initial_threshold_high
        self.threshold_low = initial_threshold_low
        self.sigma = initial_sigma
    
    def update(self, val):
        self.phase_shift = val 
        # calculate the new data
        self.edges = self.calculate_edges(self.slider_sigma.val, self.slider_threshold_low.val, self.slider_threshold_high.val)
        # self.edges = feature.canny(adjusted_height_data, sigma=self.slider_sigma.val, use_quantiles=True, low_threshold=self.slider_threshold_low.val, high_threshold=self.slider_threshold_high.val)
        
        # if resolution of data is too high compared to screen resolution then dilate the edges
        # the iterations have to be calculated only once
        if self.iterations is None:
            # get resolution of data
            yres, xres = self.data.shape
            # get resolution of screen
            xscreen, yscreen = plt.gcf().get_size_inches() * plt.gcf().dpi
            # if the resolution of the data is higher than the screen resolution then approximate the iterations for the dilation
            if xres > xscreen or yres > yscreen:
                xres_factor = int(xres/xscreen)
                yres_factor = int(yres/yscreen)
                self.iterations = max(xres_factor, yres_factor)
            else:
                self.iterations = 0
        # if the iterations are not 0 then dilate the edges
        if self.iterations != 0:
            self.edges = ndimage.binary_dilation(self.edges, iterations=self.iterations)

        # update the plot
        # mask the edges where the edges are not 1 to make onl the edges visible in the overlay with the height data
        masked_data = np.ma.masked_where(self.edges != 1, self.edges)
        self.plot.set_data(masked_data)
        # self.fig.canvas.draw_idle() # not shure what the difference is between draw_idle and draw
        plt.draw()