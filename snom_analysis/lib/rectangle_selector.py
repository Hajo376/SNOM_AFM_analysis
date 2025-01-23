from matplotlib.widgets import RectangleSelector, Button
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .snom_colormaps import SNOM_amplitude, SNOM_phase

# This is an adaptation of the example code provided by the matplotlib documentation:
# https://matplotlib.org/3.1.3/gallery/widgets/rectangle_selector.html


def select_rectangle(data, channel):
    selector = Rectangle_Selector(data, channel)
    selection = selector.selection
    return selection

class Rectangle_Selector():
    def __init__(self, data, channel):
        self.data = data
        self.channel = channel
        self.selection = None
        self.create_plot()

    def create_plot(self):
        cmap = 'gray'
        if ('Z' in self.channel) or ('MT' in self.channel):
            cmap = 'gray'
        elif ('P' or 'arg') in self.channel:
            cmap = SNOM_phase
        elif ('A' or 'abs') in self.channel:
            cmap = SNOM_amplitude
        else:
            print('Unknown channel, could not find the proper colormap!')
        self.fig, axis = plt.subplots()
        plot = plt.pcolormesh(self.data, cmap=cmap)
        axis.invert_yaxis()
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(plot, cax=cax)
        cbar.ax.get_yaxis().labelpad = 15
        label = 'data'
        title = 'Select an area via drag and drop'
        cbar.ax.set_ylabel(label, rotation=270)
        axis.set_title(title)
        axis.axis('scaled')

        def line_select_callback(eclick, erelease):
            #eclick and erelease are the press and release events
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata
            self.selection = [[round(x1), round(y1)], [round(x2), round(y2)]]
        
        def toggle_selector(event):
            print(' Key pressed.')
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                print(' RectangleSelector deactivated.')
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                print(' RectangleSelector activated.')
                toggle_selector.RS.set_active(True)
        # drawtype is 'box' or 'line' or 'none'
        toggle_selector.RS = RectangleSelector(axis, line_select_callback,
                                            useblit=True,
                                            button=[1, 3],  # don't use middle button
                                            minspanx=5, minspany=5,
                                            spancoords='pixels',
                                            interactive=True)
        self.cid = self.fig.canvas.mpl_connect('key_press_event', toggle_selector)

        accept = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(accept, 'Accept')
        button.on_clicked(self.accept)
        plt.show()

    def accept(self, value):
        # print('value:', value)
        self.fig.canvas.mpl_disconnect(self.cid)
        plt.close()
    