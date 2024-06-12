from matplotlib.widgets import RectangleSelector, Button
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from SNOM_AFM_analysis.lib.snom_colormaps import SNOM_amplitude, SNOM_phase




def Select_Rectangle(data, channel):
    selector = Rectangle_Selector(data, channel)
    selection = selector.selection
    # print('selection: ', selection)
    return selection

class Rectangle_Selector():
    def __init__(self, data, channel):
        self.data = data
        self.channel = channel
        self.selection = None
        self.Create_Plot()

    def Create_Plot(self):
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
        print("\n      click  -->  release")

        def line_select_callback(eclick, erelease):
            'eclick and erelease are the press and release events'
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata
            self.selection = [[round(x1), round(y1)], [round(x2), round(y2)]]
            # print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
            # print(" The button you used were: %s %s" % (eclick.button, erelease.button))
        
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
        # plt.connect('key_press_event', toggle_selector)
        self.cid = self.fig.canvas.mpl_connect('key_press_event', toggle_selector)

        # def mycallback(event):
        #     if toggle_selector.RS.active:
        #         toggle_selector.RS.update()
        # plt.connect('draw_event', mycallback)
        accept = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(accept, 'Accept')
        button.on_clicked(self.accept)
        plt.show()

    def accept(self, value):
        print('value:', value)
        self.fig.canvas.mpl_disconnect(self.cid)
        plt.close()
    




'''
def Do_Stuff(data):
    fig, current_ax = plt.subplots()                 # make a new plotting range
    # N = 100000                                       # If N is large one can see
    # x = np.linspace(0.0, 10.0, N)                    # improvement by use blitting!

    # plt.plot(x, +np.sin(.2*np.pi*x), lw=3.5, c='b', alpha=.7)  # plot something
    # plt.plot(x, +np.cos(.2*np.pi*x), lw=3.5, c='r', alpha=.5)
    # plt.plot(x, -np.sin(.2*np.pi*x), lw=3.5, c='g', alpha=.3)
    plt.pcolormesh(data, cmap=SNOM_phase)
    print("\n      click  -->  release")

    # drawtype is 'box' or 'line' or 'none'
    toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                                        drawtype='box', useblit=True,
                                        button=[1, 3],  # don't use middle button
                                        minspanx=5, minspany=5,
                                        spancoords='pixels',
                                        interactive=True)
    plt.connect('key_press_event', toggle_selector)
    # # Create a `matplotlib.widgets.Button` to accept the current value and close the window.
    # accept = plt.axes([0.8, 0.025, 0.1, 0.04])
    # button = Button(accept, 'Accept', hovercolor='0.975')
    # button.on_clicked(accept)
    plt.show()

def accept(value):
    print('value:', value)
    plt.close()


def line_select_callback(eclick, erelease):
    'eclick and erelease are the press and release events'
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
    print(" The button you used were: %s %s" % (eclick.button, erelease.button))


def toggle_selector(event):
    print(' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print(' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print(' RectangleSelector activated.')
        toggle_selector.RS.set_active(True)

'''
