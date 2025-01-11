import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
from matplotlib.backend_bases import MouseButton
import matplotlib.lines as mlines
# import tkinter as tk
from .snom_colormaps import SNOM_amplitude, SNOM_phase, SNOM_height

import skimage as ski

def select_profile(data, channel):
    # root = tk.Tk()
    # selector = ArraySelector(root, data, channel)
    selector = CutlineSelector(data, channel)
    plt.show()
    # root.mainloop()

    # create profile 
    profile = ski.measure.profile_line(data.T, selector.start, selector.end, linewidth=selector.width) # somehow x and y are switched, therefore transpose the array
    # also return the physical length of the profile and the integration width
    return profile, selector.start, selector.end, selector.width

class CutlineSelector:
    def __init__(self, img_array, channel):
        self.img_array = img_array
        self.channel = channel
        self.fig, self.ax = plt.subplots()
        # determine the colormap
        if ('Z' in self.channel) or ('MT' in self.channel):
            self.cmap = SNOM_height
        elif ('P' or 'arg') in self.channel:
            self.cmap = SNOM_phase
        elif ('A' or 'abs') in self.channel:
            self.cmap = SNOM_amplitude
        elif ('H' or 'height') in self.channel:
            self.cmap = SNOM_height
        else:
            self.cmap = 'gray'
            print('Unknown channel, could not find the proper colormap!')
        self.ax.imshow(img_array, cmap='gray')
        
        # Calculate initial width: 1/20th of the smallest dimension but at least 1
        min_dimension = min(img_array.shape)
        self.width = max(1, min_dimension // 20)

        self.start = None
        self.end = None
        self.dragging_start = False
        self.dragging_end = False
        self.dragging = False
        self.shift_pressed = False
        self.tolerance = 10  # Tolerance in pixels for grabbing endpoints

        # Matplotlib objects for the cutline and boundary lines
        self.main_line, = self.ax.plot([], [], 'r-', linewidth=1.5)
        self.perp_lines = []

        # Connect Matplotlib event handlers
        self.cid_click = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.cid_release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cid_key_press = self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.cid_key_release = self.fig.canvas.mpl_connect('key_release_event', self.on_key_release)

        # Adding buttons and sliders
        self.add_widgets()

    def add_widgets(self):
        # Adjust button placement
        button_ax = plt.axes([0.82, 0.025, 0.15, 0.04])
        self.button = Button(button_ax, 'Extract Cutline')
        self.button.on_clicked(self.extract_cutline)

        # Slider for adjusting width
        slider_ax = plt.axes([0.15, 0.025, 0.5, 0.04], facecolor='lightgoldenrodyellow')
        self.slider = Slider(slider_ax, 'Width', 1, 100, valinit=self.width, valstep=1)
        self.slider.on_changed(self.update_width)

    def update_width(self, val):
        self.width = int(val)
        self.update_cutline()

    def on_key_press(self, event):
        if event.key == 'shift':
            self.shift_pressed = True

    def on_key_release(self, event):
        if event.key == 'shift':
            self.shift_pressed = False

    def on_click(self, event):
        if event.button == MouseButton.LEFT and event.inaxes == self.ax:
            x, y = event.xdata, event.ydata

            # Check if clicking near the start or end point to enable dragging
            if self.start and self.is_near_point((x, y), self.start):
                self.dragging_start = True
            elif self.end and self.is_near_point((x, y), self.end):
                self.dragging_end = True
            else:
                # Start a new cutline
                self.start = (x, y)
                self.end = None
                self.dragging = True
                self.clear_cutline()

    def on_release(self, event):
        self.dragging = False
        self.dragging_start = False
        self.dragging_end = False

    def on_motion(self, event):
        if event.inaxes != self.ax:
            return

        x, y = event.xdata, event.ydata

        # Change cursor when near endpoints
        if self.start and self.is_near_point((x, y), self.start):
            self.fig.canvas.set_cursor(3)  # Resize cursor
        elif self.end and self.is_near_point((x, y), self.end):
            self.fig.canvas.set_cursor(3)  # Resize cursor
        else:
            self.fig.canvas.set_cursor(1)  # Default cursor

        if self.dragging:
            # Draw the cutline in real-time
            if self.start:
                if self.shift_pressed:
                    x0, y0 = self.start
                    # Snap to horizontal or vertical if shift is pressed
                    if abs(x - x0) > abs(y - y0):
                        y = y0
                    else:
                        x = x0
                self.end = self.clip_point_to_image(x, y)
        elif self.dragging_start:
            # Adjust the start point while keeping orientation
            self.start = self.adjust_endpoint(self.end, (x, y))
        elif self.dragging_end:
            # Adjust the end point while keeping orientation
            self.end = self.adjust_endpoint(self.start, (x, y))

        self.update_cutline()

    def is_near_point(self, p1, p2):
        """Check if point p1 is near point p2 within a certain tolerance."""
        return np.hypot(p1[0] - p2[0], p1[1] - p2[1]) < self.tolerance

    def adjust_endpoint(self, fixed_point, new_point):
        """Adjust endpoint along the fixed orientation."""
        fx, fy = fixed_point
        dx = new_point[0] - fx
        dy = new_point[1] - fy
        direction = np.array([dx, dy])
        unit_direction = direction / np.hypot(*direction)
        length = np.hypot(dx, dy)
        return (fx + unit_direction[0] * length, fy + unit_direction[1] * length)

    def clip_point_to_image(self, x, y):
        """Ensure point (x, y) is within the image bounds."""
        x = np.clip(x, 0, self.img_array.shape[1] - 1)
        y = np.clip(y, 0, self.img_array.shape[0] - 1)
        return x, y

    def clear_cutline(self):
        self.main_line.set_data([], [])
        for pline in self.perp_lines:
            pline.remove()
        self.perp_lines = []
        self.fig.canvas.draw()

    def update_cutline(self):
        if self.start is None or self.end is None:
            return

        # Update the main cutline
        x0, y0 = self.start
        x1, y1 = self.end
        self.main_line.set_data([x0, x1], [y0, y1])

        # Draw perpendicular lines to indicate width
        self.draw_perpendicular_lines(x0, y0, x1, y1)
        self.fig.canvas.draw()

    def draw_perpendicular_lines(self, x0, y0, x1, y1):
        for pline in self.perp_lines:
            pline.remove()
        self.perp_lines = []

        dx, dy = x1 - x0, y1 - y0
        length = np.hypot(dx, dy)
        if length == 0:
            return
        # self.width = length # to be able to access the width of the profile line from outside the class

        # Normalize direction
        dx /= length
        dy /= length

        offset_x = -dy * self.width / 2
        offset_y = dx * self.width / 2

        perp_line1 = mlines.Line2D([x0 + offset_x, x0 - offset_x], [y0 + offset_y, y0 - offset_y], color='r')
        perp_line2 = mlines.Line2D([x1 + offset_x, x1 - offset_x], [y1 + offset_y, y1 - offset_y], color='r')

        self.ax.add_line(perp_line1)
        self.ax.add_line(perp_line2)
        self.perp_lines = [perp_line1, perp_line2]

    def extract_cutline(self, event):
        if self.start is None or self.end is None:
            print("No cutline selected.")
            return
        print(f"Cutline from {self.start} to {self.end} extracted.")
        # close the figure
        plt.close(self.fig)
