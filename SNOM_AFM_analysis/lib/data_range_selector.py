import tkinter as tk
import numpy as np
from PIL import Image, ImageTk
import cv2

from SNOM_AFM_analysis.lib.snom_colormaps import SNOM_amplitude, SNOM_phase, SNOM_height


def select_data_range(data, channel):
    root = tk.Tk()
    selector = ArraySelector(root, data, channel)
    root.mainloop()

    # selection = selector.selection
    start = selector.start
    end = selector.end
    is_horizontal = selector.is_horizontal
    inverted = selector.inverted
    return start, end, is_horizontal, inverted

class ArraySelector:
    def __init__(self, root, data, channel):
        self.root = root
        self.array = data
        self.original_array = self.array.copy()
        self.channel = channel
        self.height, self.width = data.shape
        self.original_height, self.original_width = data.shape

        # additional parameters
        self.highlighted_image = None
        self.rect = None
        self.start = None
        self.end = None
        self.dragging = False
        self.resizing = None  # 'left', 'right', or None
        self.inverted = False  # Track if the selection is inverted
        self.is_horizontal = True  # Track selection mode
        
        # scale the data to a size usable for the canvas, min should be 300x300 pixels
        # if both axes are smaller than 300 pixels then scale the data such that the bigger axis has at least 300 pixels
        min_size = 400
        max_size = 800
        self.scaling_factor = 1
        # first check if the data is already larger than the minimum size in at least one axis
        if self.width > min_size or self.height > min_size:
            pass
        # if both axes are smaller than the minimum size then scale the data such that the bigger axis has at least min_size pixels
        elif self.width < min_size and self.height < min_size:
            if self.width > self.height:
                self.scaling_factor = min_size/self.width
            else:
                self.scaling_factor = min_size/self.height
        elif self.width < min_size:
            self.scaling_factor = min_size/self.width
        elif self.height < min_size:
            self.scaling_factor = min_size/self.height
        # check if one of the axes is larger than the maximum size
        if self.width > max_size or self.height > max_size:
            if self.width > self.height:
                self.scaling_factor = max_size/self.width
            else:
                self.scaling_factor = max_size/self.height
        elif self.width > max_size:
            self.scaling_factor = max_size/self.width
        elif self.height > max_size:
            self.scaling_factor = max_size/self.height
        self.update_scaling_factor()

        # change colormap depending on the channel
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

        canvas_width = self.width*self.scaling_factor
        canvas_height = self.height*self.scaling_factor

        self.canvas = tk.Canvas(self.root, width=canvas_width, height=canvas_height)

        self.canvas.grid(row=0, column=0, sticky='nsew')

        self.fill_canvas()

        self.canvas.bind("<ButtonPress-1>", self.on_button_press)
        self.canvas.bind("<B1-Motion>", self.on_mouse_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_button_release)
        self.canvas.bind("<Motion>", self.on_mouse_move)
        self.canvas.bind("<Configure>", self.on_windowsize_changed)
        
        self.toggle_button = tk.Button(root, text="Toggle Selection Mode", command=self.toggle_selection_mode)
        self.toggle_button.grid(row=1, column=0)
        
        self.invert_button = tk.Button(root, text="Invert Selection", command=self.invert_selection)
        self.invert_button.grid(row=2, column=0)

        self.select_button = tk.Button(root, text="Get Coordinates", command=self.get_coordinates)
        self.select_button.grid(row=3, column=0)
        
        # Center the window
        root.eval(f'tk::PlaceWindow {str(self.root)} center')

        # configure canvas to scale with window
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

    def update_scaling_factor(self):
        if self.scaling_factor != 1:
            if self.scaling_factor > 1:
                self.scaling_factor = round(self.scaling_factor) # allow for scaling up to an arbitrary integer of the original size
            else:
                self.scaling_factor = 0.5 # allow for scaling down to half the size

    def on_button_press(self, event):
        if self.resizing:
            self.dragging = True
            return
        
        if self.is_horizontal:
            self.start = event.x
            self.end = event.x
            if self.rect:
                self.canvas.delete(self.rect)
            self.rect = self.canvas.create_rectangle(self.start, 0, self.start, self.canvas.winfo_height(), outline='red', width=2)
        else:
            self.start = event.y
            self.end = event.y
            if self.rect:
                self.canvas.delete(self.rect)
            self.rect = self.canvas.create_rectangle(0, self.start, self.canvas.winfo_width(), self.start, outline='red', width=2)

    def on_mouse_drag(self, event):
        if self.dragging:
            if self.resizing == 'left':
                if self.is_horizontal:
                    self.start = max(0, event.x)
                else:
                    self.start = max(0, event.y)
            elif self.resizing == 'right':
                if self.is_horizontal:
                    self.end = min(self.width, event.x)
                else:
                    self.end = min(self.height, event.y)
        else:
            if self.is_horizontal:
                # self.end = event.x
                if event.x < 0:
                    self.end = 0
                elif event.x > self.width:
                    self.end = self.width
                else:
                    self.end = event.x
            # else:
            #     self.end = event.y
            else:
                self.end = event.y
                if event.y < 0:
                    self.start = 0
                elif event.y > self.height:
                    self.end = self.height

        # Prevent overlapping
        if self.start > self.end:
            self.start, self.end = self.end, self.start

        if self.is_horizontal:
            self.canvas.coords(self.rect, self.start, 0, self.end, self.canvas.winfo_height())
        else:
            self.canvas.coords(self.rect, 0, self.start, self.canvas.winfo_width(), self.end)

        #print start and end
        print(f"Start: {self.start}, End: {self.end}")
        self.highlight_selection()

    def on_button_release(self, event):
        self.dragging = False
        self.resizing = None

    def on_mouse_move(self, event):
        if self.rect:
            if self.is_horizontal:
                left_border = self.start
                right_border = self.end
                border_width = 5
                
                if left_border - border_width < event.x < left_border + border_width:
                    self.canvas.config(cursor="sb_h_double_arrow")
                    self.resizing = 'left'
                elif right_border - border_width < event.x < right_border + border_width:
                    self.canvas.config(cursor="sb_h_double_arrow")
                    self.resizing = 'right'
                else:
                    self.canvas.config(cursor="")
                    self.resizing = None
            else:
                top_border = self.start
                bottom_border = self.end
                border_width = 5

                if top_border - border_width < event.y < top_border + border_width:
                    self.canvas.config(cursor="sb_v_double_arrow")
                    self.resizing = 'left'
                elif bottom_border - border_width < event.y < bottom_border + border_width:
                    self.canvas.config(cursor="sb_v_double_arrow")
                    self.resizing = 'right'
                else:
                    self.canvas.config(cursor="")
                    self.resizing = None

    def fill_canvas(self):#
        self.array = cv2.resize(self.original_array, (int(self.scaling_factor*self.original_width), int(self.scaling_factor*self.original_height)),interpolation=cv2.INTER_NEAREST)
        # Normalize array to 0-255 for display
        self.array = ((self.array - np.min(self.array)) / (np.max(self.array) - np.min(self.array)) * 255).astype(np.uint8)

        self.height, self.width = self.array.shape
        # delete and redraw the image if already existing
        try:
            self.canvas.delete("all")
        except:
            pass
        self.image = Image.fromarray(np.uint8(self.cmap(self.array)*255))
        self.tk_image = ImageTk.PhotoImage(self.image)
        self.canvas.create_image(self.width/2, self.height/2, image=self.tk_image)
        self.highlight_selection()

    def on_windowsize_changed(self, event):
        # if the canvas size becomse so large that the image scaled by an integer factor is smaller than the canvas size, then we need to rescale the image by that factor
        # get the current size of the canvas
        canvas_width = self.canvas.winfo_width()
        canvas_height = self.canvas.winfo_height()
        # get the current size of the data
        height, width = self.array.shape
        # if the current scaling factor is 0.5 we need to increment it by 0.5 not 1 
        # if the current scaling factor is 1 we need to increment it by 0.5 not 1
        if self.scaling_factor == 0.5 or self.scaling_factor == 1:
            scaling_increment = 0.5
        else:
            scaling_increment = 1
        # resize the image to the canvas size if the canvas size is larger than the image times the scaling factor, where the scaling factor is an integer
        if canvas_width > self.original_width*(self.scaling_factor+scaling_increment) and canvas_height > self.original_height*(self.scaling_factor+scaling_increment):
            # increment the scaling factor
            self.scaling_factor += scaling_increment
            self.update_scaling_factor()
            self.fill_canvas()
        elif canvas_width < width or canvas_height < height:
            # decrement the scaling factor
            self.scaling_factor -= scaling_increment
            self.update_scaling_factor()
            self.fill_canvas()        

    def highlight_selection(self):
        if self.highlighted_image:
            self.canvas.delete(self.highlighted_image)

        # Create an image for the highlighting based on inverted state
        # highlighted_img = np.zeros_like(self.array)
        highlighted_img = np.ones_like(self.array).astype(np.uint8)*255
        # highlighted_img = self.array.copy()
        if self.inverted:
            highlighted_img[:, :] = 128  # Copy original array
            if self.is_horizontal:
                highlighted_img[:, self.start:self.end] = 255  # Set selected area to black
            else:
                highlighted_img[self.start:self.end, :] = 255  # Set selected area to black
        else:
            if self.is_horizontal and self.start is not None and self.end is not None:
                highlighted_img[:, self.start:self.end] = 128#self.array[:, self.start:self.end]  # Highlight selected area
            elif not self.is_horizontal and self.start is not None and self.end is not None:
                highlighted_img[self.start:self.end, :] = 128#self.array[self.start:self.end, :]  # Highlight selected area
        mask = Image.fromarray(highlighted_img)
        # create an overlay in red with 30% opacity of the highlighted area with the original image
        overlay = Image.new('RGBA', self.image.size, (255, 0, 0, 0))  # Red with 30% opacity
        combined = Image.composite(self.image.convert('RGBA'), overlay, mask)
        # convert to rgb
        # combined = combined.convert('RGB')
        # now update the image
        self.tk_image = ImageTk.PhotoImage(combined)
        self.canvas.create_image(0, 0, anchor=tk.NW, image=self.tk_image)

    def invert_selection(self):
        self.inverted = not self.inverted  # Toggle inverted state
        self.highlight_selection()  # Update highlighting

    def toggle_selection_mode(self):
        # Clear current selection
        if self.rect:
            self.canvas.delete(self.rect)
            self.rect = None
        self.start = None
        self.end = None

        self.is_horizontal = not self.is_horizontal  # Toggle selection mode
        self.highlight_selection()  # Update highlighting for the new mode

    def get_coordinates(self):
        # reduze the coordinates to the original size
        if self.scaling_factor != 1:
            self.start = int(self.start/self.scaling_factor)
            self.end = int(self.end/self.scaling_factor)-1 # '-1' due to conversion between pixel in image and array index
        if self.start is not None and self.end is not None:
            if self.is_horizontal:
                print(f"Selected horizontal coordinates: {self.start}, {self.end}")
            else:
                print(f"Selected vertical coordinates: {self.start}, {self.end}")

        # close the window
        self.root.destroy()
