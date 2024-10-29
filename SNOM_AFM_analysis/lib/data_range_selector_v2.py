import tkinter as tk
import numpy as np
from PIL import Image, ImageTk


# works but the non highlited areas are not visible
class ArraySelector:
    def __init__(self, root, array):
        self.root = root
        self.array = array
        self.height, self.width = array.shape
        
        # Normalize array to 0-255 for display
        self.array = ((array - np.min(array)) / (np.max(array) - np.min(array)) * 255).astype(np.uint8)
        self.image = Image.fromarray(self.array)
        self.tk_image = ImageTk.PhotoImage(self.image)
        
        self.canvas = tk.Canvas(root, width=self.tk_image.width(), height=self.tk_image.height())
        self.canvas.pack()
        
        self.canvas.create_image(0, 0, anchor=tk.NW, image=self.tk_image)
        
        self.rect = None
        self.start = None
        self.end = None
        self.dragging = False
        self.resizing = None  # 'left', 'right', or None
        self.inverted = False  # Track if the selection is inverted
        self.is_horizontal = True  # Track selection mode

        self.canvas.bind("<ButtonPress-1>", self.on_button_press)
        self.canvas.bind("<B1-Motion>", self.on_mouse_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_button_release)
        self.canvas.bind("<Motion>", self.on_mouse_move)
        
        self.select_button = tk.Button(root, text="Get Coordinates", command=self.get_coordinates)
        self.select_button.pack()
        
        self.invert_button = tk.Button(root, text="Invert Selection", command=self.invert_selection)
        self.invert_button.pack()

        self.toggle_button = tk.Button(root, text="Toggle Selection Mode", command=self.toggle_selection_mode)
        self.toggle_button.pack()
        
        self.highlighted_image = None

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
                self.end = event.x
            else:
                self.end = event.y

        # Prevent overlapping
        if self.start > self.end:
            self.start, self.end = self.end, self.start

        if self.is_horizontal:
            self.canvas.coords(self.rect, self.start, 0, self.end, self.canvas.winfo_height())
        else:
            self.canvas.coords(self.rect, 0, self.start, self.canvas.winfo_width(), self.end)

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
        # highlighted_img = Image.fromarray(highlighted_img).convert('RGBA')


        overlay = Image.new('RGBA', self.image.size, (255, 0, 0, 0))  # Red with 30% opacity
        combined = Image.composite(self.image.convert('RGBA'), overlay, mask)
        # convert to rgb
        # combined = combined.convert('RGB')
        # now update the image
        self.tk_image = ImageTk.PhotoImage(combined)
        self.canvas.create_image(0, 0, anchor=tk.NW, image=self.tk_image)




        '''# we dont need to delete the original image
        if self.highlighted_image:
            self.canvas.delete(self.highlighted_image)

        # Create an image for the highlighting based on inverted state
        highlighted_img = np.zeros_like(self.array)
        # highlighted_img = self.array.copy()
        if self.inverted:
            highlighted_img[:, :] = self.array  # Copy original array
            if self.is_horizontal:
                highlighted_img[:, self.start:self.end] = 0  # Set selected area to black
            else:
                highlighted_img[self.start:self.end, :] = 0  # Set selected area to black
        else:
            if self.is_horizontal and self.start is not None and self.end is not None:
                highlighted_img[:, self.start:self.end] = self.array[:, self.start:self.end]  # Highlight selected area
            elif not self.is_horizontal and self.start is not None and self.end is not None:
                highlighted_img[self.start:self.end, :] = self.array[self.start:self.end, :]  # Highlight selected area

        # Create an image for the highlighted area
        # highlighted_img = Image.fromarray(highlighted_img)
        highlighted_img = Image.fromarray(highlighted_img)
        overlay = Image.new('RGBA', highlighted_img.size, (255, 0, 0, 77))  # Red with 30% opacity
        combined = Image.alpha_composite(highlighted_img.convert('RGBA'), overlay)

        self.tk_highlighted = ImageTk.PhotoImage(combined)
        
        self.highlighted_image = self.canvas.create_image(0, 0, anchor=tk.NW, image=self.tk_highlighted)'''

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
        if self.start is not None and self.end is not None:
            if self.is_horizontal:
                print(f"Selected horizontal coordinates: {self.start}, {self.end}")
            else:
                print(f"Selected vertical coordinates: {self.start}, {self.end}")

if __name__ == "__main__":
    root = tk.Tk()
    
    # Create a random 2D numpy array to simulate image data
    array_data = np.random.rand(200, 400)  # Example size: 200 rows by 400 columns
    app = ArraySelector(root, array_data)
    
    root.mainloop()
