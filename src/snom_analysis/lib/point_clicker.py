import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from mpl_point_clicker import clicker


class ImageClicker:
    def __init__(self, data: np.ndarray, cmap="gray", message: str=None):
        if data.ndim != 2:
            raise ValueError("Input must be a 2D NumPy array")

        self.data = data
        self.cmap = cmap
        if message is None:
            self.message = "Click points, then press Accept"
        else:
            self.message = message
        self.coords = None
        self._build_ui()

    def _build_ui(self):
        self.fig, self.ax = plt.subplots()
        plt.subplots_adjust(bottom=0.2)

        self.im = self.ax.imshow(self.data, cmap=self.cmap)
        self.ax.set_title(self.message)

        self.clicker = clicker(self.ax, ["points"], markers=["x"])

        # Accept button
        ax_btn = plt.axes([0.4, 0.05, 0.2, 0.075])
        self.btn_accept = Button(ax_btn, "Accept")
        self.btn_accept.on_clicked(self._accept)

        plt.show()

    def _accept(self, event):        
        pts = self.clicker.get_positions()["points"]
        self.coords = [[round(element[0]), round(element[1])] for element in pts]
        plt.close(self.fig)
        
