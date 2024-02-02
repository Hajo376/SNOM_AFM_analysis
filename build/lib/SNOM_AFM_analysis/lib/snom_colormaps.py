from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

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


# all implemented colormaps
all_colormaps = {
    "<SNOM_amplitude>": SNOM_amplitude,
    "<SNOM_height>": SNOM_height,
    "<SNOM_phase>": SNOM_phase,
    "<SNOM_realpart>": SNOM_realpart
}
