'''Some definitions for the snom_analysis package. Mainly enums for internal referencing.'''

from enum import Enum, auto

# keep this for internal referencing
class Definitions(Enum):
    """This class keeps track of the implemented definitions."""
    vertical = auto()
    horizontal = auto()

class MeasurementTypes(Enum):
    AFM = auto()
    SNOM = auto()
    APPROACHCURVE = auto()
    SCAN3D = auto()
    SPECTRUM = auto()
    NONE = auto()

class MeasurementTags(Enum):
    """This class keeps track of the implemented measurement tags. 
    These are tags wich are measurement specific.
    Only tags which are listed here can be used.
    """
    SCAN = auto()   # scan type, afm, snom, approach curve, 2d/3d, PsHet...
    PROJECT = auto()   
    DESCRIPTION = auto()   
    DATE = auto()   
    SCANNERCENTERPOSITION = auto()   
    ROTATION = auto()   
    SCANAREA = auto()  
    PIXELAREA = auto()   
    AVERAGING = auto()   
    INTEGRATIONTIME = auto()   
    LASERSOURCE = auto()   
    DETECTOR = auto()   
    TARGETWAVELENGTH = auto()    
    DEMODULATIONMODE = auto()   
    TIPFREQUENCY = auto()   
    TIPAMPLITUTDE = auto()   
    TAPPINGAMPLITUDE = auto()   
    MODULATIONFREQUENCY = auto()   
    MODULATIONAMPLITUDE = auto()   
    MODULATIONOFFSET = auto()   
    SETPOINT = auto()   
    REGULATOR = auto()   
    TIPPOTENTIAL = auto()   
    M1ASCALING = auto()  
    QFACTOR = auto()   
    VERSION = auto()   

class ChannelTags(Enum):
    """This class keeps track of the implemented channel tags.
    These are tags which are channel specific.
    So multiple channels might have varying channel tag values but
    they will share the same measurement tag values. However, the 
    current channel tag values are always to prefer over the same measurement tag value.
    As the channel tag values change when the channel is manipulated.
    """
    PIXELAREA = auto()  
    YINCOMPLETE = auto() 
    SCANNERCENTERPOSITION = auto()   
    ROTATION = auto()   
    SCANAREA = auto()  
    XYUNIT = auto()
    ZUNIT = auto()
    WAVENUMBERSCALING = auto()
    # additional tags
    PIXELSCALING = auto()
    # additional for aachen files (.dump)
    # INTEGRATIONTIME = auto()
    # TIPFREQUENCY = auto()
    # MODULATIONFREQUENCY = auto()
    # TAPPINGAMPLITUDE = auto()
    # MODULATIONOFFSET = auto()
    # SETPOINT = auto()
    
class PlotDefinitions:
    """This class contains all the definitions for the plotting parameters.
    """
    hide_ticks = True
    figsizex = 10
    figsizey = 5
    show_titles = True
    tight_layout = True
    colorbar_width = 2 # in percent of the fig width, standard is 2
    hspace = 0.4 #standard is 0.4
    # Define Plot font sizes
    font_size_default = 8
    font_size_axes_title = 12
    font_size_axes_label = 10
    font_size_tick_labels = 8
    font_size_legend = 8
    font_size_fig_title = 12
    #definitions for color bar ranges:
    # using the same range for all channels is useful for comparison
    # make all height channels have the same range?
    height_cbar_range = False
    vmin_height = None
    vmax_height = None
    # make all amplitude channels have the same range?
    amp_cbar_range = False
    vmin_amp = None#1 # to make shure that the values will be initialized with the first plotting command
    vmax_amp = None#-1
    # phase_cbar_range = True
    # plot the full 2pi range for the phase channels no matter what the actual data range is?
    full_phase_range = True # this will overwrite the cbar
    # make all phase channels have the same range?
    shared_phase_range = False # only used if full phase range is false
    vmin_phase = None
    vmax_phase = None
    real_cbar_range = True
    vlimit_real = None
    # vmin_real = None
    # vmax_real = None
    # show plot automatically? turn to false for gui programming
    show_plot = True
    autodelete_all_subplots = True # if true old subplots will be deleted on creation of new measurement
    # matplotlib style file
    use_mplstyle = False
   