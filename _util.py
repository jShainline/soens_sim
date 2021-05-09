import numpy as np
from matplotlib import pyplot as plt
from pylab import *

def color_dictionary():

    colors = dict()    

    ## define colors
    #blues  lightest to darkest
    blueVec1 = np.array([145,184,219]); colors['blue1'] = blueVec1/256;
    blueVec2 = np.array([96,161,219]); colors['blue2'] = blueVec2/256;
    blueVec3 = np.array([24,90,149]); colors['blue3'] = blueVec3/256;
    blueVec4 = np.array([44,73,100]); colors['blue4'] = blueVec4/256;
    blueVec5 = np.array([4,44,80]); colors['blue5'] = blueVec5/256;
    #reds  lightest to darkest
    redVec1 = np.array([246,177,156]); colors['red1'] = redVec1/256;
    redVec2 = np.array([246,131,98]); colors['red2'] = redVec2/256;
    redVec3 = np.array([230,69,23]); colors['red3'] = redVec3/256;
    redVec4 = np.array([154,82,61]); colors['red4'] = redVec4/256;
    redVec5 = np.array([123,31,4]); colors['red5'] = redVec5/256;
    #greens  lightest to darkest
    greenVec1 = np.array([142,223,180]); colors['green1'] = greenVec1/256;
    greenVec2 = np.array([89,223,151]); colors['green2'] = greenVec2/256;
    greenVec3 = np.array([16,162,84]); colors['green3'] = greenVec3/256;
    greenVec4 = np.array([43,109,74]); colors['green4'] = greenVec4/256;
    greenVec5 = np.array([3,87,42]); colors['green5'] = greenVec5/256;
    #yellows  lightest to darkest
    yellowVec1 = np.array([246,204,156]); colors['yellow1'] = yellowVec1/256;
    yellowVec2 = np.array([246,185,98]); colors['yellow2'] = yellowVec2/256;
    yellowVec3 = np.array([230,144,23]); colors['yellow3'] = yellowVec3/256;
    yellowVec4 = np.array([154,115,61]); colors['yellow4'] = yellowVec4/256;
    yellowVec5 = np.array([123,74,4]); colors['yellow5'] = yellowVec5/256;
    
    #blue grays
    gBlueVec1 = np.array([197,199,202]); colors['bluegrey1'] = gBlueVec1/256;
    gBlueVec2 = np.array([195,198,202]); colors['bluegrey2'] = gBlueVec2/256;
    gBlueVec3 = np.array([142,145,149]); colors['bluegrey3'] = gBlueVec3/256;
    gBlueVec4 = np.array([108,110,111]); colors['bluegrey4'] = gBlueVec4/256;
    gBlueVec5 = np.array([46,73,97]); colors['bluegrey5'] = gBlueVec5/256;
    #red grays
    gRedVec1 = np.array([242,237,236]); colors['redgrey1'] = gRedVec1/256;
    gRedVec2 = np.array([242,235,233]); colors['redgrey2'] = gRedVec2/256;
    gRedVec3 = np.array([230,231,218]); colors['redgrey3'] = gRedVec3/256;
    gRedVec4 = np.array([172,167,166]); colors['redgrey4'] = gRedVec4/256;
    gRedVec5 = np.array([149,88,71]); colors['redgrey5'] = gRedVec5/256;
    #green grays
    gGreenVec1 = np.array([203,209,206]); colors['greengrey1'] = gGreenVec1/256;
    gGreenVec2 = np.array([201,209,204]); colors['greengrey2'] = gGreenVec2/256;
    gGreenVec3 = np.array([154,162,158]); colors['greengrey3'] = gGreenVec3/256;
    gGreenVec4 = np.array([117,122,119]); colors['greengrey4'] = gGreenVec4/256;
    gGreenVec5 = np.array([50,105,76]); colors['greengrey5'] = gGreenVec5/256;
    #yellow grays
    gYellowVec1 = np.array([242,240,236]); colors['yellowgrey1'] = gYellowVec1/256;
    gYellowVec2 = np.array([242,239,233]); colors['yellowgrey2'] = gYellowVec2/256;
    gYellowVec3 = np.array([230,225,218]); colors['yellowgrey3'] = gYellowVec3/256;
    gYellowVec4 = np.array([172,169,166]); colors['yellowgrey4'] = gYellowVec4/256;
    gYellowVec5 =np.array( [149,117,71]); colors['yellowgrey5'] = gYellowVec5/256;
    
    #pure grays (white to black)
    gVec1 = np.array([256,256,256]); colors['grey1'] = gVec1/256;
    colors['white'] = colors['grey1']
    gVec2 = np.array([242,242,242]); colors['grey2'] = gVec2/256;
    gVec3 = np.array([230,230,230]); colors['grey3'] = gVec3/256;
    gVec4 = np.array([204,204,204]); colors['grey4'] = gVec4/256;
    gVec5 = np.array([179,179,179]); colors['grey5'] = gVec5/256;
    gVec6 = np.array([153,153,153]); colors['grey6'] = gVec6/256;
    gVec7 = np.array([128,128,128]); colors['grey7'] = gVec7/256;
    gVec8 = np.array([102,102,102]); colors['grey8'] = gVec8/256;
    gVec9 = np.array([77,77,77]); colors['grey9'] = gVec9/256;
    gVec10 = np.array([51,51,51]); colors['grey10'] = gVec10/256;
    gVec11 = np.array([26,26,26]); colors['grey11'] = gVec11/256;
    gVec12 = np.array([0,0,0]); colors['grey12'] = gVec12/256;
    colors['black'] = np.array([0,0,0]);
    
    return colors

def physical_constants():

    p = dict(h = 6.62606957e-34,#Planck's constant in kg m^2/s
         hBar = 6.62606957e-34/2/np.pi,
         c = 299792458,#speed of light in meters per second
         epsilon0 = 8.854187817e-12,#permittivity of free space in farads per meter
         mu0 = 4*np.pi*1e-7,#permeability of free space in volt seconds per amp meter
         kB = 1.3806e-23,#Boltzmann's constant
         e = 1.60217657e-19,#electron charge in coulombs
         mE = 9.10938291e-31,#mass of electron in kg
         eV = 1.60217657e-19,#joules per eV
         Ry = 9.10938291e-31*1.60217657e-19**4/(8*8.854187817e-12**2*(6.62606957e-34/2/np.pi)**3*299792458),#13.3*eV;#Rydberg in joules
         a0 = 4*np.pi*8.854187817e-12*(6.62606957e-34/2/np.pi)**2/(9.10938291e-31*1.60217657e-19**2),#estimate of Bohr radius
         Phi0 = 6.62606957e-34/(2*1.60217657e-19),#flux quantum
         Phi0__pH_ns = 6.62606957e3/(2*1.60217657)
         )

    return p 

def set_plot_params(case):
    
    pp = dict()
    
    if case == 'large':
    
        pp['title_font_size'] = 10
        pp['subtitle_font_size'] = 10
        pp['axes_labels_font_size'] = 10
        pp['axes_labels_pad'] = 0 # 4
        pp['tick_labels_font_size'] = 10
        pp['legend_font_size'] = 10
        pp['nominal_linewidth'] = 0.75
        pp['fine_linewidth'] = 0.5
        pp['bold_linewidth'] = 2
        pp['nominal_markersize'] = 2
        pp['big_markersize'] = 3
        tn = 14
        pp['fig_size'] = (tn,tn/1.618)
        # pp['fig_size'] = (tn,tn/1.2)
        pp['axes_linewidth'] = 0.75
        
        pp['major_tick_width'] = 0.75
        pp['major_tick_length'] = 3
        pp['minor_tick_width'] = 0.25
        pp['minor_tick_length'] = 2
        
        pp['xmargin'] = 0 # 0.05 # space between traces and axes
        pp['ymargin'] = 0.05 # 0.05
        
        plt.rcParams['font.family'] = ['sans-serif']
        plt.rcParams['font.sans-serif'] = 'CMU Sans Serif'#'Computer Modern Sans Serif'
        plt.rcParams['figure.figsize'] = pp['fig_size']
        plt.rcParams['figure.titlesize'] = pp['title_font_size']
        plt.rcParams['figure.autolayout'] = True
        
        plt.rcParams['axes.linewidth'] = pp['axes_linewidth']
        plt.rcParams['axes.grid'] = False
        plt.rcParams['axes.titlesize'] = pp['subtitle_font_size']
        plt.rcParams['axes.labelsize'] = pp['axes_labels_font_size']
        plt.rcParams['axes.labelpad'] = pp['axes_labels_pad']
        plt.rcParams['axes.xmargin'] = pp['xmargin']
        plt.rcParams['axes.ymargin'] = pp['ymargin']
        plt.rcParams['axes.titlepad'] = 0
        
        plt.rcParams['legend.fontsize'] = pp['legend_font_size']
        plt.rcParams['legend.loc'] = 'best'
        
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['savefig.format'] = 'pdf'
        plt.rcParams['savefig.pad_inches'] = 0
        
        plt.rcParams['xtick.labelsize'] = pp['tick_labels_font_size']
        plt.rcParams['xtick.direction'] = 'out'
        plt.rcParams['xtick.major.bottom'] = True
        plt.rcParams['xtick.major.top'] = True
        plt.rcParams['xtick.major.size'] = pp['major_tick_length']
        plt.rcParams['xtick.major.width'] = pp['major_tick_width']
        plt.rcParams['xtick.minor.visible'] = True
        plt.rcParams['xtick.minor.size'] = pp['minor_tick_length']
        plt.rcParams['xtick.minor.width'] = pp['minor_tick_width']
        
        plt.rcParams['ytick.labelsize'] = pp['tick_labels_font_size']
        plt.rcParams['ytick.direction'] = 'out'
        plt.rcParams['ytick.major.left'] = True
        plt.rcParams['ytick.major.right'] = True
        plt.rcParams['ytick.major.size'] = pp['major_tick_length']
        plt.rcParams['ytick.major.width'] = pp['major_tick_width']
        plt.rcParams['ytick.minor.visible'] = True
        plt.rcParams['ytick.minor.size'] = pp['minor_tick_length']
        plt.rcParams['ytick.minor.width'] = pp['minor_tick_width']
        
        plt.rcParams['font.family'] = ['sans-serif']
        plt.rcParams['font.sans-serif'] = 'Verdana'#'Computer Modern Sans Serif'
        plt.rcParams['figure.figsize'] = [15,15/1.618]
        plt.rcParams['figure.titlesize'] = 14
        plt.rcParams['axes.titlesize'] = 16
        plt.rcParams['axes.labelsize'] = 14
        plt.rcParams['legend.fontsize'] = 12
        plt.rcParams['xtick.labelsize'] = 14
        plt.rcParams['ytick.labelsize'] = 14
        plt.rcParams['figure.autolayout'] = False
        
    elif case == 'publication':
        
        pp['title_font_size'] = 10
        pp['subtitle_font_size'] = 10
        pp['axes_labels_font_size'] = 10
        pp['axes_labels_pad'] = 0 # 4
        pp['tick_labels_font_size'] = 10
        pp['legend_font_size'] = 10
        pp['nominal_linewidth'] = 0.75
        pp['fine_linewidth'] = 0.5
        pp['bold_linewidth'] = 2
        pp['nominal_markersize'] = 2
        pp['big_markersize'] = 3
        tn = 1.1*8.6/2.54
        pp['fig_size'] = (tn,tn/1.618)
        # pp['fig_size'] = (tn,tn/1.2)
        pp['axes_linewidth'] = 0.75
        
        pp['major_tick_width'] = 0.75
        pp['major_tick_length'] = 3
        pp['minor_tick_width'] = 0.25
        pp['minor_tick_length'] = 2
        
        pp['xmargin'] = 0 # 0.05 # space between traces and axes
        pp['ymargin'] = 0.05 # 0.05
        
        plt.rcParams['font.family'] = ['sans-serif']
        plt.rcParams['font.sans-serif'] = 'CMU Sans Serif'#'Computer Modern Sans Serif'
        plt.rcParams['figure.figsize'] = pp['fig_size']
        plt.rcParams['figure.titlesize'] = pp['title_font_size']
        plt.rcParams['figure.autolayout'] = True
        
        plt.rcParams['axes.linewidth'] = pp['axes_linewidth']
        plt.rcParams['axes.grid'] = False
        plt.rcParams['axes.titlesize'] = pp['subtitle_font_size']
        plt.rcParams['axes.labelsize'] = pp['axes_labels_font_size']
        plt.rcParams['axes.labelpad'] = pp['axes_labels_pad']
        plt.rcParams['axes.xmargin'] = pp['xmargin']
        plt.rcParams['axes.ymargin'] = pp['ymargin']
        plt.rcParams['axes.titlepad'] = 0
        
        plt.rcParams['legend.fontsize'] = pp['legend_font_size']
        plt.rcParams['legend.loc'] = 'best'
        
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['savefig.format'] = 'pdf'
        plt.rcParams['savefig.pad_inches'] = 0
        
        plt.rcParams['xtick.labelsize'] = pp['tick_labels_font_size']
        plt.rcParams['xtick.direction'] = 'out'
        plt.rcParams['xtick.major.bottom'] = True
        plt.rcParams['xtick.major.top'] = True
        plt.rcParams['xtick.major.size'] = pp['major_tick_length']
        plt.rcParams['xtick.major.width'] = pp['major_tick_width']
        plt.rcParams['xtick.minor.visible'] = True
        plt.rcParams['xtick.minor.size'] = pp['minor_tick_length']
        plt.rcParams['xtick.minor.width'] = pp['minor_tick_width']
        
        plt.rcParams['ytick.labelsize'] = pp['tick_labels_font_size']
        plt.rcParams['ytick.direction'] = 'out'
        plt.rcParams['ytick.major.left'] = True
        plt.rcParams['ytick.major.right'] = True
        plt.rcParams['ytick.major.size'] = pp['major_tick_length']
        plt.rcParams['ytick.major.width'] = pp['major_tick_width']
        plt.rcParams['ytick.minor.visible'] = True
        plt.rcParams['ytick.minor.size'] = pp['minor_tick_length']
        plt.rcParams['ytick.minor.width'] = pp['minor_tick_width']
        
        plt.rcParams['font.family'] = ['sans-serif']
        plt.rcParams['font.sans-serif'] = 'Verdana'#'Computer Modern Sans Serif'
        plt.rcParams['figure.figsize'] = [tn,tn/1.618]
        plt.rcParams['figure.titlesize'] = 10
        plt.rcParams['axes.titlesize'] = 10
        plt.rcParams['axes.labelsize'] = 10
        plt.rcParams['legend.fontsize'] = 8
        plt.rcParams['xtick.labelsize'] = 10
        plt.rcParams['ytick.labelsize'] = 10
        plt.rcParams['figure.autolayout'] = False    
    
   
    fig_size = tn
    colors = color_dictionary()
    plt.rcParams['axes.prop_cycle'] = cycler('color', [colors['blue1'],colors['blue2'],colors['blue3'],colors['blue4'],colors['blue5'],
                                                        colors['blue4'],colors['blue3'],colors['blue2'],colors['blue1'],
                                                        colors['red1'],colors['red2'],colors['red3'],colors['red4'],colors['red5'],
                                                        colors['red4'],colors['red3'],colors['red2'],colors['red1'],
                                                        colors['green1'],colors['green2'],colors['green3'],colors['green4'],colors['green5'],
                                                        colors['green4'],colors['green3'],colors['green2'],colors['green1'],
                                                        colors['yellow1'],colors['yellow2'],colors['yellow3'],colors['yellow4'],colors['yellow5'],
                                                        colors['yellow4'],colors['yellow3'],colors['yellow2'],colors['yellow1']])

    return fig_size