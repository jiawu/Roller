import matplotlib.pyplot as plt

def style1():
    plt.rcParams['font.size'] = 24
    plt.rcParams['axes.labelsize'] = 24
    plt.rcParams['xtick.labelsize'] = 24
    plt.rcParams['ytick.labelsize'] = 24
    plt.rcParams['figure.autolayout'] = True

    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['savefig.pad_inches'] = 0.1

    plt.rcParams['figure.facecolor'] = "white"
    plt.rcParams['text.color'] = "black"
    plt.rcParams['axes.labelcolor'] = "black"
    plt.rcParams['legend.frameon'] = "False"
    plt.rcParams['legend.numpoints'] = 1
    plt.rcParams['legend.scatterpoints'] = 1
    plt.rcParams['legend.fontsize'] = 18
    plt.rcParams['legend.markerscale'] = 1
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.color'] = 'black'
    plt.rcParams['ytick.color'] = 'black'
    plt.rcParams['axes.axisbelow'] = True
    plt.rcParams['image.cmap'] = 'Greys'
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ["Arial", "sans-serif"]

    plt.rcParams['axes.grid'] = False
    plt.rcParams['axes.facecolor'] = "white"
    plt.rcParams['axes.edgecolor'] = "black"
    plt.rcParams['axes.linewidth'] = 1.25
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.minor.size'] = 5
    plt.rcParams['ytick.major.size'] = 0
    plt.rcParams['ytick.minor.size'] = 0