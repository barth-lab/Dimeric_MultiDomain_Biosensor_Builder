"""
Assess the dccm extracted from Rosetta and put in heatmap to demonstrate connectivity
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl 
import matplotlib.colors as mcol
import matplotlib.cm as cm

def create_discrete_heatmap(data):
    """
    Create a discrete heatmap of coupling to match with the manuscript figures
    """
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

    # create discrete heatmap
    bounds=[-1, -0.714, -0.429, -0.143, 0.143, 0.429, 0.714, 1]
    # digitisation pushes results forward because of binning
    iratio = np.digitize(data.flat, bounds).reshape(data.shape) -1
    iratio = np.flip(iratio, axis=0)
    
    # create cmap with correct colours versus paper
    colours=["#0101cfff", "#077efcff", "#9ccbffff", "#ffffffff", "#fdfe9aff", "#ff9936ff", "#c90101ff"]
    cmap = mpl.colors.ListedColormap(colours) # 7 discrete colors
    cmap = cm.get_cmap(cmap, lut=len(bounds))
    cmap_bounds = np.arange(len(bounds))
    norm = mcol.BoundaryNorm(cmap_bounds, cmap.N)

    # plot
    map_plot = plt.pcolormesh(iratio, cmap=cmap, norm=norm)
    cbar = plt.colorbar(map_plot, ticks=[0,1,2,3,4,5,6,7], orientation="vertical")
    cbar.set_ticklabels(["-1", "-0.71", "-0.43", "-0.14", "0.14", "0.43", "0.71", "1"])

    return fig, ax, cbar

filename=str(sys.argv[1])
domains_file=str(sys.argv[2])

# domains should contain the name and resid positions of each domain you want to add into the image
domains_txt = np.loadtxt(domains_file, dtype=str)
domains = []
for D in domains_txt:
    # record domain locations for each chain to be plotted on the heatmap
    # these are pdb numbers, not python numbers
    domains.append([D[0], int(D[1]), int(D[2])])

data = np.loadtxt(filename, dtype=float)

fig, ax, cbar = create_discrete_heatmap(data)

#im = ax.imshow(data, cmap='hot', interpolation='nearest')

#cbar = plt.colorbar(im)
cbar.set_label('DCCM', rotation=270)

ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)

no_resid = len(data)

# get figure metadata for arrow width
bbox = ax.get_window_extent()
# dpi used to convert from display units to inches
dpi = fig.dpi
height = bbox.height / dpi  # in inches
width = bbox.width / dpi  # in inches

# annotations for domains
for D in domains:
    start = (D[1]-1) / no_resid; end = (D[2]-1) / no_resid
    mid =  float(start + (end-start)/2)

    # horizontal labels
    ax.annotate('', xy=(start, 1.04), xytext=(end, 1.04), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops={'arrowstyle': '|-|'})
    ax.annotate(D[0], xy=(mid, 1.1), ha='center', va='center', xycoords='axes fraction')
    # used to be textcoords in second ax instead of xycoord

    # vertical labels
    ax.annotate('', xy=(-0.04, -1*(start-1)), xytext=(-0.04, -1*(end-1)), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops={'arrowstyle': '|-|'})
    ax.annotate(D[0], xy=(-0.13, -1*(mid-1)), ha='center', va='center', xycoords='axes fraction')
 
#plt.savefig("test.png")
plt.show()
