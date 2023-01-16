"""
Assess the dccm extracted from Rosetta and put in heatmap to demonstrate connectivity
Average across all dccm maps present in the folder
"""

import numpy as np
import matplotlib.pyplot as plt
import sys, os


absolute = sys.argv[1] if len(sys.argv) >= 2 else False

path = os.getcwd()
directory = os.fsencode(path)
DCCM = []
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".dccm"):
        data = np.loadtxt(filename, dtype=float)
        if absolute:
            data = np.abs(data) # take the absolute value to just check for correlations
        DCCM.append(data)
DCCM = np.mean(DCCM, axis=0)

np.savetxt("DCCM_avg.txt", DCCM)

"""
#domains_file=str(sys.argv[1])

# domains should contain the name and resid positions of each domain you want to add into the image
domains_txt = np.loadtxt(domains_file, dtype=str)
domains = []
for D in domains_txt:
    # record domain locations for each chain to be plotted on the heatmap
    # these are pdb numbers, not python numbers
    domains.append([D[0], int(D[1]), int(D[2])])

# PLOTTING #
fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
im = ax.imshow(DCCM, cmap='hot', interpolation='nearest')

cbar = plt.colorbar(im)
cbar.set_label('DCCM', rotation=270)

ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)

no_resid = len(DCCM)

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
"""