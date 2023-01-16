"""
Directly compare two DCCM maps - you need to be carful if the domain sizes are different, or if there are extra residues
In other words, your selection of residues to pick needs to match between two maps
"""

import numpy as np
import sys, os
import matplotlib.pyplot as plt
from pylab import *
import matplotlib as mpl 
import matplotlib.colors as mcol
import matplotlib.cm as cm

filename1 = sys.argv[1]
filename2 = sys.argv[2]

map1 = np.loadtxt(filename1, dtype=float)
map2 = np.loadtxt(filename2, dtype=float)

# create absolute values to compare just correlations
map1 = np.abs(map1) # must be small linker
map2 = np.abs(map2) # must be extra helix

# which residue range do you want to take from each map?
# Base this on the domains.txt file which originally comes from the verbose fasta
# These resid are unique to CSF1R
# small linker
slice1 = np.asarray([214, 235, 449, 470]) - 1 # -1 for pythonic indexing
#slice1 = np.asarray([231, 252, 483, 504]) - 1
# extra helix
slice2 = np.asarray([231, 252, 483, 504]) - 1

big_domains = [["D1_A", 1, 101],
               ["D2_A", 102, 203],
               ["D1_B", 204, 304],
               ["D2_B", 304, 405]]

TM_domains = [["TM_A", 1, 22],
              ["TM_B", 22, 43]]

D2_point = 203
map1_len = len(map1)
map2_len = len(map2)
# just compare coupling to end of relevent domains (can't compare like for like with extra helix)
map1_0_A = map1[:D2_point, slice1[0]:slice1[1]]
map1_0_B = map1[:D2_point, slice1[2]:slice1[3]]
map1_1_A = map1[int(map1_len/2):int(map1_len/2)+D2_point, slice1[0]:slice1[1]]
map1_1_B = map1[int(map1_len/2):int(map1_len/2)+D2_point, slice1[2]:slice1[3]]
map2_0_A = map2[:D2_point, slice2[0]:slice2[1]]
map2_0_B = map2[:D2_point, slice2[2]:slice2[3]]
map2_1_A = map2[int(map2_len/2):int(map2_len/2)+D2_point, slice2[0]:slice2[1]]
map2_1_B = map2[int(map2_len/2):int(map2_len/2)+D2_point, slice2[2]:slice2[3]]

# now compare
map_0_A = map1_0_A - map2_0_A
map_0_B = map1_0_B - map2_0_B
map_1_A = map1_1_A - map2_1_A
map_1_B = map1_1_B - map2_1_B

# cat to plot
map_0 = np.concatenate((map_0_A, map_0_B), axis=0)
map_1 = np.concatenate((map_1_A, map_1_B), axis=0)
map_total = np.concatenate((map_0, map_1), axis=1)

# percentage change
map1_total = np.concatenate((map1_0_A, map1_0_B, map1_1_A, map1_1_B))
map2_total = np.concatenate((map2_0_A, map2_0_B, map2_1_A, map2_1_B))
print("Total DCCM change is: %.2f"%(np.sum(map_total)))
print("Normalised DCCM change is: %.2f"%(np.sum(map_total) / (np.shape(map_total)[0]*np.shape(map_total)[1])))
non_neglig = map_total[~np.logical_and(map_total > -0.14, map_total < 0.14)]
print("Non negligable is: %.2f"%(np.sum(non_neglig) / len(non_neglig)))
print("Percentage change is: %.2f"%(np.sum(map_total)*100/np.sum(map2_total)))

# plot
fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

# create discrete heatmap
bounds=[-1, -0.714, -0.429, -0.143, 0.143, 0.429, 0.714, 1]
# digitisation pushes results forward because of binning
iratio = np.digitize(map_total.flat, bounds).reshape(map_total.shape) -1
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
cbar.set_label(r'$\Delta$ |DCCM|', rotation=270, labelpad=20)
cbar.set_ticklabels(["-1", "-0.71", "-0.43", "-0.14", "0.14", "0.43", "0.71", "1"])

ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)

# get figure metadata for arrow width
bbox = ax.get_window_extent()
# dpi used to convert from display units to inches
dpi = fig.dpi
height = bbox.height / dpi  # in inches
width = bbox.width / dpi  # in inches

# annotations for domains
no_resid = np.shape(iratio)[0]
for D in big_domains:

    start = (D[1]-1) / no_resid; end = (D[2]-1) / no_resid
    mid =  float(start + (end-start)/2)

    # vertical labels
    ax.annotate('', xy=(-0.04, -1*(start-1)), xytext=(-0.04, -1*(end-1)), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops={'arrowstyle': '|-|'})
    ax.annotate(D[0], xy=(-0.13, -1*(mid-1)), ha='center', va='center', xycoords='axes fraction')

no_resid = np.shape(iratio)[1]
for D in TM_domains:
    start = (D[1]-1) / no_resid; end = (D[2]-1) / no_resid
    mid =  float(start + (end-start)/2)

    # horizontal labels
    ax.annotate('', xy=(start, 1.04), xytext=(end, 1.04), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops={'arrowstyle': '|-|'})
    ax.annotate(D[0], xy=(mid, 1.1), ha='center', va='center', xycoords='axes fraction')

plt.savefig("Compare_CSF1R_linker_coupling")
plt.show()

"""
Old generic code for heatmap plot of raw coupling

fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
im = ax.imshow(map_total, cmap="hot", interpolation="nearest")

cbar = plt.colorbar(im)
cbar.set_label('DCCM', rotation=270)

ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)

no_resid = np.shape(map_total)[1]

# get figure metadata for arrow width
bbox = ax.get_window_extent()
# dpi used to convert from display units to inches
dpi = fig.dpi
height = bbox.height / dpi  # in inches
width = bbox.width / dpi  # in inches
plt.show()
"""
