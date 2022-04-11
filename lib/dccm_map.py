import numpy as np
import matplotlib.pyplot as plt

# record domain locations for each chain to be plotted on the heatmap
# these are pdb numbers, not python numbers
D1_A = ["D1_A", 1, 101]
D2_A = ["D2_A", 102, 202]
D3_A = ["D2_A", 214, 235]
D1_B = ["D1_B", 236, 336]
D2_B = ["D2_B", 337, 438]
D3_B = ["D3_B", 449, 470]

domains = [D1_A, D2_A, D3_A, D1_B, D2_B, D3_B]

data = np.loadtxt("dccm_matrix.txt", dtype=float)

fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
im = ax.imshow(data, cmap='hot', interpolation='nearest')

cbar = plt.colorbar(im)
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
    ax.annotate(D[0], xy=(mid, 1.1), ha='center', va='center', textcoords='axes fraction')

    # vertical labels
    ax.annotate('', xy=(-0.04, start), xytext=(-0.04, end), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops={'arrowstyle': '|-|'})
    ax.annotate(D[0], xy=(-0.13, mid), ha='center', va='center', textcoords='axes fraction')
 
    #ax.annotate(D[0], xy=(mid, 1.0), xytext=(mid, 1.10), xycoords='axes fraction', 
    #        fontsize=11, ha='center', va='bottom',
    #        bbox=dict(boxstyle='square', fc='white'),
    #        arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=1.5'%((end-start)), lw=2.0))
#plt.savefig("test.png")
plt.show()
