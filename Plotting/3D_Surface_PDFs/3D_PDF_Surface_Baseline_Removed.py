import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap
from scipy.stats import gaussian_kde
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits import mplot3d
import matplotlib.colors

### Data Import ###
### PubChem Known Sulfoxide Data ###
known_name = 'Processed_PubChem_Known_Sulfoxides_Properties.csv'
raw_known = pd.read_csv(f'{known_name}')
trimmed_known = raw_known.where(raw_known['Cm'] <= 1400).dropna().reset_index(drop=True)
known_df = trimmed_known.sample(frac=0.0025, random_state=1, ignore_index=True)

### Generated Sulfoxides ###
gen_name = 'PubChem_generated_sulfoxides_Properties.csv'
raw_gen = pd.read_csv(f'{gen_name}')
trimmed_gen = raw_gen.where(raw_gen['Cm'] <= 1400).dropna().reset_index(drop=True)
gen_df = trimmed_gen.sample(frac=0.1, random_state=1, ignore_index=True)

### X/Y Data ###
#SF Sulfoxides
a = known_df['FSP3']
b = known_df['Cm']
#Generated Sulfoxides
x = gen_df['FSP3']
y = gen_df['Cm']

# Create meshgrids
#SF Sulfoxides
aa, bb = np.mgrid[0:1:500j, 0:1000:1000j]
#Generated Sulfoxides
xx, yy = np.mgrid[0:1:500j, 100:1200:1000j]

### Calculate KDE ###
#SF Sulfoxides
sf_positions = np.vstack([aa.ravel(), bb.ravel()])
sf_values = np.vstack([a, b])
sf_kernel = gaussian_kde(sf_values)
g = np.reshape(sf_kernel(sf_positions).T, aa.shape)

#Generated Sulfoxides
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([x, y])
kernel = gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)

### Plotting ###
fig = plt.figure(figsize=(24, 14))
ax = plt.axes(projection='3d')

#Blue surface (SciFinder data)
aa2, bb2, g2 = aa.flatten(), bb.flatten(), g.flatten()
usable_points = (0.001 < g2) & (g2 < 0.012)
aa2, bb2, g2 = aa2[usable_points], bb2[usable_points], g2[usable_points]
ax.plot_trisurf(aa2, bb2, g2, cmap=cm.winter, linewidth=0.05, edgecolor='white', alpha=1)
ax.contour(aa, bb, g, zdir='y', offset=0, cmap=cm.Blues)
ax.contour(aa, bb, g, zdir='x', offset=1, cmap=cm.Blues)


#Red surface (Generated data)
xx2, yy2, f2 = xx.flatten(), yy.flatten(), f.flatten()
usable_points = (0.001 < f2) & (f2 < 0.012)
xx2, yy2, f2 = xx2[usable_points], yy2[usable_points], f2[usable_points]
ax.plot_trisurf(xx2, yy2, f2, cmap=cm.autumn, linewidth=0.03, edgecolor='white', alpha=1)
ax.contour(xx, yy, f, 10, zdir='y', offset=0, cmap=cm.Reds)
ax.contour(xx, yy, f, 10, zdir='x', offset=1, cmap=cm.Reds)

#Axes labels and formatting
ax.xaxis.set_rotate_label(False)
ax.yaxis.set_rotate_label(False)
ax.zaxis.set_rotate_label(False)
ax.set_xlabel(r'$\mathbf{Fsp^3}$', fontsize=30, fontweight='bold', labelpad=25)
ax.set_ylabel(r'$\mathbf{C_m}$', fontsize=30, fontweight='bold', labelpad=20)
ax.set_zlabel(r'$\mathbf{PDF}$', fontsize=30, fontweight='bold', labelpad=45)

#ranges for ticks
x_ticks = np.arange(0.0, 1.2, 0.2)
y_ticks = np.arange(0, 1400, 200)
z_ticks = np.arange(0.0, 0.014, 0.002)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.set_zticks(z_ticks)

#X ticks
ax.tick_params(axis='x', labelsize=20, pad=15)
x_labels = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
ax.set_xticklabels(x_labels, rotation=3, va='bottom', ha='center')

#Y ticks
ax.tick_params(axis='y', labelsize=20, pad=15)
y_labels = [0, 200, 400, 600, 800, 1000, 1200]
ax.set_yticklabels(y_labels, rotation=3, va='bottom', ha='center')

ax.tick_params(axis='z', labelsize=20, pad=15)
# z_labels = [0.000, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012]
# ax.set_zticklabels(z_labels, rotation=3, va='center', ha='left')


ax.view_init(12, 145)
# ax.set_xlim([0, 1])
# ax.set_ylim([0, 1200])
ax.set_zlim([0, 0.008])
#plt.show()

### Colorbars ###
#Color ranges
norm1 = mpl.colors.Normalize(vmin=0.0, vmax=0.012)
sm = plt.cm.ScalarMappable(cmap=cm.winter, norm=norm1)
sm.set_array([])
sf_cb = plt.colorbar(sm, shrink=0.45, aspect=11, location='left', pad=-0.1)

norm2 = mpl.colors.Normalize(vmin=0.0, vmax=0.06)
sm2 = plt.cm.ScalarMappable(cmap=cm.autumn, norm=norm2)
sm2.set_array([])
gen_cb = plt.colorbar(sm2, shrink=0.45, aspect=11, location='left', pad=.05)
gen_cb.ax.tick_params(labelsize=20, right=True, labelright=True, left=False, labelleft=False)

#cb tick formatting
sf_cb.ax.tick_params(labelsize=20)
sf_cb.set_label('Probability Density Function (PDF)', fontweight='bold', fontsize=20, labelpad=20)

#size and format of sci-not
sf_cb.ax.yaxis.get_offset_text().set_fontsize(20)
sf_cb.ax.ticklabel_format(useOffset=True, style='sci', useMathText=True)
gen_cb.ax.yaxis.get_offset_text().set_fontsize(20)
gen_cb.ax.ticklabel_format(useOffset=True, style='sci', useMathText=True)

plt.savefig('3d_pdf_surface_no_base.png', bbox_inches="tight")




