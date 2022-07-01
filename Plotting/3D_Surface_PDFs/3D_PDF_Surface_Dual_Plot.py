#util
import pandas as pd
import numpy as np
#PDF calculation
from scipy.stats import gaussian_kde
#Plotting
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.cm as cm
#plot formatting
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable


### Data Import ###
### SciFinder Sulfoxide Data ###
known_name = 'Processed_PubChem_Known_Sulfoxides_Properties.csv'
raw_known = pd.read_csv(f'{known_name}')
trimmed_known = raw_known.where(raw_known['Cm'] <= 1400).dropna().reset_index(drop=True)
known_df = trimmed_known.sample(frac=0.0025, random_state=2, ignore_index=True)

### Generated Sulfoxides ###
gen_name = 'PubChem_generated_sulfoxides_Properties.csv'
raw_gen = pd.read_csv(f'{gen_name}')
trimmed_gen = raw_gen.where(raw_gen['Cm'] <= 1400).dropna().reset_index(drop=True)
gen_df = trimmed_gen.sample(frac=0.1, random_state=2, ignore_index=True)

### X/Y Data ###
#SF Sulfoxides
a = known_df['FSP3']
b = known_df['Cm']
#Generated Sulfoxides
x = gen_df['FSP3']
y = gen_df['Cm']

### Create meshgrids ###
#SF Sulfoxides
aa, bb = np.mgrid[0:1:500j, 0:1400:1000j]
#Generated Sulfoxides
xx, yy = np.mgrid[0:1:500j, 0:1400:1000j]

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
#fig = plt.figure(figsize=(18, 18))
fig = plt.figure(figsize=(14,8))


### Blue surface (SciFinder data) ###
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax1.plot_surface(aa, bb, g, rstride=1, cstride=1, cmap=cm.winter, alpha=1)
ax1.plot_wireframe(aa, bb, g, rstride=2, cstride=2, color='white', lw=0.05)
ax1.contour(aa, bb, g, 5, zdir='x', offset=1, cmap=cm.Blues)
ax1.contour(aa, bb, g, 7, zdir='y', offset=0, cmap=cm.Blues)

#Ranges
# ax1.set_xlim([0, 1])
ax1.set_xlim([0, 1])
ax1.set_ylim([0, 1400])
ax1.set_zlim([0, 0.009])

#Labels
# ax1.set_title('PDF Plot (SciFinder Sulfoxides)', pad=0)
ax1.set_xlabel('FSP3', fontsize=10, fontweight='bold', labelpad=5)
ax1.set_ylabel('Cm', fontsize=10, fontweight='bold', labelpad=5)
ax1.set_zlabel('PDF', fontsize=10, fontweight='bold', labelpad=10)

#Orientation
ax1.view_init(7, 145)

### Red surface (Generated data) ###
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
ax2.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap=cm.autumn, alpha=1)
ax2.plot_wireframe(xx, yy, f, rstride=2, cstride=2, color='white', lw=0.05)
ax2.contour(xx, yy, f, 7, zdir='x', offset=1, cmap=cm.Reds)
ax2.contour(xx, yy, f, 5, zdir='y', offset=0, cmap=cm.Reds)

### Plot Formatting ###
### Ranges
ax2.set_xlim([0, 1])
ax2.set_ylim([0, 1400])
ax2.set_zlim([0, 0.009])

### Labels
# ax2.set_title('Surface plot of Gaussian 2D KDE (Generated Sulfoxides)', pad=0)
ax2.set_xlabel('FSP3', fontsize=10, fontweight='bold', labelpad=5)
ax2.set_ylabel('Cm', fontsize=10, fontweight='bold', labelpad=5)
ax2.set_zlabel('PDF', fontsize=10, fontweight='bold', labelpad=10)

### Orientation
ax2.view_init(7, 145)


fig.subplots_adjust(wspace=0.15, hspace=0)
#plt.show()
plt.savefig('3D_PDF_Separate_Surfaces.png', bbox_inches="tight")




