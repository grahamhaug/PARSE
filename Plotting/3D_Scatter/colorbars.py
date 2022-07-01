import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from scipy.stats import gaussian_kde
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

### Data Import ###
### SciFinder Sulfoxide Data ###
known_name = 'Processed_PubChem_Known_Sulfoxides_Properties.csv'
raw_known = pd.read_csv(f'{known_name}')
trimmed_known = raw_known.where(raw_known['Cm'] <= 1400).dropna().reset_index(drop=True)
known_df = trimmed_known.sample(frac=0.0025, random_state=1, ignore_index=True)

### Generated Sulfoxides ###
gen_name = 'PubChem_generated_sulfoxides_Properties.csv'
raw_gen = pd.read_csv(f'{gen_name}')
trimmed_gen = raw_gen.where(raw_gen['Cm'] <= 1400).dropna().reset_index(drop=True)
gen_df = trimmed_gen.sample(frac=0.1, random_state=1, ignore_index=True)


### Calculate point density ###
#scifinder set
known_x = known_df['MW']
known_y = known_df['FSP3']
known_z = known_df['Cm']
known_xyz = np.vstack([known_x, known_y, known_z])
known_dens = gaussian_kde(known_xyz)(known_xyz)
#sort points by density
idx = known_dens.argsort()
known_x, known_y, known_z, known_dens = known_x[idx], known_y[idx], known_z[idx], known_dens[idx]

#generated set
gen_x = gen_df['MW']
gen_y = gen_df['FSP3']
gen_z = gen_df['Cm']
gen_xyz = np.vstack([gen_x, gen_y, gen_z])
gen_dens = gaussian_kde(gen_xyz)(gen_xyz)
#sort points by density
idz = gen_dens.argsort()
gen_x, gen_y, gen_z, gen_dens = gen_x[idz], gen_y[idz], gen_z[idz], gen_dens[idz]


### Plotting Data ###
# #initialize a 3D plot
fig = plt.figure(figsize=(24,24))
#fig, axes = plt.subplots(1,2)

#ax = fig.add_subplot(projection='3d')
#ax.tick_params(grid_color='r')

ax1 = fig.add_subplot(111)
divider = make_axes_locatable(ax1)
cax = divider.append_axes('left', size='20%', pad=1)
ax1.remove()

my_yrs = plt.cm.YlOrRd(np.linspace(0.3,0.95,100))
my_yrs = ListedColormap(my_yrs[8:,:-1])

blue_colors = [(0, 0, 1), (0, 0.11, 0.38)] # Experiment with this
my_blues = LinearSegmentedColormap.from_list('test', blue_colors, N=100)

### Colorbars ###
#Color ranges
norm1 = mpl.colors.Normalize(vmin=0.00001, vmax=0.00008)
sm = plt.cm.ScalarMappable(cmap=my_blues, norm=norm1)
sm.set_array([])
sf_cb = fig.colorbar(sm, cax=cax, shrink=0.65, aspect=10)

#remove ticks and labels
sf_cb.set_ticks([])

plt.savefig('blue_colorbar.png', bbox_inches="tight")


fig = plt.figure(figsize=(24,24))
ax2=fig.add_subplot(111)

divider = make_axes_locatable(ax2)
cax = divider.append_axes('left', size='20%', pad=1)
#save the fig

ax2.remove()

norm2 = mpl.colors.Normalize(vmin=0.000005, vmax=0.000025)
sm2 = plt.cm.ScalarMappable(cmap=my_yrs, norm=norm2)
sm2.set_array([])
gen_cb = fig.colorbar(sm2, cax=cax, shrink=0.65, aspect=10)
# gen_cb.ax.tick_params(labelsize=25, right=True, labelright=True, left=False, labelleft=False)

#remove the ticks/labels
gen_cb.set_ticks([])

#save the fig
plt.savefig('red_colorbar.png', bbox_inches="tight")

