import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from scipy.stats import gaussian_kde
from matplotlib.colors import LinearSegmentedColormap

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
fig = plt.figure(figsize=(24,16))
ax = fig.add_subplot(projection='3d')
ax.tick_params(grid_color='r')

my_yrs = plt.cm.YlOrRd(np.linspace(0.3,0.95,100))
my_yrs = ListedColormap(my_yrs[8:,:-1])


blue_colors = [(0, 0, 1), (0, 0.11, 0.38)] # Experiment with this
my_blues = LinearSegmentedColormap.from_list('test', blue_colors, N=100)

# red_colors = [(1, 1, 0.2), (0.4, 0, 0)] # Experiment with this
# my_reds = LinearSegmentedColormap.from_list('test2', red_colors, N=100)

#plot the 3param data for each set
known_plot = ax.scatter(known_x, known_y, known_z,
            c=known_dens, cmap=my_blues, label='Known',
            linewidths=0.1, alpha=0.5, s=30)

gen_plot = ax.scatter(gen_x, gen_y, gen_z,
            c=gen_dens, cmap=my_yrs, label='Generated',
            linewidths=0.1, alpha=0.5, s=30)

#ranges for ticks
x_ticks = np.arange(0.0, 1700, 400)
y_ticks = np.arange(0.00, 1.1, 0.2)
z_ticks = np.arange(0.0, 1700, 400)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.set_zticks(z_ticks)

#Axes labels and formatting
ax.xaxis.set_rotate_label(False)
ax.yaxis.set_rotate_label(False)
ax.zaxis.set_rotate_label(False)

ax.set_xlabel(r'$\mathbf{MW}$', fontsize=40, fontweight='bold', labelpad=50)
ax.set_ylabel(r'$\mathbf{Fsp^3}$', fontsize=40, fontweight='bold', labelpad=40)
ax.set_zlabel(r'$\mathbf{C_m}$', fontsize=40, fontweight='bold', labelpad=50)

#X ticks
ax.tick_params(axis='x', labelsize=30, pad=30)
x_labels = [0, 400, 800, 1200, 1600]
ax.set_xticklabels(x_labels, rotation=3, va='bottom', ha='center')

#Y ticks
ax.tick_params(axis='y', labelsize=30, pad=25)
y_labels = [0.00, 0.2, 0.4, 0.6, 0.8, 1.0]
ax.set_yticklabels(y_labels, rotation=3, va='bottom', ha='center')

ax.tick_params(axis='z', labelsize=30)
z_labels = [0, 200, 600, 1000, 1400]
ax.set_zticklabels(z_labels, rotation=3, va='center', ha='left')

#view angle for 3D Plot
ax.view_init(15, 115)

# ### Colorbars ###
# #Color ranges
# norm1 = mpl.colors.Normalize(vmin=0.00001, vmax=0.00008)
# sm = plt.cm.ScalarMappable(cmap=my_blues, norm=norm1)
# sm.set_array([])
# sf_cb = plt.colorbar(sm, shrink=0.65, aspect=10, location='left', pad=-0.05)

# norm2 = mpl.colors.Normalize(vmin=0.000005, vmax=0.000025)
# sm2 = plt.cm.ScalarMappable(cmap=my_yrs, norm=norm2)
# sm2.set_array([])
# gen_cb = plt.colorbar(sm2, shrink=0.65, aspect=10, location='left', pad=.1)
# gen_cb.ax.tick_params(labelsize=25, right=True, labelright=True, left=False, labelleft=False)

# #cb tick formatting
# sf_cb.ax.tick_params(labelsize=25)
# sf_cb.set_label('Probability Density Function (PDF)', fontsize=25, labelpad=20)

# #size and format of sci-not
# sf_cb.ax.yaxis.get_offset_text().set_fontsize(20)
# sf_cb.ax.ticklabel_format(useOffset=True, style='sci', useMathText=True)
# gen_cb.ax.yaxis.get_offset_text().set_fontsize(20)
# gen_cb.ax.ticklabel_format(useOffset=True, style='sci', useMathText=True)

#Legend parameters
# import matplotlib.lines as mlines
# sfLeg = mlines.Line2D([], [], color='blue', marker='o', markersize=12, ls='', label='SciFinder')
# genLeg = mlines.Line2D([], [], color='red', marker='o', markersize=12, ls='', label='Generated')
# ax.legend(handles=[sfLeg, genLeg], fontsize=20, edgecolor='k', loc='upper right', bbox_to_anchor=(0.9, 0.76))

#save the fig
plt.savefig('150-si-altangle2.png', bbox_inches="tight")

### For Multiplotting ###
# first = 10
# # second = list(range(0, 160, 10))
# # #np.arange(0, 360, 10)
# second = [110, 130, 160, 210, 240, 270]

# for thing in second:
#     ax.view_init(first, thing)
#     plt.savefig(f'truncated-density-final_{thing}x{first}.png', bbox_inches="tight")
