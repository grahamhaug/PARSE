import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from scipy.stats import gaussian_kde
from matplotlib.colors import LinearSegmentedColormap

### Data Import ###
### SciFinder Sulfoxide Data ###
sf_name = 'properties_scifinder_sulfoxides.csv'
sf_df = pd.read_csv(f'{sf_name}')

### Generated Sulfoxides ###
gen_name = 'generated_sulfoxides_calcprops.csv'
gen_df = pd.read_csv(f'{gen_name}')
trimmed_df = gen_df.where(gen_df['Cm'] <= 1400).dropna().reset_index(drop=True)


### Calculate point density ###
#scifinder set
sf_x = sf_df['MW']
sf_y = sf_df['FSP3']
sf_z = sf_df['Cm']
sf_xyz = np.vstack([sf_x, sf_y, sf_z])
sf_dens = gaussian_kde(sf_xyz)(sf_xyz)
#sort points by density
idx = sf_dens.argsort()
sf_x, sf_y, sf_z, sf_dens = sf_x[idx], sf_y[idx], sf_z[idx], sf_dens[idx]

#generated set
gen_x = trimmed_df['MW']
gen_y = trimmed_df['FSP3']
gen_z = trimmed_df['Cm']
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

my_yrs = plt.cm.YlOrRd(np.linspace(0.2,1,100))
my_yrs = ListedColormap(my_yrs[8:,:-1])

blue_colors = [(0, 0, 1), (0, 0.11, 0.38)] # Experiment with this
my_blues = LinearSegmentedColormap.from_list('test', blue_colors, N=100)

#plot the 3param data for each set
sf_plot = ax.scatter(sf_x, sf_y, sf_z,
            c=sf_dens, cmap=my_blues, label='SciFinder',
            linewidths=0.1, alpha=0.5, s=15)

gen_plot = ax.scatter(gen_x, gen_y, gen_z,
            c=gen_dens, cmap=my_yrs, label='Generated',
            linewidths=0.1, alpha=0.5, s=15)

#ranges for ticks
x_ticks = np.arange(0.0, 1700, 400)
y_ticks = np.arange(0.00, 1.1, 0.2)
z_ticks = np.arange(0.0, 1600, 200)
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
z_labels = [0, 200, 400, 600, 800, 1000, 1200, 1400]
ax.set_zticklabels(z_labels, rotation=3, va='center', ha='right')

#view angle for 3D Plot
ax.view_init(15, 220)

### Colorbars ###
#Color ranges
norm1 = mpl.colors.Normalize(vmin=0.00001, vmax=0.00008)
sm = plt.cm.ScalarMappable(cmap=my_blues, norm=norm1)
sm.set_array([])
sf_cb = plt.colorbar(sm, shrink=0.45, aspect=30, location='top', pad=-0.11)

norm2 = mpl.colors.Normalize(vmin=0.000005, vmax=0.000025)
sm2 = plt.cm.ScalarMappable(cmap=my_yrs, norm=norm2)
sm2.set_array([])
gen_cb = plt.colorbar(sm2, shrink=0.45, aspect=30, location='top', pad=.005)
tick_locator = ticker.MaxNLocator(nbins=7)
#sf_cb.locator = tick_locator
#gen_cb.locator = tick_locator
#gen_cb.update_ticks()
gen_cb.ax.tick_params(labelsize=25, bottom=True, labelbottom=True,
                     top=False, labeltop=False,
                      left=False, labelleft=False)

#cb tick formatting
sf_cb.ax.tick_params(labelsize=25)
sf_cb.set_label('Probability Density Function (PDF)', fontsize=25, labelpad=20)

#size and format of sci-not
#sf_cb.ax.yaxis.get_offset_text().set_fontsize(20)
#sf_cb.ax.ticklabel_format(useOffset=True, style='sci', useMathText=True)
#gen_cb.ax.yaxis.get_offset_text().set_fontsize(20)
sf_cb.ax.xaxis.offsetText.set_visible(False)
gen_cb.ax.xaxis.offsetText.set_visible(False)
#gen_cb.ax.ticklabel_format(useOffset=True, style='sci', useMathText=True)



#save the fig
plt.savefig('final_colorbars_above.png', bbox_inches="tight")


