import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

"""
- Generates a DSR plot from the first two NPRs for two series of data
- reads in the input csv files from NPR calculation script
"""

### Data Import ###
### Pubchem Sulfoxide Data ###
sf_name = 'known_sulfoxide_nprs.csv'
#pull in data as DF
sf_df = pd.read_csv(f'{sf_name}')

### Generated Sulfoxide Data ###
gen_name = 'generated_sulfoxide_nprs.csv'
gen_df = pd.read_csv(f'{gen_name}')

#pull the NPRs from prepared column data; store to lists
#scifinder data
sf_x1s = sf_df['X1s'].to_list()
sf_y1s = sf_df['Y1s'].to_list()
#generated data
gen_x1s = gen_df['X1s'].to_list()
gen_y1s = gen_df['Y1s'].to_list()

#Initialize plotting
plt.figure(figsize=(18.2,16))
plt.rcParams['axes.linewidth'] = 3.0
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False

#Format the Disk/Linear/Sphere triangle
#vertices of the triangle plot
x1, y1 = [0.5, 0], [0.5, 1]
x2, y2 = [0.5, 1], [0.5, 1]
x3, y3 = [0, 1], [1, 1]
#draw gray lines connecting the above
plt.plot(x1,y1,x2,y2,x3,y3,
    c='gray', lw=3)
#add some text blocks to each corner
plt.text(0, 1.02, s='Linear', fontsize=55, ha='left',va='center', fontweight='bold')
plt.text(1, 1.02, s='Sphere', fontsize=55, ha='right',va='center', fontweight='bold')
plt.text(0.5, 0.472, s='Disk', fontsize=55, ha='center',va='bottom', fontweight='bold')


#plot the gen data
gen_plot = plt.scatter(gen_x1s, gen_y1s,
            color=(1, 0.24, 0.06), label='Generated',
            linewidths=0.1, alpha=0.5, s=65)


#plot the sf data
sf_plot = plt.scatter(sf_x1s, sf_y1s,
            c='b', label='SciFinder',
            linewidths=0.1, alpha=0.7, s=65)


#axis formatting
plt.xlabel('NPR1', fontsize=70, fontweight='bold', labelpad=15)
plt.ylabel('NPR2', fontsize=70, fontweight='bold', labelpad=15)

plt.tick_params(axis='x', labelsize=55, width=8.0, length=10.0)
plt.tick_params(axis='y', labelsize=55, width=8.0, length=10.0)

plt.savefig('PMI_Plot_65_custom_orange3_a7.png', bbox_inches="tight")




