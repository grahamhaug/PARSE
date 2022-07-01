import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scipy
import scipy.stats as sts

"""
- Plot histograms of the min/avg/max pairwise distance between sulfinates
- for k-means and randomly selected sulfinates 
- Reads two csvs (the summary output files containing PD data for kmeans/random gen)
- Outputs a figure of 6 subplots with a histogram for each metric
"""

#import data
import_file_name = 'Random_Sampling_Performance-min.csv'
incoming_data = pd.read_csv(f'{import_file_name}')

import_file_name2 = 'KMeans_Sampling_Performance-minmaxed.csv'
incoming_data2 = pd.read_csv(f'{import_file_name2}')

#minimum data
x1 = incoming_data['PD_Min']
mu1 = np.mean(x1)
sigma1 = np.std(x1)
#kmeans
z1 = incoming_data2['PD_Min']

#average data
x2 = incoming_data['PD_Avg']
mu2 = np.mean(x2) #mean
sigma2 = np.std(x2) #stdev
#kmeans
z2 = incoming_data2['PD_Avg']

#maximum data
x3 = incoming_data['PD_Max']
mu3 = np.mean(x3) #mean
sigma3 = np.std(x3) #stdev
#kmeans
z3 = incoming_data2['PD_Max']

#plot histograms
#fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(6,10))
plt.figure(figsize=(12,10))
#fig.tight_layout()
num_bins = 80

#plot the minimum pairwise data
#Random Min. distance
plt.subplot(3,2,1)
#random data with curve
n, bins, patches = plt.hist(x1, 40, range=[-0.1, 0.3], label='Random Selection',
                          facecolor='blue', alpha=0.25, density=True)
# y = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
#      np.exp(-0.5 * (1 / sigma1 * (bins - mu1))**2))
# plt.plot(bins, y, c='navy')
plt.xlabel('Standardized Euclidean Distance')
plt.ylabel('Probability Density')
plt.title('Minimum Pairwise Distance: Randomly Selected Sulfinates')
plt.legend()

#K-means Min. PD
plt.subplot(3,2,2)
n, bins, patches = plt.hist(z1, 30, range=[-0.1, 0.3], label='Clustered Selection',
                          facecolor='green', alpha=0.5, density=True)
plt.xlabel('Standardized Euclidean Distance')
plt.ylabel('Probability Density')
plt.title('Minimum Pairwise Distance: KMeans-Selected Sulfinates')
plt.legend()

#Random Avg. Pd
plt.subplot(3,2,3)
n, bins, patches = plt.hist(x2, num_bins, range=[0.2, 0.8], label='Random Selection',
                          facecolor='blue', alpha=0.25, density=True)
# y = ((1 / (np.sqrt(2 * np.pi) * sigma2)) *
#      np.exp(-0.5 * (1 / sigma2 * (bins - mu2))**2))
# plt.plot(bins, y, c='navy')
plt.xlabel('Standardized Euclidean Distance')
plt.ylabel('Probability Density')
plt.title('Average Pairwise Distance: Randomly Selected Sulfinates')
plt.legend()

#K-means Avg. PD
plt.subplot(3,2,4)
n, bins, patches = plt.hist(z2, num_bins, range=[0.2, 0.8], label='Clustered Selection',
                          facecolor='green', alpha=0.5, density=True)
plt.xlabel('Standardized Euclidean Distance')
plt.ylabel('Probability Density')
plt.title('Average Pairwise Distance: KMeans-Selected Sulfinates')
plt.legend()

#Random Max PD
plt.subplot(3,2,5)
n, bins, patches = plt.hist(x3, 35, range=[0.8, 1.6], label='Random Selection',
                          facecolor='blue', alpha=0.25, density=True)
# y = ((1 / (np.sqrt(2 * np.pi) * sigma3)) *
#      np.exp(-0.5 * (1 / sigma3 * (bins - mu3))**2))
# plt.plot(bins, y, c='navy')
plt.xlabel('Standardized Euclidean Distance')
plt.ylabel('Probability density')
plt.title('Maximum Pairwise Distance: Randomly Selected Sulfinates')
plt.legend()

#K-Means Max PD
plt.subplot(3,2,6)
n, bins, patches = plt.hist(z3, 35, range=[0.8, 1.6], label='Clustered Selection', 
                         facecolor='green', alpha=0.5, density=True)
plt.xlabel('Standardized Euclidean Distance')
plt.ylabel('Probability Density')
plt.title('Maximum Pairwise Distance: KMeans-Selected Sulfinates')
plt.legend()

#add a little vertical space between subplots
plt.subplots_adjust(hspace=0.45)

plt.show()