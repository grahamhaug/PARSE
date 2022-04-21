import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm
import math

"""
Example of how Gaussian KDE works
"""
  
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(-5,18)
ax.set_ylim(0,3)

#example list of points
x_points = [1, 4, 5, 6, 10, 12]
y = [0] * len(x_points)

#plot the points
plt.scatter(x_points, y, color='navy')

#draw a gaussian kernel at each point
coords = {}
for point in x_points:
	mu = point
	variance = 1
	sigma = math.sqrt(variance)
	domain = np.linspace(mu - 3*sigma, mu + 3*sigma, 25)
	fx = stats.norm.pdf(domain, mu, sigma)
	plt.plot(domain, fx, color='cornflowerblue')
	print(len(domain))

	#for summing the gaussian kernels
	# dict key = x, dict value = summed y's
	for i in range(len(domain)):
		x = domain[i]
		y = fx[i]
		if x in coords.keys():
			# add to current sum
			coords[x] = coords[x] + y
		else: # key doesn't exist yet
			coords[x] = y

#plot KDE example (sum of kernels)
gauss_df = pd.DataFrame({'gauss_x':list(coords.keys()),'gauss_fx':list(coords.values())})
plt.plot(gauss_df['gauss_x'], gauss_df['gauss_fx'], 'navy')

#plot formatting
ticks = x_points
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(False)
ax.get_yaxis().set_ticks([])
plt.xticks(x_points, x_points)

#show the plot
plt.show()
