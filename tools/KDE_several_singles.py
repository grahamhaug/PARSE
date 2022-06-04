import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm
import math
  
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(-5,18)
ax.set_ylim(0,3)

x_points = [1, 4, 5, 6, 10, 12]
y = [0] * len(x_points)

plt.scatter(x_points, y, color='navy')

for point in x_points:
	mu = point
	variance = 0.5
	sigma = math.sqrt(variance)
	x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
	plt.plot(x, stats.norm.pdf(x, mu, sigma), color='cornflowerblue')

ticks = x_points
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(False)
ax.get_yaxis().set_ticks([])
plt.xticks(x_points, x_points)

plt.show()

