import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm
import math
  
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(-5,18)
ax.set_ylim(0,3)

# x_points = [1, 4, 5, 6, 10, 12]
# y = [0] * len(x_points)

# plt.scatter(x_points, y, color='navy')


mu = 1
variance = 1
sigma = math.sqrt(variance)
x = np.linspace(mu - 3*sigma, mu + 3*sigma, 25)
fx = stats.norm.pdf(x, mu, sigma)
plt.plot(x, fx, color='cornflowerblue')
gauss_df = pd.DataFrame(columns=['X', 'Fx'])
gauss_df['X'] = x
gauss_df['Fx'] = fx

print(gauss_df)




# print(list_of_xs)
# print(list_of_ys)

# np_ys = np.array(list_of_ys)
# summed_ys = np_ys.sum(axis=0)
# print("")
# print(summed_ys)

# plt.plot(x_points, summed_ys, 'k')

# ticks = x_points
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(False)
# ax.get_yaxis().set_ticks([])
#plt.xticks(x_points, x_points)

plt.show()



