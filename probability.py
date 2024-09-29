import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm



x_axis = np.arange(-5, 5, 0.01)

mean = 0
sd = 1

plt.plot(x_axis, norm.pdf(x_axis, mean, sd))
plt.show()