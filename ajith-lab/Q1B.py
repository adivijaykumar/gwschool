import numpy as np
import scipy
import matplotlib.pyplot as plt


def central(func, x, h):
    return (func(x + h) - func(x - h)) / (2 * h)

t, h_plus, h_cross = np.loadtxt('nr_data', unpack=True)
phi = np.arctan(h_cross / h_plus)
plt.plot(t, phi)
plt.show()
