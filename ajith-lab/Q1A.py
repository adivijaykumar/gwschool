import numpy as np
import scipy
import matplotlib.pyplot as plt


def forward(func, x, h):
    return (func(x + h) - func(x)) / h


def backward(func, x, h):
    return (func(x) - func(x - h)) / h


def central(func, x, h):
    return (func(x + h) - func(x - h)) / (2 * h)


def func(x):
    return np.exp(x) * np.sin(x)


def der_func(x):
    return np.exp(x) * (np.sin(x) + np.cos(x))

h = 0.1
f = []
b = []
c = []
x = np.arange(0, 2 * 3.14, 0.1)
print x
deriv = der_func(x)

fig = plt.figure(1)
fig1 = plt.figure(2)

i = 1
for h in [0.5, 0.1, 0.01]:
    f = forward(func, x, h)
    b = backward(func, x, h)
    c = central(func, x, h)

    ax = fig.add_subplot(1, 3, i)
    ax.plot(x, f, label='forward')
    ax.plot(x, b, label='backward')
    ax.plot(x, c, label='central')
    ax.plot(x, deriv, 'ro')
    ax.legend()

    ax1 = fig1.add_subplot(1, 3, i)
    ax1.plot(x, np.abs(f - deriv), label='forward')
    ax1.plot(x, np.abs(b - deriv), label='backward')
    ax1.plot(x, np.abs(c - deriv), label='central')
    ax1.legend()
    i = i + 1

plt.show()
