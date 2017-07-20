import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import scipy

e, p = np.loadtxt('sly4.dat', comments='#', usecols=(2, 3), unpack=True)

G = 6.67 * 10**-8
c = 3 * 10**10

e = G / c**4 * e
p = G / c**4 * p
fun = interp1d(np.log(p), np.log(e), fill_value="extrapolate")

e_p = np.exp(fun(np.log(p)))

integrand = 1 / (e_p + p)

h_p = scipy.integrate.cumtrapz(integrand, p, initial=0)

p_h = interp1d(h_p, p, fill_value="extrapolate")


e_h1 = np.exp(fun(np.log(p_h(h_p))))
e_h = interp1d(h_p, e_h1, fill_value="extrapolate")

print 'done'


def dydt(theta, h):
    m, r = theta
    r_prime = (-1) * r * (r - 2 * m) / (m + 4 * np.pi * r**3 * p_h(h))
    m_prime = 4 * np.pi * r**2 * e_h(h) * r_prime

    return [m_prime, r_prime]

sol = scipy.integrate.odeint(dydt, [2, 6], h_p)

m_arr = []
r_arr = []
for arr in sol:
    print arr[01]
    m_arr.append(arr[0])
    r_arr.append(arr[1])

plt.plot(m_arr, r_arr)
plt.show()
