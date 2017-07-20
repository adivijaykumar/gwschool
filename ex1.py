import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as pyplot
import scipy

e, p = np.loadtxt('sly4.dat', comments='#', usecols=(2, 3), unpack=True)

fun = interp1d(np.log(p), np.log(e), fill_value="extrapolate")

e_p = np.exp(fun(np.log(p)))

integrand = 1 / (e_p + p)

h_p = scipy.integrate.cumtrapz(integrand, p, initial=0)

p_h = interp1d(h_p, p, fill_value="extrapolate")


e_h = np.exp(fun(np.log(p_h)))

print 'done'

def dydt(theta, h):
    m, r = theta
    r_prime = (-1) * r * (r - 2 * m) / (m + 4 * np.pi * r**3 * p_h)
    m_prime = 4 * np.pi * r**2 * e_h * r_prime

    return [m_prime, r_prime]
