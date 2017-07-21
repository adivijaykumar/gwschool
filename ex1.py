import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import scipy

e, p = np.loadtxt('sly4.dat', comments='#', usecols=(2, 3), unpack=True)

G = 6.67 * 10**-8
c = 3 * 10**10

e = 1.6199 * 10**(-18) * e
p = 1.8063 * 10**(-39) * p

fun = interp1d(np.log(p), np.log(e), fill_value="extrapolate")
# e_p = interp1d(p, e, fill_value="extrapolate")
e_p = np.exp(fun(np.log(p)))
e_p = interp1d(p, e_p, fill_value="extrapolate")


def integrand(z):
    return 1 / (e_p(z) + z)

#h_p = []
# for x in p:
#    h_p.append(scipy.integrate.quad(integrand, 0, x)[0])

h_p = scipy.integrate.cumtrapz(integrand(p), p, initial=0)
p_h = interp1d(h_p, p, fill_value="extrapolate")


# e_h1 = np.exp(fun(np.log(p_h(h_p))))
e_h = interp1d(h_p, e_p(p), fill_value="extrapolate")

print 'done'


def dydt(theta, h):
    m, r = theta
    r_prime = (-1) * r * (r - 2 * m) / (m + 4 * np.pi * r**3 * p_h(h))
    m_prime = 4 * np.pi * r**2 * e_h(h) * r_prime
#    m_prime = 4 * np.pi * r**2 * e_p(p_h(h)) * r_prime

    return [m_prime, r_prime]

m_arr = []
r_arr = []

h_c = 1.0
h_arr = np.linspace(0.1, h_c, num=40)
for h in h_arr:
    sol = scipy.integrate.odeint(
        dydt, [2 * 10 ** (-31), 10**(-10)], np.arange(h, 0, -0.001))

#    if np.amax(sol[:, 0]) == sol[:, -1]:
#        print 'yes'
    m_arr.append(np.amax(sol[:, 0]))
    r_arr.append((10.0/6.67)*np.amax(sol[:, 1]))

print m_arr
plt.plot(r_arr, m_arr)
plt.show()
