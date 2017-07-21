import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import scipy
from scipy.misc import derivative

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

# h_p = []
# for x in p:
#    h_p.append(scipy.integrate.quad(integrand, 0, x)[0])

h_p = scipy.integrate.cumtrapz(integrand(p), p, initial=0)
p_h = interp1d(h_p, p, fill_value="extrapolate")


# e_h1 = np.exp(fun(np.log(p_h(h_p))))
e_h = interp1d(h_p, e_p(p), fill_value="extrapolate")


def dydt(theta, h):
    m, r = theta
    r_prime = (-1) * r * (r - 2 * m) / (m + 4 * np.pi * r**3 * p_h(h))
    m_prime = 4 * np.pi * r**2 * e_h(h) * r_prime
#    m_prime = 4 * np.pi * r**2 * e_p(p_h(h)) * r_prime

    return [m_prime, r_prime]

h_c = 0.7818
h = np.arange(h_c, 0, -0.001)
sol = scipy.integrate.odeint(
    dydt, [2 * 10 ** (-31), 10**(-10)], h)

m_h = sol[:, 0]
r_h = sol[:, 1]

m_r = interp1d(r_h, m_h, fill_value="extrapolate")
h_r = interp1d(r_h, h, fill_value="extrapolate")

e_r = interp1d(r_h, e_h(h_r(r_h)), fill_value="extrapolate")
p_r = interp1d(r_h, p_h(h_r(r_h)), fill_value="extrapolate")

print p_r(2)
print 'starting GR'
def dedp(r):
	dx = 1e-6
	return (e_r(r+dx) - e_r(r-dx))/(p_r(r+dx)-p_r(r-dx))
#    derivative(e_r, r, dx=1e-6) / derivative(p_r, r, dx=1e-6)

M = 2.05823201401
R = 10.19147416 * 0.667

print dedp(0)

def dydt1(y, r):
    f = 1 - M / R
    #temp = y*(y-1)+2/f*(1-3*m_r(r)/r-2*np.pi*r**2*(e_r(r)+3*p_r(r)))*y-1/f*(6-4*np.pi*r**2*(e_r(r)+p_r(r)*(3+dedp(r))))
    temp = y*(y-1)+2/f*(1-3*m_r(r)/r-2*np.pi*r**2*(e_r(r)+3*p_r(r)))*y-1/f*(6-4*np.pi*r**2*(e_r(r)+p_r(r)*(3+dedp(r))))
    #temp = y * (y - 1) + 2 / f * (1 - 3 * m_r(r) / r - 2 * np.pi * r**2 * (e_r(r)) +3 * p_r(r)) - 1 / f * (6 - 4 * np.pi * r**2 * (e_r(r)) + p_r(r)) * (3 + dedp(r))
    temp = (-1) * temp / r
    return temp

sol = scipy.integrate.odeint(dydt1, 2, np.arange(0.0001, R, 0.001))
plt.plot(np.arange(0.0001, R, 0.001),sol)
plt.show()
