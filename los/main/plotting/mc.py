import numpy as np
import numpy.random as rand

a, b, c = 1, 30, 2

us = rand.random(100 * 1000)
vs = rand.random(100 * 1000)

ths = 2 * np.pi * us
phis = np.arccos(2 * vs - 1)

def r(th, phi):
    return np.sqrt(1 / (np.cos(th)**2*np.sin(phi)**2 / a**2 +
                        np.sin(th)**2*np.sin(phi)**2 / b**2 + 
                        np.cos(phi)**2 / c**2))
    

print 4 * np.pi / 3 * a * b * c
rs = r(ths, phis)
rs = rs**3
print 4 * np.pi / 3 * np.mean(rs)
