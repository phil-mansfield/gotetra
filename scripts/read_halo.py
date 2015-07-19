from __future__ import division

import numpy as np
import sys

NUM_ROWS = 5000

H100 = 0.7
H70  = H100 / 0.7

H0_mks   = 3.24086e-18
G_mks    = 6.67259e-11
mpc_mks  = 3.08560e+22
msun_mks = 1.98900e+30

omega_M = 0.27
omega_L = 0.73
omega_B = 0.0469
omega_R = 0.0

def hubble_frac(z):
    return np.sqrt(omega_R*(1.0+z)**4.0 +
                   omega_M*(1.0+z)**3.0 +
                   omega_L)

def rho_c_mks(z):
    H = hubble_frac(z) * H0_mks * H100
    return 3.0 * H * H / (8.0 * np.pi * G_mks)

def rho_c(z):
    return rho_c_mks(z) * mpc_mks * mpc_mks * mpc_mks / msun_mks

# Repeat after me: X mpc/h, Y mpc -> X / h = Y.
def rho_c_h(z):
    mpc_mks_h, msun_mks_h = mpc_mks / H100, msun_mks / H100
    return rho_c_mks(z) * mpc_mks_h * mpc_mks_h * mpc_mks_h / msun_mks_h

def rho_m(z):
    return rho_c(z) * omega_M * (1+z)**3.0

def rho_m_h(z):
    return rho_c_h(z) * omega_M * (1+z)**3.0

def m_to_r(m, rho):
    v = m / rho
    return (v * 3 / 4 / np.pi)**(1/3)

file = sys.argv[1]
# Things I want:
# scale - 0
# id - 1
# R_vir - 11
# R_s - 12
# R_s_Klypin - 34
# M_200m - 36
# M_200c -  37
# M_500c - 38
# M_2500c - 39
# x - 17
# y - 18
# z - 19
#cols = np.loadtxt(file, usecols=(0,1,11,12,34,36,37,38,38,17,18,19))
cols = np.loadtxt(file, usecols=(37,17,18,19,11))
print "# Total halos in file: %d" % len(cols)
exit(0)
cols = ((cols[cols[:, 0].argsort()])[-NUM_ROWS:])[::-1]

class Struct(object): pass
ss = [Struct() for i in xrange(len(cols))]

for (s, col) in zip(ss, cols):
    #s.scale, s.id, s.R_vir, s.R_s = col[0], col[1], col[2]/1e3, col[3]/1e3
    #s.R_s_Klypin, s.M_200m, s.M_200c, s.M_500c = col[4]/1e3, col[5], col[6], col[7]
    #s.M_2500c, s.x, s.y, s.z = col[8], col[9], col[10], col[11]
    s.M_200c = col[0]
    s.x, s.y, s.z = col[1], col[2], col[3]
    s.R_vir = col[4] / 1000


print "# %4s %9s %9s %9s %9s %9s %9s" % ("N", "M200c", "R200c", "X", "Y", "Z", "Rvir")
#print ("# %4s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s" %
#       ("N", "X", "Y", "Z", "R_vir", "R_s", "R_s_Kly", "R_200m", "R_200c", "R_500c", "R_2500c"))
for (i, s) in enumerate(ss):
    #z = 1 / (1 + 1/s.scale)
    #rho_200m, rho_200c = 200 * rho_m_h(z), 200 * rho_c_h(z)
    #rho_500c, rho_2500c = 500 * rho_m_h(z), 2500 * rho_c_h(z)
    #print ("  %4d %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g" % 
    #       (i, s.x, s.y, s.z, 
    #        s.R_vir, s.R_s, s.R_s_Klypin,
    #        m_to_r(s.M_200m, rho_200m),
    #        m_to_r(s.M_200c, rho_200c),
    #        m_to_r(s.M_500c, rho_500c),
    #        m_to_r(s.M_2500c, rho_2500c)))
    rho_200c = 200 * rho_c_h(0)
    print ("  %4s %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g" %
           (i, s.M_200c, m_to_r(s.M_200c, rho_200c), s.x, s.y, s.z, s.R_vir))
