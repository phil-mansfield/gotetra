from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats

def c_to_gamma(c_vir):
    # Only valid for z = 0.
    return 3.43 - 2.74 * np.log10(c_vir)

def gamma_to_r_rat(gamma):
    return 0.54 * (1 + 0.53 * 0.27) * (1 + 1.36 * np.exp(-gamma / 3.04))

def gamma_to_m_rat(gamma):
    return 0.59 * (1 + 0.35 * 0.27) * (1 + 0.92 * np.exp(-gamma / 4.54))

fnames = sys.argv[1:]
g_m_sp, g_r_sp, g_m200, g_r200, g_rvir, g_rs, g_G1, g_G25, g_G5 = [], [], [], [], [], [], [], [], []
for fname in fnames:
    cols = map(list, zip(*np.loadtxt(fname)))
    m_sp, r_sp, m200, r200, rvir, rs, G1, G25, G5 = cols[2:]
    g_m_sp += m_sp
    g_r_sp += r_sp
    g_m200 += m200
    g_r200 += r200
    g_rvir += rvir
    g_rs += rs
    g_G1 += G1
    g_G25 += G25
    g_G5 += G5
    
g_m_sp = np.array(g_m_sp)
g_r_sp = np.array(g_r_sp)
g_m200 = np.array(g_m200)
g_r200 = np.array(g_r200)
g_rvir = np.array(g_rvir)
g_rs = np.array(g_rs)
g_G1 = np.array(g_G1)
g_G25 = np.array(g_G25)
g_G5 = np.array(g_G5)

plt.figure()
m_rat = g_m_sp / g_m200
mask = (m_rat < 2.2) & (m_rat > 0.8)
#mask = m_rat < 5
plt.plot(g_m200[mask], m_rat[mask], "b.")
plt.xscale("log")
plt.ylabel(r"$M_{\rm sp} / M_{200\rm m}$")
plt.xlabel(r"$M_{200\rm m}$ [$M_\odot/h$]")
plt.figure()
r_rat = g_r_sp / g_r200
plt.plot(g_r200[mask], r_rat[mask], "b.")
plt.ylabel(r"$R_{\rm sp} / R_{200\rm m}$")
plt.xlabel(r"$R_{200\rm m}$ [Mpc/$h$]")
plt.xscale("log")

c = g_rvir / g_rs
gamma = c_to_gamma(c)
gamma_low, gamma_high = 0, 6
gamma_ref = np.linspace(gamma_low, gamma_high, 100)

plt.figure()
#plt.plot(g_G1[mask], r_rat[mask], "b.")
#plt.plot(g_G25[mask], r_rat[mask], "r.")
plt.plot(g_G5[mask], r_rat[mask], "r.")
plt.plot(gamma_ref, gamma_to_r_rat(gamma_ref), "k", lw=3)
vals, edges, _ = stats.binned_statistic(
    g_G5[mask], r_rat[mask], "mean", 20, (gamma_low, gamma_high),
)
sqrs, edges, _ = stats.binned_statistic(
    g_G5[mask], r_rat[mask]**2, "mean", 20, (gamma_low, gamma_high),
)
std = np.sqrt(sqrs - vals**2)
plt.plot((edges[1:] + edges[:-1]) / 2, vals,  "r", lw=3)
plt.plot((edges[1:] + edges[:-1]) / 2, vals + std,  "r", lw=1)
plt.plot((edges[1:] + edges[:-1]) / 2, vals - std,  "r", lw=1)

plt.ylabel(r"$R_{\rm sp} / R_{200\rm m}$")
plt.xlabel(r"$\Gamma$")
plt.xlim(gamma_low, gamma_high)
plt.figure()
#plt.plot(g_G1[mask], m_rat[mask], "b.")
#plt.plot(g_G25[mask], m_rat[mask], "r.")
plt.plot(g_G5[mask], m_rat[mask], "r.")
plt.plot(gamma_ref, gamma_to_m_rat(gamma_ref), "k", lw=3)
vals, edges, _ = stats.binned_statistic(
    g_G5[mask], m_rat[mask], "mean", 20, (gamma_low, gamma_high),
)
sqrs, edges, _ = stats.binned_statistic(
    g_G5[mask], m_rat[mask]**2, "mean", 20, (gamma_low, gamma_high),
)
std = np.sqrt(sqrs - vals**2)
plt.plot((edges[1:] + edges[:-1]) / 2, vals,  "r", lw=3)
plt.plot((edges[1:] + edges[:-1]) / 2, vals + std,  "r", lw=1)
plt.plot((edges[1:] + edges[:-1]) / 2, vals - std,  "r", lw=1)
plt.ylabel(r"$M_{\rm sp} / M_{200\rm m}$")
plt.xlabel(r"$\Gamma$")
plt.xlim(gamma_low, gamma_high)

"""
plt.figure()
plt.title(r"$M_{200 \rm m} \approx 10^{11}\ M_\odot/h$")
plt.xlabel(r"$M_{\rm sp} / M_{200\rm m}$")
m_rat = g_m_sp / g_m200
m_min, m_max = 1, 2
plt.hist(m_rat[m_rat < 2], bins=10, range=(m_min, m_max))
plt.xlim(m_min, m_max)
plt.figure()
plt.title(r"$M_{200 \rm m} \approx 10^{11}\ M_\odot/h$")
plt.xlabel(r"$R_{\rm sp} / R_{200\rm m}$")
r_min, r_max = 1, 2
plt.hist(g_r_sp / g_r200, bins=10, range=(r_min, r_max))
plt.xlim(r_min, r_max)
"""
plt.show()
