from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import scipy.signal as signal

dir = "data"
subs = ["sub1", "sub2", "sub4"]

med_base = path.join(dir, "h1000_med_%s.dat")
vol_base = path.join(dir, "h1000_vol_%s.dat")
med_prof_base = path.join(dir, "h1000_med_prof_%s.dat")
mean_prof_base = path.join(dir, "h1000_mean_prof_%s.dat")
sph_prof_base = path.join(dir, "h1000_sph_prof_%s.dat")
tet_prof_base = path.join(dir, "h1000_tet_prof_%s.dat")

r_sp_med_1, r200m, gamma = map(np.array, zip(*np.loadtxt(med_base % subs[0], usecols=(3, 6, 8))))
r_sp_med_2 = np.loadtxt(med_base % subs[1], usecols=(3,))
r_sp_med_4 = np.loadtxt(med_base % subs[2], usecols=(3,))

"""
plt.figure()
plt.title("Radius from Shell Median")
plt.plot(r_sp_med_2, r_sp_med_1 / r_sp_med_2 - 1, "r.", label=r"(R_2 - R_1) / R_2")
plt.plot(r_sp_med_4, r_sp_med_2 / r_sp_med_4 - 1, "b.", label=r"(R_4 - R_2) / R_4")
plt.plot([0.15, 0.5], [0, 0], "k")

plt.xlim(0.15, 0.5)
plt.ylim(-0.15, +0.15)

plt.figure()
plt.title("Radius from Shell Volume")
"""
r_sp_vol_1 = np.loadtxt(vol_base % subs[0], usecols=(3,))
r_sp_vol_2 = np.loadtxt(vol_base % subs[1], usecols=(3,))
r_sp_vol_4 = np.loadtxt(vol_base % subs[2], usecols=(3,))

"""
plt.plot(r_sp_vol_2, r_sp_vol_1 / r_sp_vol_2 - 1, "r.", label=r"(R_2 - R_1) / R_2")
plt.plot(r_sp_vol_4, r_sp_vol_2 / r_sp_vol_4 - 1, "b.", label=r"(R_4 - R_2) / R_4")
plt.plot([0.15, 0.5], [0, 0], "k")

plt.xlim(0.15, 0.5)
plt.ylim(-0.15, +0.15)
"""
def convert_rows(rows):
    n = (len(rows[0]) - 2) // 2
    rs = [row[2:n+2] for row in rows]
    rhos = [row[n+2:] for row in rows]
    return rs, rhos

prof_rs, med_prof_rhos_1 = convert_rows(np.loadtxt(med_prof_base % subs[0]))
_, med_prof_rhos_2 = convert_rows(np.loadtxt(med_prof_base % subs[1]))
_, med_prof_rhos_4 = convert_rows(np.loadtxt(med_prof_base % subs[2]))
_, mean_prof_rhos_1 = convert_rows(np.loadtxt(mean_prof_base % subs[0]))
_, mean_prof_rhos_2 = convert_rows(np.loadtxt(mean_prof_base % subs[1]))
_, mean_prof_rhos_4 = convert_rows(np.loadtxt(mean_prof_base % subs[2]))
sph_prof_rs, sph_prof_rhos_1 = convert_rows(np.loadtxt(sph_prof_base % subs[0]))
#_, sph_prof_rhos_2 = convert_rows(np.loadtxt(sph_prof_base % subs[1]))
#_, sph_prof_rhos_4 = convert_rows(np.loadtxt(sph_prof_base % subs[2]))
_, tet_prof_rhos_1 = convert_rows(np.loadtxt(tet_prof_base % subs[0]))
_, tet_prof_rhos_2 = convert_rows(np.loadtxt(tet_prof_base % subs[1]))
_, tet_prof_rhos_4 = convert_rows(np.loadtxt(tet_prof_base % subs[2]))

for i in [18, 24, 28]:
    plt.figure(i + 1)
    plt.title((r"Profile Comparison: $f_{12} = %.2g$, $f_{24} = %.2g$ " + 
               r"$R_{\rm sp} = %.3g$") %
              ((r_sp_vol_1[i] - r_sp_vol_2[i]) / r_sp_vol_2[i],
               (r_sp_vol_2[i] - r_sp_vol_4[i]) / r_sp_vol_4[i],
               r_sp_vol_1[i] / r200m[i]))
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\rho / \rho_{\rm m}$")
    plt.ylabel(r"$r\ [{\rm Mpc}/h]$")

    plt.plot(prof_rs[i], med_prof_rhos_1[i], "r", lw=2)
    plt.plot(prof_rs[i], med_prof_rhos_2[i], "g", lw=2)
    plt.plot(prof_rs[i], med_prof_rhos_4[i], "b", lw=2)

    plt.plot(prof_rs[i], 4*mean_prof_rhos_1[i], "r", lw=2)
    plt.plot(prof_rs[i], 4*mean_prof_rhos_2[i], "g", lw=2)
    plt.plot(prof_rs[i], 4*mean_prof_rhos_4[i], "b", lw=2)

    plt.plot(sph_prof_rs[i], sph_prof_rhos_1[i], "k", lw=2)
    plt.plot(sph_prof_rs[i], 4*sph_prof_rhos_1[i], "k", lw=2)
    plt.plot(sph_prof_rs[i], 16*sph_prof_rhos_1[i], "k", lw=2)
    #plt.plot(prof_rs[i], 16*sph_prof_rhos_2[i], "g", lw=3)
    #plt.plot(prof_rs[i], 16*sph_prof_rhos_4[i], "b", lw=3)

    plt.plot(prof_rs[i], 16*tet_prof_rhos_1[i], "r", lw=2)
    plt.plot(prof_rs[i], 16*8*tet_prof_rhos_2[i], "g", lw=2)
    plt.plot(prof_rs[i], 16*tet_prof_rhos_4[i], "b", lw=2)

    #lo, hi = plt.ylim()
    lo, hi = 10**-1, 10**4.5
    plt.plot([r200m[i], r200m[i]], [lo, hi], "k", label=r"$R_{\rm 200m}$")
    plt.plot([r_sp_vol_1[i], r_sp_vol_1[i]], [lo, hi],
             "r", label=r"$R_{\rm sp, 1}$")
    plt.plot([r_sp_vol_2[i], r_sp_vol_2[i]], [lo, hi],
             "g", label=r"$R_{\rm sp, 2}$")
    plt.plot([r_sp_vol_4[i], r_sp_vol_4[i]], [lo, hi],
             "b", label=r"$R_{\rm sp, 4}$")
    plt.legend()
    plt.ylim(lo, hi)
    """
    plt.figure()
    dlr = np.log10(prof_rs[i][1]) - np.log10(prof_rs[i][0])
    w = 61
    drhos_1 = signal.savgol_filter(
        np.log10(med_prof_rhos_1[i]), w,
        polyorder=4, deriv=1, delta=dlr,
    )
    drhos_2 = signal.savgol_filter(
        np.log10(med_prof_rhos_2[i]), w,
        polyorder=4, deriv=1, delta=dlr,
    )
    drhos_4 = signal.savgol_filter(
        np.log10(med_prof_rhos_4[i]), w,
        polyorder=4, deriv=1, delta=dlr,
    )
    plt.plot(prof_rs[i], drhos_1, "r", lw=3)
    plt.plot(prof_rs[i], drhos_2, "g", lw=3)
    plt.plot(prof_rs[i], drhos_4, "b", lw=3)

    lo, hi = plt.ylim()
    lo, hi = -8, 0
    plt.plot([r200m[i], r200m[i]], [lo, hi], "k", label=r"$R_{\rm 200m}$")
    plt.plot([r_sp_vol_1[i], r_sp_vol_1[i]], [lo, hi],
             "r", label=r"$R_{\rm sp, 1}$")
    plt.plot([r_sp_vol_2[i], r_sp_vol_2[i]], [lo, hi],
             "g", label=r"$R_{\rm sp, 2}$")
    plt.plot([r_sp_vol_4[i], r_sp_vol_4[i]], [lo, hi],
             "b", label=r"$R_{\rm sp, 4}$")
    plt.xscale("log")
    plt.legend(loc="lower right")
    plt.ylim(lo, hi)
    """
plt.show()
