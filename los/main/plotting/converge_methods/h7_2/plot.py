from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import scipy.signal as signal

import colossus.Cosmology as cosmology

params = {"flat":True, "H0":70, "Om0":0.27,
          "Ob0":0.0469, "sigma8":0.82, "ns":0.95}
cosmo = cosmology.setCosmology("meowCosmo", params)

dir = "data"
subs = ["sub1", "sub2", "sub4", "sub8"]

h_num = 7

sph_file = "data/h%d_sph_prof.dat" % h_num
tet_prof_base = path.join(dir, "h%d_tet_prof_%%s.dat" % h_num)
tril_prof_base = path.join(dir, "h%d_tril_prof_%%s.dat" % h_num)
cube_prof_base = path.join(dir, "h%d_cube_prof_%%s.dat" % h_num)
vol_base = path.join(dir, "h%d_vol_%%s.dat" % h_num)

def convert_rows(rows):
    n = (len(rows[0]) - 2) // 2
    rs = [row[2:n+2] for row in rows]
    rhos = [row[n+2:] for row in rows]
    return rs, rhos

def mass_prof(rs, rhos):
    prev_r = 0.0

    factor = cosmo.rho_m(0) * 1e9 * 4 * np.pi / 3
    outer_m = factor * rhos * rs**3
    inner_m = factor * rhos * np.append([0], rs[:-1]**3)
    dm = outer_m - inner_m

    return np.cumsum(dm)
        

sph_prof_rs, sph_prof_rhos = convert_rows(np.loadtxt(sph_file))

prof_rs, tet_prof_rhos_1 = convert_rows(np.loadtxt(tet_prof_base % subs[0]))
_, tet_prof_rhos_2 = convert_rows(np.loadtxt(tet_prof_base % subs[1]))
_, tet_prof_rhos_4 = convert_rows(np.loadtxt(tet_prof_base % subs[2]))
_, tet_prof_rhos_8 = convert_rows(np.loadtxt(tet_prof_base % subs[3]))

_, tril_prof_rhos_1 = convert_rows(np.loadtxt(tril_prof_base % subs[0]))
_, tril_prof_rhos_2 = convert_rows(np.loadtxt(tril_prof_base % subs[1]))
_, tril_prof_rhos_4 = convert_rows(np.loadtxt(tril_prof_base % subs[2]))
_, tril_prof_rhos_8 = convert_rows(np.loadtxt(tril_prof_base % subs[3]))


_, cube_prof_rhos_1 = convert_rows(np.loadtxt(cube_prof_base % subs[0]))
_, cube_prof_rhos_2 = convert_rows(np.loadtxt(cube_prof_base % subs[1]))
_, cube_prof_rhos_4 = convert_rows(np.loadtxt(cube_prof_base % subs[2]))
_, cube_prof_rhos_8 = convert_rows(np.loadtxt(cube_prof_base % subs[3]))

r_sp_1, r200m, gamma = map(np.array, zip(*np.loadtxt(vol_base % subs[0], usecols=(3, 6, 8))))
r_sp_2 = np.loadtxt(vol_base % subs[1], usecols=(3,))
r_sp_4 = np.loadtxt(vol_base % subs[2], usecols=(3,))
r_sp_8 = np.loadtxt(vol_base % subs[3], usecols=(3,))

for i in xrange(len(prof_rs)):
    m200m = 4 * np.pi / 3 * cosmo.rho_m(0) * 1e9 * r200m[i]**3 * 200
    f_12 = (r_sp_1[i] - r_sp_2[i]) / (r_sp_2[i])
    f_24 = (r_sp_2[i] - r_sp_4[i]) / (r_sp_4[i])

    plt.figure()

    plt.title((r"$M_{\rm 200m} = %.1g\ M_\odot/h$ $N_{\rm 200m} = %.0g$ " +
               r"$\Gamma = %.1f$ $f_{12}=%.2f$ $f_{24}=%.2f$") %
              (m200m, float(m200m/1.7e7), gamma[i], f_12, f_24))

    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel(r"$\rho / \rho_{\rm m}$")
    plt.xlabel(r"$r\ [{\rm Mpc}/h]$")
    
    plt.plot(prof_rs[i], tet_prof_rhos_8[i], ":m", lw=1, label="gotetra, sub-8")
    plt.plot(prof_rs[i], tet_prof_rhos_4[i], ":b", lw=1, label="gotetra, sub-4")
    plt.plot(prof_rs[i], tet_prof_rhos_2[i], ":g", lw=1, label="gotetra, sub-2")
    plt.plot(prof_rs[i], tet_prof_rhos_1[i], ":r", lw=1, label="gotetra, sub-1")

    plt.plot(prof_rs[i], tril_prof_rhos_8[i], "--m", lw=2)
    plt.plot(prof_rs[i], tril_prof_rhos_4[i], "--b", lw=2)
    plt.plot(prof_rs[i], tril_prof_rhos_2[i], "--g", lw=2)
    plt.plot(prof_rs[i], tril_prof_rhos_1[i], "--r", lw=2)

    plt.plot(prof_rs[i], cube_prof_rhos_8[i], "m", lw=2)
    plt.plot(prof_rs[i], cube_prof_rhos_4[i], "b", lw=2)
    plt.plot(prof_rs[i], cube_prof_rhos_2[i], "g", lw=2)
    plt.plot(prof_rs[i], cube_prof_rhos_1[i], "r", lw=2)

    plt.plot(sph_prof_rs[i], sph_prof_rhos[i], "k", lw=2)

    lo, hi = plt.ylim()
    plt.plot([r200m[i], r200m[i]], [lo, hi], "--k",
             lw=2, label=r"$R_{\rm 200m}$")
    plt.ylim(lo, hi)

    plt.plot([r_sp_1[i], r_sp_1[i]], [lo, hi], "r")
    plt.plot([r_sp_2[i], r_sp_2[i]], [lo, hi], "b")
    plt.plot([r_sp_4[i], r_sp_4[i]], [lo, hi], "g")
    plt.legend(loc="upper right")
    
    plt.figure()

    plt.plot(prof_rs[i], tet_prof_rhos_8[i] / sph_prof_rhos[i],
             ":m", lw=2, label="gotetra, sub-8")
    plt.plot(prof_rs[i], tet_prof_rhos_4[i] / sph_prof_rhos[i],
             ":b", lw=2, label="gotetra, sub-4")
    plt.plot(prof_rs[i], tet_prof_rhos_2[i] / sph_prof_rhos[i],
             ":g", lw=2, label="gotetra, sub-2")
    plt.plot(prof_rs[i], tet_prof_rhos_1[i] / sph_prof_rhos[i],
             ":r", lw=2, label="gotetra, sub-1")

    plt.plot(prof_rs[i], tril_prof_rhos_8[i] / sph_prof_rhos[i], "--m", lw=2)
    plt.plot(prof_rs[i], tril_prof_rhos_4[i] / sph_prof_rhos[i], "--b", lw=2)
    plt.plot(prof_rs[i], tril_prof_rhos_2[i] / sph_prof_rhos[i], "--g", lw=2)
    plt.plot(prof_rs[i], tril_prof_rhos_1[i] / sph_prof_rhos[i], "--r", lw=2)

    plt.plot(prof_rs[i], cube_prof_rhos_8[i] / sph_prof_rhos[i], "m", lw=2)
    plt.plot(prof_rs[i], cube_prof_rhos_4[i] / sph_prof_rhos[i], "b", lw=2)
    plt.plot(prof_rs[i], cube_prof_rhos_2[i] / sph_prof_rhos[i], "g", lw=2)
    plt.plot(prof_rs[i], cube_prof_rhos_1[i] / sph_prof_rhos[i], "r", lw=2)

    plt.xscale("log")
    plt.ylim(0.4, 1.6)
    r_min, r_max = plt.xlim()
    lo, hi = plt.ylim()
    plt.plot([r_min, r_max], [1, 1], "k", lw=2)
    plt.plot([r_sp_1[i], r_sp_1[i]], [lo, hi], "r")
    plt.plot([r_sp_2[i], r_sp_2[i]], [lo, hi], "b")
    plt.plot([r_sp_4[i], r_sp_4[i]], [lo, hi], "g")
    plt.ylim(lo, hi)

    plt.ylabel(r"$\rho / \rho_{\rm particle}$")
    plt.xlabel(r"$r\ [{\rm Mpc}/h]$")
    plt.title((r"$M_{\rm 200m} = %.1g\ M_\odot/h$ $N_{\rm 200m} = %.0g$ " +
               r"$\Gamma = %.1f$ $f_{12}=%.2f$ $f_{24}=%.2f$") %
              (m200m, float(m200m/1.7e7), gamma[i], f_12, f_24))

plt.show()
