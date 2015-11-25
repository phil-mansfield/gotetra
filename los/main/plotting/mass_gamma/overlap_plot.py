from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys
import numpy.random as rand
import os.path as path
import time

import colossus.Cosmology as cosmology

m_low, m_high = 1e11, 1e16
hist_bins = 20

Ls = [63, 125, 250, 500, 2000]

lo_exps = [10, 10, 10, 10, 10]
# (Corresponds to: upper_lims = 10**np.array([13.25, 14.4, 14.6, 15.1, 15.6]))
hi_exps = [19, 20, 18, 16, 12]

dir = "N_all"

cs = ["y", "g", "r", "b", "k"]
fs = ["L63_mass_gamma.dat", "L125_mass_gamma.dat", "L250_mass_gamma.dat",
      "L500_mass_gamma.dat", "L2000_mass_gamma.dat"]
fs = [path.join(dir, f) for f in fs]

params = {"flat":True, "H0":70, "Om0":0.27,
          "Ob0":0.0469, "sigma8":0.82, "ns":0.95}
cosmo = cosmology.setCosmology("meowCosmo", params)

def read(fname):
    ids, snaps, Gs, Ms = map(
        lambda xs: np.array(xs), zip(*np.loadtxt(fname, usecols=(0,1,2,3))),
    )
    ids = np.array(ids, dtype=int)
    snaps = np.array(snaps, dtype=int)
    return ids, snaps, Gs, Ms

def mp(L, N):
    return L**3 * cosmo.rho_m(0) / 1e9 * N**3

for i in xrange(len(Ls)):
    print("L%d" % int(Ls[i]))
    t0 = time.time()
    ids, snaps, gs, ms = read(fs[i])
    t1 = time.time()

    ns = ms / mp(Ls[i], 1024)
    lo, hi = lo_exps[i], hi_exps[i]

    ir_idx = (ns >= 2**lo) & (ns <= 2**hi)
    #print(len(ms), np.sum(ir_idx))
    #print(2**lo, 2**hi)
    #print(ms)
    #print(np.sum(ns >= 2**lo), np.sum(ns <= 2**hi))
    ns = ns[ir_idx]
    ms = ms[ir_idx]
    gs = gs[ir_idx]


    arange = np.arange(0, len(ms))
    n_lim = 25
    def filter_bin(m_low, m_high, g_low, g_high):
        mask = ((ms < m_high) & (ms > m_low) & 
                (gs < g_high) & (gs > g_low))
        idxs = arange[mask]
        if len(idxs) <= n_lim: return mask
        idxs = rand.choice(idxs, n_lim, replace=False)
        mask = np.zeros(len(mask), dtype=bool)
        mask[idxs] = True
        return mask

    mask = np.zeros(len(ms), dtype=bool)
    g_min, g_max, dg  = 0, 7, 1
    g_bins = int(round((g_max - g_min) / dg))
    for j in xrange(lo, hi):
        m_low, m_high = mp(Ls[i], 1024) * 2**j, mp(Ls[i], 1024) * 2**(j+1)
        for g_bin in xrange(g_bins):
            g_low, g_high = g_min + g_bin*dg, g_min + (g_bin + 1)*dg
            mask = mask | filter_bin(m_low, m_high, g_low, g_high)

    print("|full| = %d, |ir| = %d, |filtered| = %d" %
          (len(ir_idx), len(ns), np.sum(mask)))

    plt.figure(0)
    plt.hist(np.log2(ns[mask]), histtype="step", lw=3, range=(lo, hi),
             label=r"$L$ = %g" % Ls[i], color=cs[i])
    plt.hist(np.log2(ns), histtype="step", lw=1, range=(lo, hi), color=cs[i])

    plt.figure(1)
    plt.hist(np.log10(ms[mask]), histtype="step", lw=3, range=(10, 16),
             label=r"$L$ = %g" % Ls[i], color=cs[i])
    plt.hist(np.log10(ms), histtype="step", lw=1, range=(10, 16), color=cs[i])

    plt.figure(2)
    plt.hist(gs[mask], histtype="step", lw=3, range=(0, 7),
             label=r"$L$ = %g" % Ls[i], color=cs[i], bins=7)
    plt.hist(gs, histtype="step", lw=1, range=(0, 7),
             bins=7, color=cs[i])

    t2 = time.time()

    plt.figure(3)
    plt.plot(gs[mask], ms[mask], ".", c=cs[i])
    plt.title("L%d" % Ls[i])
    plt.ylabel(r"$M_{\rm 200m}$")
    plt.xlabel(r"$\Gamma$")
    plt.yscale("log")

    id_name = "out/L%d_overlap_id.dat" % Ls[i]

    ids, snaps = ids[mask], snaps[mask]
    with open(id_name, "w+") as fp:
        for i in xrange(len(ids)):
            print(ids[i], snaps[i], file=fp)

    print("t1 = %.3g t2 = %.3g" % (t1 - t0, t2 - t1))

plt.figure(0)
plt.xlabel(r"$\log_2n$")
plt.yscale("log")
plt.legend()

plt.figure(1)
plt.xlabel(r"$\log_{10}M_{\rm 200m}$")
plt.yscale("log")
plt.legend()

plt.figure(2)
plt.xlabel(r"$\log_{10}\Gamma$")
plt.yscale("log")
plt.legend()

plt.show()
