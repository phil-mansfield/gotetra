from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import colossus.Cosmology as cosmology
import scipy.stats as stats

dir = "data"
fs = ["h63_rad.dat", "h125_rad.dat", "h250_rad.dat",
      "h500_rad.dat"]
fs = [path.join(dir, f) for f in fs]
Ls = [62.5, 125, 250, 500]

lo_exps = [10, 10, 10, 10]
hi_exps = [19, 20, 18, 16]

params = {"flat":True, "H0":70, "Om0":0.27,
          "Ob0":0.0469, "sigma8":0.82, "ns":0.95}
cosmo = cosmology.setCosmology("meowCosmo", params)

cs = ["m", "g", "r", "b"]

for i in xrange(len(Ls)):
    #if i == 0: continue
    mp = cosmo.rho_m(0) * 1e9 * Ls[i]**3 / 1024**3
    rows = zip(*np.loadtxt(fs[i], usecols=(2, 3, 4, 5, 6, 7, 8)))
    m_sp, r_sp, r_min, r_max, r_200m, m_200c, gamma = map(np.array, rows)

    plt.figure(0)
    plt.plot(gamma, r_sp / r_200m, ".", c=cs[i], alpha=0.5)

for i in xrange(len(Ls)):
    #if i == 0: continue
    mp = cosmo.rho_m(0) * 1e9 * Ls[i]**3 / 1024**3
    rows = zip(*np.loadtxt(fs[i], usecols=(2, 3, 4, 5, 6, 7, 8)))
    m_sp, r_sp, r_min, r_max, r_200m, m_200c, gamma = map(np.array, rows)

    plt.figure(0)
    meds, edges, _ = stats.binned_statistic(
        gamma, r_sp / r_200m, "median", range=(0, 7),
    )
    plt.plot((edges[1:] + edges[:-1])/2, meds, lw=3, c=cs[i])
"""
n_L63_range = range(10, 24)

def plot_mass_range(L63_n_min, L63_n_max):
    plt.figure()
    plt.xlim(0, 7)

    m_min = cosmo.rho_m(0) * (62.5/1024)**3 * L63_n_min * 1e9
    m_max = cosmo.rho_m(0) * (62.5/1024)**3 * L63_n_max * 1e9
    
    print(m_min, m_max)

    for i in xrange(len(Ls)):
        rows = zip(*np.loadtxt(fs[i], usecols=(3, 6, 8)))
        r_sp, r_200m, gamma = map(np.array, rows)

        mp = cosmo.rho_m(0) * (Ls[i]/1024)**3 * 1e9

        m_200m = (4 * np.pi / 3) * r_200m**3 * cosmo.rho_m(0) * 200e9
        mask = (m_200m <= m_max*8**i) & (m_200m >= m_min*8**i)
        if np.sum(mask) < 10: continue

        gamma = gamma[mask]
        r_sp = r_sp[mask]
        r_200m = r_200m[mask]

        plt.plot(gamma, r_sp/r_200m, ".", c=cs[i], alpha=0.5)
        meds, edges, _ = stats.binned_statistic(
            gamma, r_sp/r_200m, "median", range=(0, min(np.max(gamma), 7)),
            
        )
        plt.plot((edges[1:] + edges[:-1])/2, meds, lw=3, c=cs[i],
                 label=(r"L%g ($2^{%d}$ - $2^{%d}$)" %
                        (Ls[i], (np.log2(L63_n_min) - 3*i),
                         (np.log2(L63_n_max) - 3*i))))

    plt.xlim(0, 7)
    plt.ylim()
    plt.title(r"$M_{\rm 200m} \in [%.3g,\ %.3g]$" % (m_min, m_max))
    plt.legend(loc="upper right")

for i in xrange(10, 26, 2):
    plot_mass_range(2**i, 2**(i+2))
"""
plt.show()
