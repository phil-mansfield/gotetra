from __future__ import division

import matplotlib.pyplot as plt
import numpy as np
import os.path as path
import colossus.Cosmology as cosmology
import scipy.stats as stats

dir = "data"
sub1 = ["h2000_rad.dat", "h500_rad.dat", "h250_rad.dat",
        "h125_rad.dat", "h63_rad.dat", "h63_rad_100_sub1.dat"]
sub2 = ["h2000_rad_sub2.dat", "h500_rad_sub2.dat", "h250_rad_sub2.dat",
        "h125_rad_sub2.dat", "h63_rad_sub2.dat", "h63_rad_100_sub2.dat"]
sub4 = ["h2000_rad_sub4.dat", "h500_rad_sub4.dat", "h250_rad_sub4.dat",
        "h125_rad_sub4.dat", "h63_rad_sub4.dat", "h63_rad_100_sub4.dat"]
Ls = [2000, 500, 250, 125, 62.5, 62.5]
cs = ["k", "b", "r", "g", "y", "m"]
sub1 = [path.join(dir, f) for f in sub1]
sub2 = [path.join(dir, f) for f in sub2]
sub4 = [path.join(dir, f) for f in sub4]

params = {"flat":True, "H0":70, "Om0":0.27,
          "Ob0":0.0469, "sigma8":0.82, "ns":0.95}
cosmo = cosmology.setCosmology("meowCosmo", params)

def r_to_m(r): return cosmo.rho_m(0) * 200 * r**3 * 4 * np.pi / 3

alpha = 0.5
for i in xrange(len(Ls)):
    r_sp_idx = 3
    r_sp1, r_200m, gamma = map(
        np.array, zip(*np.loadtxt(sub1[i], usecols=(r_sp_idx, 6, 8))),
    )
    r_sp2 = np.loadtxt(sub2[i], usecols=(r_sp_idx,))
    r_sp4 = np.loadtxt(sub4[i], usecols=(r_sp_idx,))
    m_200m = r_to_m(r_200m)

    mp = (cosmo.rho_m(0) * Ls[i]**3) / 1024**3
    n = m_200m / mp
    mask = n > 5e4

    plt.figure(0)
    rat_1_2 = (r_sp1 - r_sp2) / r_sp2
    rat_2_4 = (r_sp2 - r_sp4) / r_sp4
    plt.plot(n/8, rat_1_2, ".", alpha=alpha, c=cs[i], label="L%g" % Ls[i])

    meds, edges, _ = stats.binned_statistic(np.log10(n), rat_1_2, "median")
    ns, _, _ = stats.binned_statistic(np.log10(n), rat_1_2, "count")
    mids = 10**((edges[1:] + edges[:-1]) / 2)
    plt.plot(mids[ns > 0]/8, meds[ns > 0], c=cs[i], lw=3)

    plt.xlabel(r"$N_{\rm 200m}$")
    plt.ylabel(r"$(R_{\rm sp}(N) - R_{\rm sp}(N\times 8)) / R_{\rm sp}(N)$")
    plt.xscale("log")

    plt.figure(5)
    los, edges, _ = stats.binned_statistic(
        np.log10(n), rat_1_2, lambda xs: np.percentile(xs, 50-68/2),
    )
    his, edges, _ = stats.binned_statistic(
        np.log10(n), rat_1_2, lambda xs: np.percentile(xs, 50+68/2),
    )
    plt.fill_between(mids[ns > 0] / 8, his[ns > 0], los[ns > 0],
                     facecolor=cs[i], alpha=0.3)
    plt.plot(mids[ns > 0]/8, meds[ns > 0], lw=3, label="L%g" % Ls[i], c=cs[i])

    plt.xlabel(r"$N_{\rm 200m}$")
    plt.ylabel(r"$(R_{\rm sp}(N) - R_{\rm sp}(N\times 8)) / R_{\rm sp}(N)$")
    plt.xscale("log")
    plt.legend()

    plt.figure(1)
    plt.plot(gamma[n > 5e4], rat_1_2[n > 5e4], ".",
             alpha=alpha, c=cs[i], label="L%g" % Ls[i])
    plt.xlabel(r"$\Gamma$")
    plt.ylabel(r"$(R_{\rm sp, 1} - R_{\rm sp, 2}) / R_{\rm sp, 1}$")
    if np.sum(n > 5e4) <= 10: continue
    meds, edges, _ = stats.binned_statistic(
        gamma[n > 5e4], rat_1_2[n > 5e4], "median",
    )
    ns, _, _ = stats.binned_statistic(gamma[n > 5e4], rat_1_2[n > 5e4], "count")
    mids = (edges[1:] + edges[:-1]) / 2
    plt.plot(mids[ns > 0], meds[ns > 0], c=cs[i], lw=3)

    plt.figure(2)
    plt.plot(rat_2_4[mask], rat_1_2[mask], ".", alpha=alpha,
             c=cs[i], label="L%g" % Ls[i])

    range_mask = (rat_2_4 < 0.3) & (rat_2_4 > -0.3)
    if np.sum(mask) <= 0: continue
    meds, edges, _ = stats.binned_statistic(
        rat_2_4[mask & range_mask], rat_1_2[mask & range_mask], "median",
    )
    ns, _, _ = stats.binned_statistic(rat_2_4, rat_1_2, "count")
    mids = (edges[1:] + edges[:-1]) / 2
    plt.plot(mids[ns > 0], meds[ns > 0], c=cs[i], lw=3)
    plt.xlabel(r"$(R_{\rm sp, 1} - R_{\rm sp, 2}) / R_{\rm sp, 1}$")
    plt.ylabel(r"$(R_{\rm sp, 2} - R_{\rm sp, 4}) / R_{\rm sp, 2}$")
    

plt.figure(0)
plt.legend()
plt.ylim(-0.2, 0.3)
lo, hi = plt.xlim()
plt.plot([lo, hi], [0, 0], "--k", lw=2)
plt.xlim(lo, hi)

plt.figure(1)
plt.xlim(0, 7)
plt.ylim(-0.2, 0.3)
lo, hi = plt.xlim()
plt.plot([lo, hi], [0, 0], "--k", lw=2)
plt.xlim(lo, hi)
plt.legend()
plt.figure(2)
plt.legend(loc="lower left")
plt.ylim(-0.2, 0.3)
plt.xlim(-0.2, 0.3)
lo, hi = plt.xlim()
plt.plot([lo, hi], [0, 0], "--k", lw=2)
lo, hi = plt.ylim()
plt.plot([0, 0], [lo, hi], "--k", lw=2)
plt.xlim(lo, hi)
plt.ylim(lo, hi)
plt.plot([-0.2, 0.3], [-0.2 / 2, 0.3 / 2], "k")

edges = [2**11, 2**14, 2**17, 2**20]
labels = [r"$N_{\rm 200m} \in [2^{11}, 2^{14})$",
          r"$N_{\rm 200m} \in [2^{14}, 2^{17})$",
          r"$N_{\rm 200m} \in [2^{17}, 2^{20})$",
          r"$N_{\rm 200m} \geq 2^{20}$"]
rats = [[] for _ in edges]
gammas = [[] for _ in edges]
for i in xrange(len(Ls)):
    r_sp, r_200m, gamma = map(
        np.array, zip(*np.loadtxt(sub2[i], usecols=(3, 6, 8))),
    )
    m_200m = r_to_m(r_200m)
    mp = (cosmo.rho_m(0) * Ls[i]**3) / 1024**3
    n = m_200m / mp
    idxs = np.digitize(n, edges)

    for i in xrange(len(edges)):
        rats[i] = np.append(rats[i], (r_sp / r_200m)[idxs == i])
        gammas[i] = np.append(gammas[i], gamma[idxs == i])

cs = ["g", "r", "b", "k"]
alpha=0.5
rats = map(np.array, rats)
gammas = map(np.array, gammas)

bins = 10
for i in xrange(len(edges)):
    plt.figure(3)
    meds, edges, _ = stats.binned_statistic(
        gammas[i], rats[i], "median", range=(0,7), bins=bins,
    )
    ns, _, _ = stats.binned_statistic(
        gammas[i], rats[i], "count", range=(0, 7), bins=bins,
    )
    mids = (edges[1:] + edges[:-1]) / 2
    plt.plot(mids[ns>1], meds[ns>1], c=cs[i], lw=3, label=labels[i])

    ps1, _, _ = stats.binned_statistic(
        gammas[i], rats[i], lambda xs: np.percentile(xs, 50 + 68/2),
        range=(0,7), bins=bins,
    )
    ms1, _, _ = stats.binned_statistic(
        gammas[i], rats[i], lambda xs: np.percentile(xs, 50 - 68/2),
        range=(0,7), bins=bins,
    )

    plt.fill_between(mids[ns>1], ms1[ns>1], ps1[ns>1], color=cs[i], alpha=0.3)

    plt.figure(4)
    mult = 1.07**(len(cs) - (i+1) + 1)
    plt.plot(mids[ns>1], meds[ns>1]*mult, c=cs[i], lw=3, label=labels[i])
    plt.fill_between(
        mids[ns>1], ms1[ns>1]*mult, ps1[ns>1]*mult, color=cs[i], alpha=0.3,
    )
    #plt.plot(mids[ns>1], ms2[ns>1]*mult, c=cs[i])
    #plt.plot(mids[ns>1], ps2[ns>1]*mult, c=cs[i])

plt.figure(3)
plt.xlim(0, 7)
plt.ylim(0.6, 1.8)
plt.legend()
plt.ylabel(r"$R_{\rm sp} / R_{\rm 200m}$")
plt.xlabel(r"$\Gamma$")

plt.figure(4)
plt.xlim(0, 7)
plt.ylim(0.6, 1.8)
plt.legend()
plt.ylabel(r"$R_{\rm sp} / R_{\rm 200m}$")
plt.xlabel(r"$\Gamma$")

plt.show()
