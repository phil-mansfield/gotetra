from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import colossus.Cosmology as cosmology
import scipy.stats as stats
import scipy.interpolate as intr
import numpy.random as rand
import scipy.special as sp

dir = "correction_data"
sub1 = ["h500_rad_sub1.dat", "h250_rad_sub1.dat",
        "h125_rad_sub1.dat", "h63_rad_sub1.dat"]
sub2 = ["h500_rad_sub2.dat", "h250_rad_sub2.dat",
        "h125_rad_sub2.dat", "h63_rad_sub2.dat"]
sub4 = ["h500_rad_sub4.dat", "h250_rad_sub4.dat",
        "h125_rad_sub4.dat", "h63_rad_sub4.dat"]

Ls = [500, 250, 125, 62.5]
cs = ["b", "r", "g", "m"]
sub1 = [path.join(dir, f) for f in sub1]
sub2 = [path.join(dir, f) for f in sub2]
sub4 = [path.join(dir, f) for f in sub4]

params = {"flat":True, "H0":70, "Om0":0.27,
          "Ob0":0.0469, "sigma8":0.82, "ns":0.95}
cosmo = cosmology.setCosmology("meowCosmo", params)

def r_to_m(r): return cosmo.rho_m(0) * 1e9 * 200 * r**3 * 4 * np.pi / 3

alpha = 0.5

low_n = 5e4
converge_n = 1e6

for i in xrange(len(Ls)):
    r_sp1, r_200m, gamma = map(
        np.array, zip(*np.loadtxt(sub1[i], usecols=(3, 6, 8))),
    )
    r_sp2 = np.loadtxt(sub2[i], usecols=(3,))
    r_sp4 = np.loadtxt(sub4[i], usecols=(3,))
    m_200m = r_to_m(r_200m)

    mp = (cosmo.rho_m(0) * Ls[i]**3 * 1e9) / 1024**3

    rat_1_2 = (r_sp1 - r_sp2) / r_sp2
    rat_2_4 = (r_sp2 - r_sp4) / r_sp4
    
    n1 = m_200m / mp
    n2 = m_200m / mp / 8
    n4 = m_200m / mp / 64
    mask1 = n1 > low_n
    mask2 = n2 > low_n
    mask4 = n4 > low_n

    range2 = (4 - np.log10(8), 8 - np.log10(8))
    range4 = (4 - np.log10(64), 8 -np.log10(64))

    def relation(n, in_val, func, range):
        if np.sum(n > 1e3) <= 0:
            return None, None, None, False

        vals, edges, _ = stats.binned_statistic(
            np.log10(n), in_val, func, bins=8, range=range,
        )
        counts, _, _ = stats.binned_statistic(
            np.log10(n), in_val, "count", bins=8, range=range,
        )

        mids = 10**((edges[1:] + edges[:-1]) / 2)
        count_lim = 5
        if np.sum(counts > count_lim) <= 1:
            return None, None, None, False

        i1d = intr.interp1d(
            np.log10(mids[counts > count_lim]),
            vals[counts > count_lim], kind="linear",
        )

        i_func = lambda xs: i1d(np.log10(xs))
        return mids[counts > count_lim], vals[counts > count_lim], i_func, True

    def bootstrap(val, func):
        if len(val) == 0: return np.nan
        meds = []
        for i in xrange(1000):
            idxs = rand.choice(
                np.arange(len(val)), size=(len(val)), replace=True,
            )
            meds.append(np.median(val[idxs]))
        return func(meds)

    def contour_plot(ns, rats, range):
        _, los, lo_i, _ = relation(
            ns, rats, lambda n: np.percentile(n, 16), range,
        )
        _, his, hi_i, _ = relation(
            ns, rats, lambda n: np.percentile(n, 84), range,
        )
        ns, meds, med_i, ok = relation(ns, rats, np.median, range)
        if not ok: return
        plt.fill_between(ns, los, his, alpha=0.3, facecolor=cs[i])
        plt.plot(ns, meds, lw=3, c=cs[i], label="L%g" % Ls[i])
        plt.xscale("log")
        plt.xlim(5e2, 1e7)
        plt.ylim(-0.05, 0.15)
        plt.plot([5e2, 1e7], [0, 0], "--k", lw=2)
        plt.legend()
        plt.xlabel(r"$N_{\rm 200m}$")


    def bootstrap_err_plot(ns, rats, range):
        _, b_his, b_hi_i, ok = relation(
            ns, rats, lambda x: bootstrap(
                x, lambda y: np.percentile(y, [50+68/2]),
            ), range)
        if not ok: return
        _, b_los, b_lo_i, ok = relation(
            ns, rats, lambda x: bootstrap(
                x, lambda y: np.percentile(y, [50-68/2]),
            ), range)
        if not ok: return
        ns, b_meds, b_med_i, ok = relation(
            ns, rats, lambda x: bootstrap(x, np.median), range,
        )
        if not ok: return
    
        plt.fill_between(ns, b_los, b_his, facecolor=cs[i], alpha=0.3)
        plt.plot(ns, b_meds, c=cs[i], lw=3, label="L%g" % Ls[i])
        plt.xscale("log")
        plt.ylim(-0.05, 0.15)
        plt.xlim(5e2, 1e7)
        plt.plot([5e2, 1e7], [0, 0], "--k", lw=2)


    plt.figure(0)
    contour_plot(n2, rat_1_2, range2)
    plt.ylabel(r"$(R_{\rm sp}(N_0) - R_{\rm sp}(N_1)) / R_{\rm sp}(N_1)$")
    plt.xlabel(r"$N_2$")
    plt.plot(n2, rat_1_2, ".", alpha=0.5, c=cs[i])
    plt.title("Scatter Around Medians: One subsampling")

    plt.figure(2)
    contour_plot(n4, rat_2_4, range4)
    plt.ylabel(r"$(R_{\rm sp}(N_1) - R_{\rm sp}(N_2)) / R_{\rm sp}(N_2)$")
    plt.xlabel("$N_3$")
    plt.plot(n4, rat_2_4, ".", alpha=0.5, c=cs[i])
    plt.title("Scatter Around Medians: Two subsamplings")

    plt.figure(3)
    bootstrap_err_plot(n2, rat_1_2, range2)
    plt.ylabel(r"$(R_{\rm sp}(N_0) - R_{\rm sp}(N_1)) / R_{\rm sp}(N_1)$")
    plt.xlabel(r"$N_2$")
    plt.title("Error on Medians: One subsampling")

    exp = 3 - i

    lxs = np.linspace(2.5, 7, 100)
    y0 = 0.05
    x0 = 4
    m =  -0.02
    a = 0.012
    corr = y0 + m * (lxs - x0) + a*exp
    x_cut = 5e3

    plt.plot(10**lxs, corr, c=cs[i])
    plt.plot([x_cut, x_cut], [-0.05, 0.15], "r--")

    plt.figure(5)
    lx_conv = -(y0 + a*exp) / m + x0
    print lx_conv 
    Fs = []
    def F(lx):
        mult = 1.0
        while lx < lx_conv:
            mult *= 1 + y0 + a*exp + m*(lx - x0)
            lx += np.log10(8)
        return mult
    for lx in lxs: Fs.append(F(lx))
    plt.plot(10**lxs, np.array(Fs) - 1, c=cs[i], lw=3)
    plt.xscale("log")
    plt.xlim(5e2, 1e7)
    plt.ylabel(r"$(R_{\rm sp}(N_{\rm converge}) - R_{\rm sp}(N))/R_{\rm sp}(N)$")
    plt.xlabel("$N$")
    fit = sp.gamma(np.maximum((lx_conv - lxs), 1))
    print fit

    plt.plot(10**lxs, fit, c=cs[i])

    plt.figure(6)
    plt.plot(mp * 10**lxs, np.array(Fs) - 1, c=cs[i], lw=3)
    plt.xscale("log")
    #plt.xlim(5e2, 1e7)
    plt.ylabel(r"$(R_{\rm sp}(N_{\rm converge}) - R_{\rm sp}(N))/R_{\rm sp}(N)$")
    plt.xlabel(r"$M_{\rm 200m}$")

    plt.figure(4)
    bootstrap_err_plot(n4, rat_2_4, range4)
    plt.ylabel(r"$(R_{\rm sp}(N_1) - R_{\rm sp}(N_2)) / R_{\rm sp}(N_2)$")
    plt.xlabel("$N_3$")
    plt.title("Error on Medians: Two subsamplings")

    if i == 0: continue
    corr = y0 + m * (lxs - x0) + a*(exp + 1)
    plt.plot([x_cut, x_cut], [-0.05, 0.15], "r--")
    plt.plot(10**lxs, corr, c=cs[i])

plt.show()
