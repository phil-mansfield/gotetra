import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

cs = ["y", "g", "r", "b", "k"]
Ls = [62.5, 125, 250, 500, 2000]
fnames = ["vmax_63.dat", "vmax_125.dat", "vmax_250.dat", "vmax_500.dat", "vmax_2000.dat"]
lims = [(40, 50), (80, 100), (150, 170), (270, 320), (600, 650)]

plt.figure(0)
plt.ylim(1, 5e6)
#plt.figure(1)
#plt.ylim(1, 1e4)

frac = 0.1

for c, fname, lim, L in zip(cs, fnames, lims, Ls):
    ids, snaps, vmax, m200b, pid = map(np.array, zip(*np.loadtxt(fname)))

    n = int(len(ids) * frac)
    ids = ids[:n]
    snaps = snaps[:n]
    vmax = vmax[:n]
    m200b = m200b[:n]
    pid = pid[:n]

    vals, edges = np.histogram(np.log10(vmax[pid < 0]), bins=40)
    vals = np.cumsum(vals[::-1])[::-1]
    mids = 10**((edges[1:] + edges[:-1]) / 2)

    plt.figure(0)

    plt.plot(mids, vals, "%s" % c, lw=3, label="%g Mpc/$h$" % L)

    vals, edges = np.histogram(np.log10(vmax[pid > 0]), bins=40)
    vals = np.cumsum(vals[::-1])[::-1]
    mids = 10**((edges[1:] + edges[:-1]) / 2)
    plt.plot(mids, vals, "--%s" % c, lw=3)

    full_lim, sub_lim = lim

    ylo, yhi = plt.ylim()
    plt.plot([full_lim, full_lim], [ylo, yhi], c)
    plt.plot([sub_lim, sub_lim], [ylo, yhi], "--%s" % c)
    plt.ylim(ylo, yhi)

    plt.figure(1)
    meds, edges, _ = stats.binned_statistic(np.log10(vmax), m200b, "median")
    sig2ps, edges, _ = stats.binned_statistic(
        np.log10(vmax), m200b, lambda xs: np.percentile(xs, 2.5),
    )
    sig2ms, edges, _ = stats.binned_statistic(
        np.log10(vmax), m200b, lambda xs: np.percentile(xs, 97.5),
    )
    sig4ps, edges, _ = stats.binned_statistic(
        np.log10(vmax), m200b, lambda xs: np.percentile(xs, 0.01),
    )
    sig4ms, edges, _ = stats.binned_statistic(
        np.log10(vmax), m200b, lambda xs: np.percentile(xs, 99.99),
    )
    mids = 10**((edges[1:] + edges[:-1]) / 2)
    plt.plot(mids, meds, lw=3, c=c)
    plt.plot(mids, sig2ps, lw=1, c=c)
    plt.plot(mids, sig2ms, lw=1, c=c)
    #plt.plot(mids, sig4ps, lw=1, c=c)
    #plt.plot(mids, sig4ms, lw=1, c=c)

    #ylo, yhi = plt.ylim()
    #plt.plot([full_lim, full_lim], [ylo, yhi], c)
    #plt.plot([sub_lim, sub_lim], [ylo, yhi], "--%s" % c)

plt.figure(0)
plt.legend(loc="lower left")
plt.xscale("log")
plt.yscale("log")

plt.xlabel(r"$V_{\rm max}$")
plt.ylabel(r"$N(> V_{\rm max})$")

plt.figure(1)
plt.ylabel(r"$M_{\rm 200b}$")
plt.xlabel(r"$V_{\rm max}$")
plt.xscale("log")
plt.yscale("log")

plt.show()
