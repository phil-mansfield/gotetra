from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys
import numpy.random as rand
import os.path as path

m_low, m_high = 1e11, 1e16
hist_bins = 20

lims = [1e11, 10**12.5, 10**13.5, 10**14.25, 2e15, 1e16]
pms = [1.7e7, 1.4e8, 1.1e9, 8.8e9, 5.6e11]
Ls = [63, 125, 250, 500, 2000]

dir = "N10_000"

cs = ["y", "g", "r", "b", "k"]
fs = ["L63_mass_gamma.dat", "L125_mass_gamma.dat", "L250_mass_gamma.dat",
      "L500_mass_gamma.dat", "L2000_mass_gamma.dat"]
fs = [path.join(dir, f) for f in fs]

def read(fname):
    return map(lambda xs: np.array(xs),
               zip(*np.loadtxt(fname, usecols=(2,3))))

def read_ids(fname):
    return map(lambda xs: np.array(xs, dtype=int),
               zip(*np.loadtxt(fname, usecols=(0,1))))

def plot_hist(ms, c):
    plt.hist(ms, bins=hist_bins, color=c,
             range=(np.log10(m_low), np.log10(m_high)))

sig1, sig2 = 68, 95
ps = [50 - sig2/2, 50 - sig1/2, 50, 50 + sig1/2, 50 + sig2/2]

def plot_percentile_fill(ms, gs, c):
    sig2m, edges, _ = stats.binned_statistic(
        ms, gs, lambda xs: np.percentile(xs, ps[0]), bins=hist_bins,
        range=(np.log10(m_low), np.log10(m_high)),
    )
    sig1m, edges, _ = stats.binned_statistic(
        ms, gs, lambda xs: np.percentile(xs, ps[1]), bins=hist_bins,
        range=(np.log10(m_low), np.log10(m_high)),
    )
    med, edges, _ = stats.binned_statistic(
        ms, gs, lambda xs: np.percentile(xs, ps[2]), bins=hist_bins,
        range=(np.log10(m_low), np.log10(m_high)),
    )
    sig1p, edges, _ = stats.binned_statistic(
        ms, gs, lambda xs: np.percentile(xs, ps[3]), bins=hist_bins,
        range=(np.log10(m_low), np.log10(m_high)),
    )
    sig2p, edges, _ = stats.binned_statistic(
        ms, gs, lambda xs: np.percentile(xs, ps[4]), bins=hist_bins,
        range=(np.log10(m_low), np.log10(m_high)),
    )

    mids = (edges[1:] + edges[:-1]) / 2
    plt.fill_between(mids, sig2m, sig2p, facecolor=c, alpha=0.3)
    plt.fill_between(mids, sig1m, sig1p, facecolor=c, alpha=0.3)

def plot_percentile_line(ms, gs, c):
    sig2m, edges, _ = stats.binned_statistic(
        ms, gs, lambda xs: np.percentile(xs, ps[0]), bins=hist_bins,
        range=(np.log10(m_low), np.log10(m_high)),
    )
    sig1m, edges, _ = stats.binned_statistic(
        ms, gs, lambda xs: np.percentile(xs, ps[1]), bins=hist_bins,
        range=(np.log10(m_low), np.log10(m_high)),
    )
    med, edges, _ = stats.binned_statistic(
        ms, gs, lambda xs: np.percentile(xs, ps[2]), bins=hist_bins,
        range=(np.log10(m_low), np.log10(m_high)),
    )
    sig1p, edges, _ = stats.binned_statistic(
        ms, gs, lambda xs: np.percentile(xs, ps[3]), bins=hist_bins,
        range=(np.log10(m_low), np.log10(m_high)),
    )
    sig2p, edges, _ = stats.binned_statistic(
        ms, gs, lambda xs: np.percentile(xs, ps[4]), bins=hist_bins,
        range=(np.log10(m_low), np.log10(m_high)),
    )

    mids = (edges[1:] + edges[:-1]) / 2
    plt.plot(mids, sig2m, c=c, lw=2)
    plt.plot(mids, sig1m, c=c, lw=2)
    plt.plot(mids, med, c=c, lw=2)
    plt.plot(mids, sig1p, c=c, lw=2)
    plt.plot(mids, sig2p, c=c, lw=2)


def plot_lims():
    for lim in lims:
        lo, hi = plt.ylim()
        plt.plot(np.log10([lim, lim]), [lo, hi], "k", lw=3)
        plt.ylim(lo, hi)

for (fname, c) in zip(fs, cs)[::-1]:
    ms, gs = read(fname)
    plt.figure(0)
    plot_hist(np.log10(ms), c)
    plt.figure(1)
    plot_percentile_fill(np.log10(ms), gs, c)
for (i, (fname, c)) in enumerate(zip(fs, cs)):
    ms, gs = read(fname)
    plot_percentile_line(np.log10(ms), gs, c)
    plt.plot(np.log10([lims[i], lims[i+1]]), [-0.5, -0.5], "--%s" % c, lw=3,
             label=r"$%g\ {\rm Mpc}/h$" % Ls[i])

plt.figure(0)
plt.ylim(0, 600)
plot_lims()
plt.xlabel(r"$\log_{10}M/M_\odot$")
plt.ylabel(r"$N$")

plt.figure(1)
plot_lims()
plt.xlabel(r"$\log_{10}M/M_\odot$")
plt.ylabel(r"$\Gamma$")
plt.xlim(np.log10(m_low), np.log10(m_high))
plt.legend(loc="upper left", frameon=False)

plt.figure(2)
for i in xrange(len(cs)):
    plt.plot([Ls[i], Ls[i]], [lims[i] / pms[i], lims[i+1] / pms[i]],
             c=cs[i], lw=3, label=r"$%g\ {\rm Mpc}/h$" % Ls[i])
    plt.plot([Ls[i], Ls[i]], [lims[i] / pms[i], lims[i+1] / pms[i]],
             "o%s" % cs[i])

plt.legend(frameon=False, loc="lower right")
plt.xlim(0, 2200)
plt.yscale("log")
plt.xlabel(r"$L\ [{\rm Mpc}/h]$")
plt.ylabel(r"Particles per halo")

# Now let's actually handle splitting up the boxes

plt.figure()
for i in xrange(len(cs)):
    ms, gs = read(fs[i])
    mask = (ms < lims[i + 1]) & (ms >= lims[i])
    gs = gs[mask]
    ms = ms[mask]
    plt.plot(ms, gs, "o", c=cs[i], label=r"$%g\ {\rm Mpc}/h$" % Ls[i])

plt.ylim(-1, 7)
plt.xlabel(r"$M/M_\odot$")
plt.ylabel(r"$\Gamma$")
plt.xscale("log")
plt.legend(frameon=False, loc="upper left")

def mask(mask, *args):
    out = [None] * len(args)
    for i in xrange(len(args)):
        out[i] = args[i][mask]
    return out

def filter_bin(m_low, m_high, g_low, g_high, ms, gs, n_lim=50):
    mask = ((ms < m_high) & (ms >= m_low) & 
            (gs < g_high) & (gs >= g_low))
    idxs = np.arange(len(ms))[mask]
    if len(idxs) <= n_lim: return mask
    idxs = rand.choice(idxs, n_lim, replace=False)
    mask = np.zeros(len(mask), dtype=bool)
    mask[idxs] = True
    return mask

ms, gs, fis = [], [], []
ids, snaps = np.zeros(0, dtype=int), np.zeros(0, dtype=int)
for i in xrange(len(cs)):
    f_ms, f_gs = read(fs[i])
    f_ids, f_snaps = read_ids(fs[i])
    f_fis = i * np.ones(len(f_ms))

    f_ms, f_gs, f_ids, f_snaps, f_fis = mask(
        (f_ms < lims[i+1]) & (f_ms >= lims[i]),
        f_ms, f_gs, f_ids, f_snaps, f_fis,
    )

    ms = np.append(ms, f_ms)
    gs = np.append(gs, f_gs)
    ids = np.append(ids, f_ids)
    snaps = np.append(snaps, f_snaps)
    fis = np.append(fis, f_fis)

ms = np.log10(ms)

m_min, m_max = 11, 16
g_min, g_max = 0, 7
m_bins = 20
g_bins = 7
dm = (m_max - m_min) / m_bins
dg = (g_max - g_min) / g_bins
n_bin = 20

tot_mask = np.zeros(len(ms), dtype=bool)
for mi in xrange(m_bins):
    m_lo, m_hi = mi*dm + m_min, (mi+1)*dm + m_min
    for gi in xrange(g_bins):
        g_lo, g_hi = gi*dg + g_min, (gi+1)*dg + g_min
        tot_mask = tot_mask | filter_bin(
            m_lo, m_hi, g_lo, g_hi, ms, gs, n_bin,
        )

plt.figure()
plt.hist(gs, range=(g_min, g_max), bins=g_bins, color="red")
ms, gs, ids, snaps, fis = mask(
    tot_mask, ms, gs, ids, snaps, fis,
)
plt.hist(gs, range=(g_min, g_max), bins=g_bins, color="blue")
plt.ylabel("$N$")
plt.xlabel(r"$\Gamma$")
plt.ylim(0, 500)

plt.figure()
plt.hist(ms, bins=hist_bins, range=(np.log10(m_low), np.log10(m_high)))
plt.ylabel("$N$")
plt.xlabel(r"$\log_{10} M/M_\odot$")

plt.figure()
plt.plot(10**ms, gs, "o")
plt.ylim(-1, 7)
plt.xscale("log")
plt.ylabel(r"$\Gamma$")
plt.xlabel(r"$M/M_\odot$")

def write(L, ids, snaps, ms, gs):
    mg_name = "out/L%d_unif_mass_gamma.dat" % L
    id_name = "out/L%d_unif_id.dat" % L

    with open(id_name, "w+") as fp:
        for i in xrange(len(ids)):
            print(ids[i], snaps[i], file=fp)

    with open(mg_name, "w+") as fp:
        for i in xrange(len(ids)):
            print(ids[i], snaps[i], ms[i], gs[i], file=fp)

for i in xrange(len(cs)):
    fmask = fis == i
    write(Ls[i], ids[fmask], snaps[fmask], 10**ms[fmask], gs[fmask])

plt.show()
