from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import sys
import numpy.random as rand
import scipy.stats as stats

fname = sys.argv[1]

ids, snaps, ms, gs = map(np.array, zip(*np.loadtxt(fname)))

n_lim = 25
#n_lim = 1000000
m_min = 12
m_max = 16
g_min = -1
g_max = 6

g_bins = 7
m_bins = 8

dm = (m_max - m_min) / m_bins
dg = (g_max - g_min) / g_bins

arange = np.arange(0, len(ms))
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
for m_bin in xrange(m_bins):
    m_low, m_high = 10**(m_min + m_bin*dm), 10**(m_min + (m_bin + 1)*dm)
    for g_bin in xrange(g_bins):
        g_low, g_high = g_min + g_bin*dg, g_min + (g_bin + 1)*dg
        mask = mask | filter_bin(m_low, m_high, g_low, g_high)

gs = gs[mask]
ms = ms[mask]

plt.scatter(ms, gs)
plt.ylabel(r"$\Gamma$")
plt.ylim(-1, 6)
plt.xlabel(r"$M$ $[M_{\odot}/h]$")
plt.xscale("log")
plt.xlim(0, 2e15)

def array_print(fmt, xs):
    for x in xs: print(fmt % x, end="")
    print()

n, edges = np.histogram(np.log10(ms), range=(12, 16), bins=8)
print("Mass Histogram:")
array_print("%6d", n)
array_print("%6g", edges)
n, edges = np.histogram(gs, range=(-1, 6), bins=7)
print("Gamma Histogram:")
array_print("%6d", n)
array_print("%6g", edges)

print(len(ids[mask]))

th_ids = ids[mask]
th_snaps = snaps[mask]
th_ms = ms
th_gs = gs

def plot_binned_range(xs, ys, c, bins, label, low, high):
    ys = np.array(ys)
    lxs = np.log10(xs)
    mid, edges, _ = stats.binned_statistic(
        lxs, ys, "median", bins=bins, range=(low, high),
    )
    low_curve, _, _ = stats.binned_statistic(
        lxs, ys, lambda y: np.percentile(y, 16),
        bins=bins, range=(low, high),
    )
    high_curve, _, _ = stats.binned_statistic(
        lxs, ys, lambda y: np.percentile(y, 84),
        bins=bins, range=(low, high),
    )

    bin_xs = 10**((edges[1:] + edges[:-1]) / 2)

    plt.plot(bin_xs, high_curve, c, lw=1)
    plt.plot(bin_xs, low_curve, c, lw=1)
    plt.fill_between(bin_xs, low_curve, high_curve, facecolor=c, alpha=0.3)
    #plt.plot(xs, ys, "%s." % c)
    plt.plot(bin_xs, mid, c, lw=3, label=label)

plot_binned_range(ms, gs, "m", 10, "Range", 13.5, 16)

"""
with open("thin_id.dat", "w+") as fp:
    for (id, snap) in zip(th_ids, th_snaps):
        print("%10d %3d" % (id, snap), file=fp)

with open("thin_gamma_mass.dat", "w+") as fp:
    for (id, snaps, ms, gs) in zip(th_ids, th_snaps, th_ms, th_gs):
        print("%10d %3d %10g %12g" % (id, snap, ms, gs), file=fp)
"""
plt.xlim(10**m_min, 10**m_max)
plt.show()
