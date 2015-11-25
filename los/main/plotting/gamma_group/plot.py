from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.signal as signal
import scipy.stats as stats

def loadtxt(fname):
    with open(fname) as fp: s = fp.read()
    lines = s.split("\n")
    lines = [[float(tok) for tok in line.split(" ") if tok != ""] for line in lines]
    return np.array([np.array(line) for line in lines if len(line) > 0])

mass_rows = np.loadtxt(sys.argv[1])
mean_rows = np.loadtxt(sys.argv[2])
med_rows = np.loadtxt(sys.argv[3])
sub_rows = loadtxt(sys.argv[4])


def to_prof(rows):
    n = (len(rows[0]) - 2) / 2
    return (np.array([row[2:n+2] for row in rows]),
            np.array([row[n+2:] for row in rows]))

m_sp, r_sp, m_200, r_200, r_vir, r_s, gamma, m_200c = map(np.array, zip(*mass_rows))[2:]

mean_rs, mean_rhos = to_prof(mean_rows)
med_rs, med_rhos = to_prof(med_rows)

m_lim = 3e12
PLOT_PROFS = False
PLOT_MEAN = False
PLOT_MEDIAN = False

m_sp = m_sp[m_200 >= m_lim]
r_sp = r_sp[m_200 >= m_lim]
r_200 = r_200[m_200 >= m_lim]
gamma = gamma[m_200 >= m_lim]
sub_rows = sub_rows[m_200 >= m_lim]
#gamma = 1/gamma

med_rhos = med_rhos[m_200 >= m_lim]
med_rs = med_rs[m_200 >= m_lim]
mean_rhos = mean_rhos[m_200 >= m_lim]
mean_rs = mean_rs[m_200 >= m_lim]

m_200c = m_200c[m_200 >= m_lim]
m_200 = m_200[m_200 >= m_lim]

msubs = []
for i in xrange(len(sub_rows)):
    vals = sub_rows[i][2:]
    n = len(vals) // 2
    ms = vals[n:]
    msubs.append(np.max(ms))
msubs = np.array(msubs)

bins = 6
bin_low, bin_high = 0, 6
dbin = (bin_high - bin_low) / bins

r_bins = [[] for _ in range(bins)]
r_200_bins = [[] for _ in range(bins)]
med_rho_bins = [[] for _ in range(bins)]
mean_rho_bins = [[] for _ in range(bins)]

def gamma_to_r_rat(gamma):
    return 0.54 * (1 + 0.53 * 0.27) * (1 + 1.36 * np.exp(-gamma / 3.04))

for i in xrange(len(m_sp)):
    if gamma[i] >= bin_high or gamma[i] <= bin_low: continue
    bin = int((gamma[i] - bin_low) / dbin)

    r, rho = mean_rs[i], mean_rhos[i]
    r_bins[bin].append(r)
    mean_rho_bins[bin].append(rho)
    r_200_bins[bin].append(r_200[i])

    if PLOT_PROFS:
        plt.figure(bin + bins)
        plt.title(r"$\Gamma \in [%.1f,\ %.1f]$" %
                  (bin * dbin + bin_low, (bin + 1) * dbin + bin_low))
        plt.plot(np.log10(r / r_200[i]), rho, "k")

    r, rho = med_rs[i], med_rhos[i]

    med_rho_bins[bin].append(rho)

def log_mean(vals):
    return np.median(vals, axis=0)
    #return np.exp(np.sum(np.log(vals), axis=0) / len(vals))

def combine(i):
    r_nil = r_bins[i][0] / r_200_bins[i][0]
    return r_nil, log_mean(mean_rho_bins[i]), log_mean(med_rho_bins[i])

def find_min(xs, vals):
    min_i = signal.argrelextrema(vals, np.less)
    min_i = min_i[(min_i > 0) & (min_i < (len(vals) - 1))]
    min_i_i = np.argmin(vals[min_i])
    return xs[min_i[min_i_i]]

gamma_sp = []
med_sp = []
mean_sp = []

for i in xrange(bins):
    r_nil, rho_mean, rho_med = combine(i)

    if PLOT_PROFS:
        plt.figure(i)
        plt.title(r"$\Gamma \in [%.1f,\ %.1f]$" %
                  (i * dbin + bin_low, (i + 1) * dbin + bin_low))

    dr = np.log10(r_nil[1]) - np.log10(r_nil[0])
    drdr_mean = signal.savgol_filter(
        np.log10(rho_mean), 61, 4, deriv=1, delta=dr,
    )
    drdr_med = signal.savgol_filter(
        np.log10(rho_med), 61, 4, deriv=1, delta=dr,
    )

    med_min = find_min(r_nil, drdr_med)
    mean_min = find_min(r_nil, drdr_mean)
    
    gamma_sp.append((i + 0.5) * dbin)
    med_sp.append(med_min)
    mean_sp.append(mean_min)

    if PLOT_PROFS:
        plt.ylim(-8, 0)
        plt.plot(np.log10([mean_min, mean_min]), [-8, 0], "r", lw=3)
        plt.plot(np.log10(r_nil), drdr_mean, "r", lw=3, label="Stacked Mean")
        plt.plot(np.log10([med_min, med_min]), [-8, 0], "b", lw=3)
        plt.plot(np.log10(r_nil), drdr_med, "b", lw=3, label="Stacked Median")
        plt.legend(loc="lower left", frameon=False)
        plt.xlabel(r"$\log_{10}\ R_{\rm sp} / R_{\rm 200m}$")
        plt.ylabel(r"$d\log\rho / d\log r$")

        plt.figure(i + bins)
        plt.plot(np.log10([mean_min, mean_min]), [1e-1, 1e3], "r", lw=3)
        plt.plot(np.log10(r_nil), rho_mean, "r", lw=3, label="Stacked Mean")
        plt.plot(np.log10([med_min, med_min]), [1e-1, 1e3], "b", lw=3)
        plt.plot(np.log10(r_nil), rho_med, "b", lw=3, label="Stacked Median")
        plt.legend(loc="lower left", frameon=False)
        plt.yscale("log")
        plt.ylim(1e-1, 1e3)
        plt.xlabel(r"$\log_{10}\ R_{\rm sp} / R_{\rm 200m}$")
        plt.ylabel(r"$\rho / \rho_{\rm m}$")


#plt.figure()
#plt.hist(gamma, bins=bins, range=(bin_low, bin_high))

med_r_sp = []
for i in xrange(len(med_rs)):
    r, rho = med_rs[i], med_rhos[i]
    dr = np.log10(r[1]) - np.log10(r[0])
    drdr = signal.savgol_filter(
        np.log10(rho), 61, 4, deriv=1, delta=dr,
    )
    sp = find_min(r, drdr)
    med_r_sp.append(sp / r_200[i])

plt.figure()

def plot_binned_range(xs, ys, c, bins, label, low=bin_low, high=bin_high):
    ys = np.array(ys)
    """
    mean, edges, _ = stats.binned_statistic(
        xs, ys, "mean", bins=bins, range=(bin_low, bin_high),
    )
    sqr, _, _ = stats.binned_statistic(
        xs, ys*ys, "mean", bins=bins, range=(bin_low, bin_high),
    )
    std = np.sqrt(sqr - mean*mean)
    bin_xs = (edges[1:] + edges[:-1]) / 2

    low, high, mid = mean-std, mean+std, mean
    """

    mid, edges, _ = stats.binned_statistic(
        xs, ys, "median", bins=bins, range=(low, high),
    )
    low_curve, _, _ = stats.binned_statistic(
        xs, ys, lambda y: np.percentile(y, 16),
        bins=bins, range=(low, high),
    )
    high_curve, _, _ = stats.binned_statistic(
        xs, ys, lambda y: np.percentile(y, 84),
        bins=bins, range=(low, high),
    )

    bin_xs = (edges[1:] + edges[:-1]) / 2

    plt.plot(bin_xs, high_curve, c, lw=1)
    plt.plot(bin_xs, low_curve, c, lw=1)
    plt.fill_between(bin_xs, low_curve, high_curve, facecolor=c, alpha=0.3)
    plt.plot(xs, ys, "%s." % c)
    plt.plot(bin_xs, mid, c, lw=3, label=label)

if PLOT_MEDIAN:
    plot_binned_range(gamma, med_r_sp, "g", 8, "LoS Median $R_{\rm sp}$")
cs = ["k", "b", "g", "r"]
lims = [(0, 0.15), (0.15, 0.3), (0.3, 0.6), (0.6, 2.0)]
plot_binned_range(gamma, (r_sp / r_200), "b",
                  8, r"gotetra $R_{\rm sp}$")
#for (c, lim_set) in zip(cs, lims):
#    lo_lim, hi_lim = lim_set
#    lim_mask = ((msubs / m_200c) > lo_lim) & ((msubs / m_200c) <= hi_lim)
#    plot_binned_range(gamma[lim_mask], (r_sp / r_200)[lim_mask], c, 4,
#                      (r"$M_{\rm 200c,sub}/M_{\rm 200c} \in [%g, %g]$" %
#                       (lo_lim, hi_lim)))
"""
for (c, lim_set) in zip(cs, lims):
    lo_lim, hi_lim = lim_set
    lim_mask = ((msubs / m_200c) > lo_lim) & ((msubs / m_200c) <= hi_lim)
    plt.plot(gamma[lim_mask], (r_sp / r_200)[lim_mask], ".%s" % c)
"""
if PLOT_MEAN:
    plt.plot(gamma_sp, mean_sp, "r", lw=2)
    plt.plot(gamma_sp, mean_sp, "ro", label=r"Stacked $R_{\rm sp}$")
#plt.plot(gamma_sp, med_sp, "b", lw=2)
#plt.plot(gamma_sp, med_sp, "bo", label="Stacked Median")

plt.xlim(bin_low, bin_high)
gamma_range = np.linspace(bin_low, bin_high, 100)
fit_range = gamma_to_r_rat(gamma_range)
plt.plot(gamma_range, fit_range, "k", lw=3, label="M+15 Fit")

plt.legend(loc="upper right")
plt.xlabel(r"$\Gamma$")
plt.ylabel(r"$R_{\rm sp} / R_{\rm 200m}$")
plt.ylim(0.5, 2.5)

plt.figure()
plot_binned_range(np.log10(r_200), r_sp / r_200, "m", 8, "Shell",
                  low=-0.6, high=0.5)
plt.ylabel(r"$R_{\rm sp} / R_{\rm 200 m}$")
plt.xlabel(r"$R_{\rm 200 m}$")
plt.legend()

plt.figure()
plot_binned_range(np.log10(m_200), m_sp / m_200, "m", 8, "Shell",
                  low=12, high=15)
plt.ylabel(r"$M_{\rm sp} / M_{\rm 200 m}$")
plt.xlabel(r"$M_{\rm 200 m}$")
plt.legend()

plt.figure()
plt.plot(m_200, r_sp, "m.")
plt.xscale("log")

plt.figure()
plt.plot(gamma, msubs / m_200c, ".m")

plt.show()
