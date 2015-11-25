import numpy as np
import sys
import matplotlib.pyplot as plt
import random
import scipy.stats as stats
import os.path as path

m_lim = 1e11

fnames = sys.argv[1:]
ids, snaps, r_sps, r_200s, gammas = [], [], [], [], []
for fname in fnames:
    f_ids, f_snaps, f_r_sps, f_m_200s, f_r_200s, f_gammas = map(
        np.array,zip(*np.loadtxt(fname, usecols=(0, 1, 3, 4, 5, 8)))
    )
    mask = f_m_200s > m_lim

    ids.append(f_ids[mask])
    snaps.append(f_snaps[mask])
    r_sps.append(f_r_sps[mask])
    r_200s.append(f_r_200s[mask])
    gammas.append(f_gammas[mask])

cs = ["r", "b", "g", "m", "k"]
for fname, id, r_sp, r_200, gamma, c in zip(fnames, ids, r_sps, r_200s, gammas, cs):
    #plt.figure()
    plt.ylim(0.6, 2.2)
    plt.xlim(0, 6)
    plt.plot(gamma, r_sp / r_200, ".", c=c)
    
    low_bin = np.floor(np.min(gamma))
    high_bin = np.ceil(np.max(gamma))
    n_bins = int(high_bin) - int(low_bin)

    meds, edges, _ = stats.binned_statistic(
        gamma, r_sp / r_200, "median", bins=n_bins, range=(low_bin, high_bin),
    )
    mids = (edges[1:] + edges[:-1]) / 2
    plt.plot(mids, meds, lw=3, c=c, label=r"%s" % (path.basename(fname)))

    plt.xlabel(r"$\Gamma$")
    plt.ylabel(r"$R_{\rm sp} / R_{\rm 200m}$")

plt.legend(loc="upper right")

n_gamma, n_ratio = 3, 3
low_gamma, low_ratio = 0, 0.6
d_gamma, d_ratio = 2, 0.4
def print_box_ids(fname, ids, snaps, ratios, gammas):
    bins = [[[] for _ in xrange(n_gamma)] for _ in xrange(n_ratio)]
    for i in xrange(len(ratios)):
        id, snap, ratio, gamma = ids[i], snaps[i], ratios[i], gammas[i]
        i_gamma = np.array((gamma - low_gamma) / d_gamma, dtype=int)
        i_ratio = np.array((ratio - low_ratio) / d_ratio, dtype=int)
        if i_ratio >= n_ratio: continue
        if i_gamma >= n_gamma: continue
        bins[i_ratio][i_gamma].append((id, snap, ratio, gamma))
    print "# %s" % fname
    for i_ratio in xrange(n_ratio):
        for i_gamma in xrange(n_gamma):
            if len(bins[i_ratio][i_gamma]) == 0: continue
            print "%d %d" % random.choice(bins[i_ratio][i_gamma])[:2]

for fname, id, snap, r_sp, r_200, gamma in zip(fnames, ids, snaps, r_sps, r_200s, gammas):
    ratio = r_sp / r_200
    print_box_ids(fname, id, snap, ratio, gamma)

plt.show()
