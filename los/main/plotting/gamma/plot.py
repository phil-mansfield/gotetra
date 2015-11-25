from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import sys
import deriv
import scipy.signal as signal

plot_type = "deriv"
w = 31
w2 = w // 2

def split(ids, snaps, mvirs, scales):
    split_is = []
    for (i, id) in enumerate(ids):
        if id == -1: split_is.append(i)
    split_is = [-1] + split_is + [len(ids)]

    seg_ids, seg_snaps, seg_mvirs, seg_scales = [], [], [], []
    for i in xrange(len(split_is) - 1):
        start, end = split_is[i] + 1, split_is[i + 1]
        seg_ids.append(ids[start: end])
        seg_snaps.append(snaps[start: end])
        seg_mvirs.append(mvirs[start: end])
        seg_scales.append(scales[start: end])

    return seg_ids, seg_snaps, seg_mvirs, seg_scales

def med_filter(xs):
    med = np.zeros(len(xs))
    for i in xrange(len(xs)):
        med[i] = np.median(xs[max(0, i-w2): i+w2])
    return med

fname = sys.argv[1]
ids, snaps, mvirs, scales = map(np.array, zip(*np.loadtxt(fname)))
ids, snaps, mvirs, scales = split(ids, snaps, mvirs, scales)
for i, (id, snap, mvir, scale) in enumerate(zip(ids, snaps, mvirs, scales)):
    if plot_type == "mass":
        plt.figure()
        plt.plot(scale, mvir, "o", lw=2)
        plt.xlim(0.05, 1)
        plt.ylim(1e8, 1e12)
        plt.xscale("log")
        #plt.yscale("log")
    if plot_type == "deriv":
        plt.figure()
        gamma = deriv.vector_deriv(np.log10(scale), np.log10(mvir))
        plt.plot([0.05, 1], [0, 0], "k")
        plt.plot([1/1.5, 1/1.5], [-2, 10], "k")
        plt.plot(scale, gamma, "r", lw=3, label="FD")
        da = (np.log10(scale[-1]) - np.log10(scale[0])) / (len(scale) - 1)
        savgol_gamma = (
            signal.savgol_filter(np.log10(mvir), w, 4, deriv=1, delta=da,
                                 mode="nearest")
        )
        plt.plot(scale, savgol_gamma, "b", lw=3, label="Sav-Gol")
        plt.plot(scale, signal.savgol_filter(gamma, w, 4, mode="nearest"),
                 "g", lw=3, label="FD + Sav-Gol")
        plt.plot(scale, med_filter(gamma), "k", lw=3, label="Median")

        dk_gamma = (np.log10(mvir)[w:] - np.log10(mvir)[:-w]) / (da * w)
        plt.plot(scale[w:], dk_gamma, "m", lw=4, label="DK14")

        plt.xscale("log")
        plt.ylim(-2, 10)
        plt.xlabel(r"$a$")
        plt.ylabel(r"$\Gamma$")
        plt.legend(loc="lower left")
        plt.xlim(0.05, 1)
        print i, savgol_gamma[-1], signal.savgol_filter(gamma, w, 4, mode="nearest")[-1], med_filter(gamma)[-1], (np.log10(mvir[-1]) - np.log10(mvir[-w2])) / (-np.log10(scale[-w2]))

plt.show()
