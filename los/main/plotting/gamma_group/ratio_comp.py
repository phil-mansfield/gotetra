import numpy as np
import sys
import matplotlib.pyplot as plt

fnames = sys.argv[1:]
ids, s_r_sps, r_200s, gammas = [], [], [], []
for fname in fnames:
    ids, r_sps, r_200s, gammas = map(
        np.array,zip(*np.loadtxt(fname, usecols=(0, 3, 5, 8)))
    )
    s_r_sps.append(r_sps)

#plt.plot(gammas, s_r_sps[0] / r_200s, ".k")
for i in xrange(len(ids)):
    gs = [gammas[i] for _ in fnames]
    rats = [s_r_sps[j][i] / r_200s[i] for j in xrange(len(fnames))]
    plt.plot(gs, rats, "k", lw=2)
    plt.plot(gs[0], rats[0], "ob")
    plt.plot(gs[-1], rats[-1], "or")

plt.ylim(0.6, 2.2)

plt.show()
