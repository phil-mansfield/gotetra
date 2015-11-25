import numpy as np
import matplotlib.pyplot as plt

def loadtxt(fname):
    with open(fname) as fp: s = fp.read()
    lines = s.split("\n")
    lines = [[float(tok) for tok in line.split(" ") if tok != ""] for line in lines]
    return [np.array(line) for line in lines if len(line) > 0]

subs = loadtxt("sub.dat")
full_subs = loadtxt("sub_full.dat")
info = np.loadtxt("mass.dat")

for i in xrange(len(subs)):
    m200, r200, gamma = info[i][2], info[i][3], info[i][4]
    
    plt.figure()
    vals = subs[i][2:]
    full_vals = full_subs[i][2:]
    n = len(vals) / 2
    rs = vals[:n] / r200
    ms = vals[n:]
    full_rs = full_vals[:n] / r200
    rs = rs[ms > 5e8]
    ms = ms[ms > 5e8]
    full_rs = full_rs[ms > 5e8]

    vals, edges = np.histogram(
        rs, weights=1/(2 * np.pi * rs*rs), bins=int(2*len(rs)**0.3333)
    )
    full_vals, full_edges = np.histogram(
        full_rs, weights=1/(2 * np.pi * full_rs*full_rs),
        bins=int(2*len(full_rs)**0.3333)
    )
    mids = (edges[1:] + edges[:-1]) / 2
    full_mids = (edges[1:] + edges[:-1]) / 2
    plt.plot(full_mids, full_vals, "or", lw=3)
    plt.plot(mids, vals, "ob", lw=3)

    plt.title(r"$M_{\rm 200m} = %g$ $\Gamma = %g$" % (m200, gamma))
    plt.xscale("log")
    plt.yscale("log")

plt.show()
