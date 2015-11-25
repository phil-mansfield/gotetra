import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

rows = np.loadtxt("prof.dat")
if len(rows.shape) == 1:
    rows = [rows]

cs = ["r", "g", "b", "k", "m", "brown", "orange", "gray", "cyan"]
plt.figure(figsize=(10,8))
max_mult = 1.25
max_mult = 5.0
lo, hi = 1e-1, 1e3

for (i, row) in enumerate(rows):
    id, snap = int(row[0]), int(row[1])
    r200m, r200c, r500c, r2500c, m200m = row[-5:]
    data = row[2:-5]
    rs = data[:len(data)/2]
    rhos = data[len(data)/2:]

    rhos = rhos[rs < r200m * max_mult]
    rs = rs[rs < r200m * max_mult]

    #if i // 3 == 2 or (i // 3 == 1 and (i % 3 == 0 or i % 3 == 2)):
    #    plt.plot([r500c, r500c], [4, hi], cs[i // 3])
    #else:
    if i % 3 == 0 and i // 3 == 2:
        #plt.plot([r500c, r500c], [lo, hi], "r")
        plt.fill_between([r500c, r200m], [lo, lo], [hi, hi],
                         facecolor="r", alpha=0.25)
        plt.plot(rs, rhos, "r", lw=3,
                 label=r"$M_{\rm 200m} = $%.1g $M_\odot/h$" % m200m)
    else:
        continue
        plt.plot(rs, rhos, cs[i // 3], lw=3)

plt.xscale("log")
plt.yscale("log")

plt.ylim(lo, hi)
plt.xlim(0.03, 0.9)

ax = plt.gca()
#ax.xaxis.set_major_formatter(ScalarFormatter())
#ax.xaxis.set_minor_formatter(ScalarFormatter())
#ax.set_xticks([0.01, 0.02, 0.04, 0.08, 0.2, 0.4, 0.8])
#ax.set_xticklabels([0.01, 0.02, 0.04, 0.08, 0.2, 0.4, 0.8])

#plt.xlim(0.03, 0.25)

plt.rc('text', usetex=True)
plt.rc('font',size=19)

plt.xlabel(r"$r$ [Mpc/$h$]")
plt.ylabel(r"$\rho$ [$\rho_{\rm m}$]")
#plt.legend(frameon=False, loc="lower left", fontsize=16)
plt.show()
