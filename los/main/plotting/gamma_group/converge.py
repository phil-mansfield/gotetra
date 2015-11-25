import matplotlib.pyplot as plt
import numpy as np
import sys

file1024 = sys.argv[1]
file512 = sys.argv[2]
file256 = sys.argv[3]
file128 = sys.argv[4]

m1024, r1024, m200, r200 = map(np.array, zip(*np.loadtxt(file1024, usecols=(2, 3, 4, 5))))
m512, r512 = map(np.array, zip(*np.loadtxt(file512, usecols=(2, 3))))
m256, r256 = map(np.array, zip(*np.loadtxt(file256, usecols=(2, 3))))
m128, r128 = map(np.array, zip(*np.loadtxt(file128, usecols=(2, 3))))

ms, rs, ns = [m128, m256, m512], [r128, r256, r512], [128, 256, 512]
cs = ["r", "b", "g"]

plt.figure(0)
plt.xlim(100, 1000)
plt.plot([100, 1000], [0, 0], "k")
plt.figure(1)
plt.xlim(100, 1000)
plt.ylim(0.005, 0.015)
plt.plot([100, 1000], [0, 0], "k")
plt.figure(2)
plt.xlim(100, 1000)
plt.plot([100, 1000], [0, 0], "k")
plt.figure(3)
plt.xlim(100, 1000)
plt.ylim(0.005, 0.015)
plt.plot([100, 1000], [0, 0], "k")

for (i, (m, r, n, c)) in enumerate(zip(ms, rs, ns, cs)):
    dm = m - m1024
    low2, low, mid, high, high2 = np.percentile(
        dm / m1024, [4.5, 16, 50, 84, 95.5],
    )
    low, high = mid - low, high - mid
    low2, high2 = mid - low2, high2 - mid

    plt.figure(0)
    plt.errorbar([n], [mid], yerr=[[low], [high]],
                 fmt="o", c=c, label="%d" % n)
    plt.errorbar([n], [mid], yerr=[[low2], [high2]],
                 fmt="o", c=c)

    plt.figure(1)
    plt.errorbar([n], np.std(dm / m1024), fmt="o", label="%d" % n)

    dr = r - r1024
    low2, low, mid, high, high2 = np.percentile(
        dr / r1024, [2, 16, 50, 84, 98],
    )
    low, high = mid - low, high - mid
    low2, high2 = mid - low2, high2 - mid

    plt.figure(2)
    plt.errorbar([n], [mid], yerr=[[low], [high]],
                 fmt="o", c=c, label="%d" % n)
    plt.errorbar([n], [mid], yerr=[[low2], [high2]],
                 fmt="o", c=c)

    plt.figure(3)
    plt.errorbar([n], np.std(dr / r1024), fmt="o", label="%d" % n)
 
plt.figure(0)
plt.xscale("log")
plt.legend(loc="lower right", frameon=False)
plt.xlabel("")
plt.ylabel(r"Median Error: $(M_{{\rm sp},n} - M_{{\rm sp}, 1024})/M_{{\rm sp}, 1024}$")

plt.figure(1)
plt.xscale("log")
plt.legend(loc="lower left", frameon=False)
plt.xlabel("")
plt.ylabel(r"RMS Error: $\sqrt{\langle ((M_{{\rm sp},n} - M_{{\rm sp}, 1024}) / M_{{\rm sp}, 1024})^2\rangle }$")

plt.figure(2)
plt.xscale("log")
plt.legend(loc="lower right", frameon=False)
plt.xlabel("")
plt.ylabel(r"Median Error: $(R_{{\rm sp},n} - R_{{\rm sp}, 1024})/R_{{\rm sp}, 1024}$")

plt.figure(3)
plt.xscale("log")
plt.legend(loc="upper right", frameon=False)
plt.xlabel("")
plt.ylabel(r"RMS Error: $\sqrt{\langle ((R_{{\rm sp},n} - R_{{\rm sp}, 1024}) / R_{{\rm sp}, 1024})^2\rangle }$")

plt.show()
