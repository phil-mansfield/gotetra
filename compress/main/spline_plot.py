import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as intr
import sys

#INTR = "linear"
INTR = "cubic"
H_MIN, H_MAX = 3.25, 10
MID_X = 6.82

L = 62.5
x0s, xs, ys, zs = map(np.array, zip(*np.loadtxt(sys.argv[1])))
n = len(x0s)
x0s, xs, ys, zs = x0s[:n-2], xs[:n-2], ys[:n-2], zs[:n-2]

flag = (x0s > H_MIN) & (x0s < H_MAX)
x0s = x0s[flag]
xs = xs[flag]
ys = ys[flag]
zs = zs[flag]

skip_cs = {
    1: "r",
    2: "g",
    4: "b",
    8: "m",
}

def plot_interp(vals, skip):
    c = skip_cs[skip]
    idxs = np.arange(0, len(vals), skip)

    if idxs[-1] != len(vals) - 1:
        idxs = np.append(idxs, [len(vals) - 1])

    f_lin = intr.interp1d(x0s[idxs], vals[idxs], INTR, assume_sorted=True)
    ev_x0s = np.linspace(x0s[0], x0s[-1], len(x0s) * 20)
    ev_xs = f_lin(ev_x0s)
    plt.plot(ev_x0s, ev_xs, c, lw=2, label="Skip = %d" % skip)
    lo, hi = np.min(ev_xs), np.max(ev_xs)
    xlo, xhi = plt.xlim()
    plt.plot([xlo, xhi], [lo, lo], "--", c=c)
    plt.plot([xlo, xhi], [hi, hi], "--", c=c)
    xlo, xhi = plt.xlim()
    return hi - MID_X, MID_X - lo

plt.figure()
plt.plot(x0s, xs, "k")
plt.plot(x0s, xs, "ok")
hi1, lo1 = plot_interp(xs, 1)
hi2, lo2 = plot_interp(xs, 2)
hi4, lo4 = plot_interp(xs, 4)
hi8, lo8 = plot_interp(xs, 8)
xlo, xhi = plt.xlim()
plt.plot([xlo, xhi], [MID_X, MID_X], "k")
plt.xlim(xlo, xhi)

print "%6s %6s %6s" % ("", "f_hi", "f_lo")
print "%6s %6.3f %6.3f" % ("f_12", (hi1 - hi2) / hi2, (lo1 - lo2) / lo2)
print "%6s %6.3f %6.3f" % ("f_24", (hi2 - hi4) / hi4, (lo2 - lo4) / lo4)
print "%6s %6.3f %6.3f" % ("f_48", (hi4 - hi8) / hi8, (lo4 - lo8) / lo8)

plt.legend(loc="lower right")
"""
plt.figure()
plt.plot(x0s, ys, "k")
plt.plot(x0s, ys, "ok")
plot_interp(ys)
plt.legend()
plt.figure()
plt.plot(x0s, zs, "k")
plt.plot(x0s, zs, "ok")
plot_interp(zs)
plt.legend()
"""
plt.show()
