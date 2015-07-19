import sys
import scipy.interpolate as intr
import numpy as np
import matplotlib.pyplot as plt

xs = [0, 1, 1.5, 2, 3, 4, 5]
ys = [2, 1, 1, 0, 2, 3, 1]

sp = intr.UnivariateSpline(xs, ys, s=0)
fine_xs = np.linspace(0, 5, 51)
fine_ys = sp(fine_xs)

go_xs, go_ys = zip(*np.loadtxt(sys.argv[1], usecols=(0, 1)))

plt.plot(fine_xs, fine_ys, "b", lw=3)
plt.plot(go_xs, go_ys, "r", lw=3)
plt.plot(xs, ys, "bo")
plt.show()
