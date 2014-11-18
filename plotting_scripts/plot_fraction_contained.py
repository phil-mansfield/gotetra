import sys
import numpy as np
import matplotlib.pyplot as plt
import parse

fracs, radii = parse.get_col(sys.argv[1], 0, 1)
inv_fracs = 1.0 - fracs

plot_fs, plot_rs = zip(*[(f, r) for (f, r) in zip(inv_fracs, radii) if f > 0])

radius_lim = (125.0) / 8 * 1.5
inner_lim = radius_lim / 3.0
inner_corner = inner_lim * np.sqrt(2)

plt.plot(plot_rs, np.log10(plot_fs), "k", linewidth=2)
plt.plot([0, plot_rs[0]], [np.log10(plot_fs[0]), np.log10(plot_fs[0])],
         "k--", linewidth=2)

plt.plot([inner_lim, inner_lim], [plt.ylim()[0], plt.ylim()[1]], "b--",
         label="radius of inner subbox edge")
plt.plot([inner_corner, inner_corner], [plt.ylim()[0], plt.ylim()[1]], "g--",
         label="radius of inner subbox corner")
plt.plot([radius_lim, radius_lim], [plt.ylim()[0], plt.ylim()[1]], "r--",
         label="radius of outer subbox edge")

plt.xlabel(r"$r$ [Mpc]")
plt.ylabel(r"$\log_{10}(f(r) - 1)$")
plt.xlim((0, radii[-1]))
plt.legend(loc=2)

w = 125 / 8.0 * 3.0

plt.title(
    "Fraction contained within %.1f Mpc subbox of 125 Mpc box at z = 0." % w
)

plt.show()
