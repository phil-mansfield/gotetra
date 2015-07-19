import numpy as np
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if len(sys.argv) != 2:
    print "Correct usage: $ %s caustic_file" % sys.argv[1]
    exit(1)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

rows = np.loadtxt(sys.argv[1])
rs = np.array(zip(*cols)[-1])
max_r = np.max(rs)
for row in rows:
    x, y, z, r = row
    line_r = max_r if np.isnan(r) else r
    term_x, term_y, term_z = x * line_r, y * line_r, z * line_r
    ax.plot([0, term_x], [0, term_y], [0, term_z], c="k")
    if not np.isnan(r):
        ax.scatter([term_x], [term_y], [term_z], c="r")

plt.show()


