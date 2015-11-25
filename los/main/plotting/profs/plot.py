import numpy as np
import matplotlib.pyplot as plt
import sys

rows = np.loadtxt(sys.argv[1])
for row in rows:
    rs = row[2:402]
    rhos = row[402:]
    plt.plot(rs, rhos, lw=3)

plt.yscale("log")
plt.xscale("log")
plt.show()
