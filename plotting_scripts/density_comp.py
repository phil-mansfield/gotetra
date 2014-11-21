import parse
import numpy as np
import matplotlib.pyplot as plt
import sys


ngp_rhos, cic_rhos, mc10_rhos, mc100_rhos, mc1000_rhos = parse.get_col(sys.argv[1], 0, 1)
N = int(sys.argv[2])

vals = [N/2]
rhos = [ngp_rhos, cic_rhos, mc10_rhos, mc100_rhos, mc1000_rhos]
for rho in rhos:
    rho = np.reshape(rhos, (N, N, N))
names = ["NGP", "CIC", "MC10", "MC100", "MC1000"]

cmap_name = "cubehelix"

for val in vals:
    for (rho, name)in zip(rhos, names):
        plt.figure()
        plt.title("%s" % name)
        plt.imshow(np.log10(rhos[:, :, val]), cmap=cmap_name)
        plt.colorbar()

        plt.figure()
        plt.title("%s" % name)
        plt.imshow(np.log10(rhos[:, val, :]), cmap=cmap_name)
        plt.colorbar()

        plt.figure()
        plt.title("%s" % name)
        plt.imshow(np.log10(rhos[:, val, :]), cmap=cmap_name)
        plt.colorbar()

plt.show()
