import parse
import numpy as np
import matplotlib.pyplot as plt
import sys


rhos_mat = np.loadtxt(sys.argv[1]).T
N = int(sys.argv[2])

ngp_rhos, cic_rhos, mc10_rhos, mc100_rhos = rhos_mat[0], rhos_mat[1], rhos_mat[2], rhos_mat[3]

vals = [N/2]
rhos = [ngp_rhos, cic_rhos, mc10_rhos, mc100_rhos]

names = ["NGP", "CIC", "MC10", "MC100"]

cmap_name = "cubehelix"

for val in vals:
    for (rho, name)in zip(rhos, names):
        rho = np.clip(np.reshape(rho, (N, N, N)), 10**-2.5, 10**3)

        plt.figure()
        plt.title("%s" % name)
        plt.imshow(np.mean(np.log10(rho[:, :, val-1:val+2]), axis=2), cmap=cmap_name)
        plt.colorbar()

plt.show()
