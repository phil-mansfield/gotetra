import parse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path as path

rhos_mat = np.loadtxt(sys.argv[1]).T
N = int(sys.argv[2])

ngp_rhos, cic_rhos, center_rhos, mc100_rhos, seq100_rhos, mc2500_rhos = rhos_mat[0], rhos_mat[1], rhos_mat[2], rhos_mat[3], rhos_mat[4], rhos_mat[5]

vals = [(N * 3)/4]
rhos = [ngp_rhos, cic_rhos, center_rhos, mc100_rhos, seq100_rhos, mc2500_rhos]
names = ["NGP", "CIC", "Center", "MC100", "Sobol 100", "MC 5000"]

ref = rho = np.clip(np.reshape(mc2500_rhos / np.mean(mc2500_rhos), (N, N, N), order="C"), 0.1, 1000)

cmap_name = "cubehelix"
#cmap_name = "gist_heat"
#cmap_name = "CMRmap"

#dmap_name = "bwr"
dmap_name = "RdBu"

i = 1
for val in vals:
    for (rho, name)in zip(rhos, names):
        rho = np.clip(np.reshape(rho / np.mean(rho), (N, N, N), order="C"), 0.1, 1000)
        #rho = np.reshape(rho, (N, N, N), order="C")

        #plt.figure()
        #plt.title("%s" % name)
        #plt.imshow(np.log10(rho[val, :, :]), cmap=cmap_name)
        #plt.colorbar()

        plt.figure()
        plt.title("%s" % name)
        diff = rho[val, :, :] / ref[val, :, :] - 1.0
        lim = np.min(diff)
        plt.imshow(np.clip(diff, lim, -lim), cmap=dmap_name)
        print name, np.mean((rho[val, :, :] / ref[val, :, :] - 1.0)**2), np.mean((rho - ref)**2)
        plt.colorbar()

plt.show()
