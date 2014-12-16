import parse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path as path

#rhos_mat = np.loadtxt(sys.argv[1]).T
N = int(sys.argv[2])
mc100_rhos = parse.get_col(sys.argv[1], 0)
#ngp_rhos, cic_rhos, center_rhos, mc10_rhos, mc100_rhos, mc1000_rhos = rhos_mat[0], rhos_mat[1], rhos_mat[2], rhos_mat[3], rhos_mat[4], rhos_mat[5]

#vals = [N/4, (N * 3)/4]
#vals = [87]
vals = range(10, N - 10)
#rhos = [ngp_rhos, cic_rhos, center_rhos, mc10_rhos, mc100_rhos, mc1000_rhos]
#names = ["NGP", "CIC", "Center", "MC10", "MC100", "MC1000"]

rhos = [mc100_rhos]
names = ["MC100"]

cmap_name = "cubehelix"
#cmap_name = "gist_heat"
#cmap_name = "CMRmap"

i = 1
for val in vals:
    for (rho, name)in zip(rhos, names):
        rho = np.clip(np.reshape(rho, (N, N, N), order="C"), 0.1, 1000)
        #rho = np.reshape(rho, (N, N, N), order="C")

        #plt.figure()
        # plt.title("%s" % name)
        print val
        plt.title("%d" % val)
        plt.imshow(np.log10(rho[val, :, :]), cmap=cmap_name)
        #plt.colorbar()

        plt.savefig("%s_%03d.png" % (path.join(sys.argv[3], sys.argv[1].split("/")[-1].split(".")[0]), i))
        i += 1
