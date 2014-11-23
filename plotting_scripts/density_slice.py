import parse
import numpy as np
import matplotlib.pyplot as plt
import sys


rhos = np.loadtxt(sys.argv[1])
N = int(sys.argv[2])
rhos = np.reshape(rhos, (N, N, N)) / np.mean(rhos)
# rhos = np.clip(rhos, np.min(1e-2), np.max(rhos))

rhos = np.clip(np.reshape(rhos, (N, N, N)), 10**-2.5, 10**3)

vals = (np.array(range(9)) + 1) * N/10
cmap_name = "cubehelix" #"gist_heat" #"cubehelix"

for val in vals:
    plt.figure()
    plt.title(r"N = %d" % val)
    plt.imshow(np.log10(rhos[:, :, val]), cmap=cmap_name)
    plt.colorbar()

plt.show()
