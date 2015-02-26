import sys
import matplotlib.pyplot as plt
import numpy as np


THREADS = 14
UNIT_BUFFER_LEN = 2**12

def min_memory(N, G):
    return 3 * N * UNIT_BUFFER_LEN * 4 + G * 8 * THREADS + 64**3 * 4 * 3

def min_throughput(N, G):
    return np.maximum(min_memory(N, G), 3 * N * UNIT_BUFFER_LEN * 12)

# This is wrong.
def optimal_p(G):
    return 2e4 * G / (341 * 342 * 342)

n_axis = 10**np.linspace(0, 7, 50)
g_axis = 10**np.linspace(0, 9, 50)
ns, gs = np.meshgrid(n_axis, g_axis)

mems = min_memory(ns, gs)
throughs = min_throughput(ns, gs)

plt.pcolor(np.log10(gs), np.log10(ns), np.log10(throughs))
y_min, y_max = plt.ylim()
x_min, x_max = plt.xlim()
plt.plot(np.log10(g_axis), np.log10(optimal_p(g_axis)), "k",
         linewidth=2)
plt.ylim((y_min, y_max))
plt.xlim((x_min, x_max))
plt.colorbar()
plt.ylabel(r"$\log_{10} N$")
plt.xlabel(r"$\log_{10} G$")
plt.title(r"$\log_{10} T$ [bytes]")

plt.figure()
plt.pcolor(np.log10(gs), np.log10(ns), np.log10(mems),
           vmin=np.min(np.log10(throughs)),
           vmax=np.max(np.log10(throughs)))
y_min, y_max = plt.ylim()
x_min, x_max = plt.xlim()
plt.plot(np.log10(g_axis), np.log10(optimal_p(g_axis)), "k",
         linewidth=2)
plt.ylim((y_min, y_max))
plt.xlim((x_min, x_max))
plt.ylabel(r"$\log_{10} N$")
plt.xlabel(r"$\log_{10} G$")
plt.title(r"$\log_{10} M$ [bytes]")
plt.colorbar()
plt.show()
