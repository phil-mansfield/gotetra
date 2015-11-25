import numpy as np
import matplotlib.pyplot as plt
import os.path as path

dir = "data"
phases = ["h100_phase_sub1.dat", "h100_phase_sub2.dat",
          "h100_phase_sub4.dat"]
phases = [path.join(dir, filename) for filename in phases]

rads = ["h100_rad_sub1.dat", "h100_rad_sub2.dat",
        "h100_rad_sub4.dat"]
rads = [path.join(dir, filename) for filename in rads]

for target_row in xrange(len(np.loadtxt(rads[0]))):
    for (phase, rad) in zip(phases, rads):
        row = np.loadtxt(phase)[target_row]
        id, snap, v_low, v_high, r_low, r_high, v_width, r_width = row[:8]
        grid = row[8:]
        grid = np.reshape(grid, (v_width, r_width))
        
        row = np.loadtxt(rad)[target_row]
        id, snap, _, r_sp, r_min, r_max, r200m, m200c, gamma = row

        plt.figure()
        plt.title(r"$M_{\rm 200c} = %.1g$ $\Gamma = %.1g$" % (m200c, gamma))
        plt.imshow(grid, cmap="afmhot", extent=[r_low, r_high, v_low, v_high],
                   aspect="auto", origin="lower")
        plt.xlabel(r"$R$ [Mpc/$h$]")
        plt.ylabel(r"$v_{\rm r}$ [km/s]")
        plt.colorbar()
        lo, hi = plt.ylim()
        plt.plot([r_sp, r_sp], [lo, hi], "w", lw=3)
        plt.plot([r_min, r_min], [lo, hi], "w", lw=1)
        plt.plot([r_max, r_max], [lo, hi], "w", lw=1)
        plt.ylim(lo, hi)

plt.show()
    
