import numpy as np
import matplotlib.pyplot as plt
import deriv

plot_type = "mass"

H0 = 1.0 / 14e9

def H(a):
    return np.sqrt(0.3 / a**3 + 0.7) * H0

def dMdt(M, a):
    return deriv.vector_deriv(a, M) * H(a) * a

def is_valid(x):
    return (~np.isinf(x)) & (~np.isnan(x)) & (x > 0)

z0 = 0.5
a0 = 1 / (1 + z0)

g_snaps, g_m_200m, g_a = zip(*np.loadtxt("fast_mass.dat", usecols=(1, 2, 3)))
z0s = []
for (i, snap) in enumerate(g_snaps):
    if int(snap) == 100: z0s.append(i + 1)

starts = [0] + z0s[:-1]
ends = z0s

for (start, end) in zip(starts, ends):
    a = np.array(g_a[start:end])
    m_200m = np.array(g_m_200m[start:end])

    valid = is_valid(a)
    a = a[valid]
    m_200m = m_200m[valid]
    plt.figure()

    if plot_type == "mass":
        plt.plot(np.log10(a), m_200m, lw=3, label=r"$M_{\rm 200m}$")
        plt.yscale("log")
        plt.xlabel(r"$\log_{10}a$", fontsize=18)
        plt.ylabel(r"$\log_{10}M$ [$M_\odot$]", fontsize=18)
        plt.legend(loc="lower left", frameon=False, fontsize=12)
        ylo, yhi = plt.ylim()
        plt.ylim(ylo, yhi)
        plt.plot(np.log10([a0, a0]), [ylo, yhi], "k", lw=3)
        plt.grid()
    elif plot_type == "gamma":
        la, lm_200m = np.log10(a), np.log10(m_200m)
        gam_200m = deriv.vector_deriv(la, lm_200m)

        plt.plot(la, gam_200m, lw=3, label=r"$R_{\rm 200m}$")

        plt.xlabel(r"$\log_{10}a$", fontsize=18)
        plt.ylabel(r"$\Gamma$", fontsize=18)
        plt.yscale("symlog")
        plt.legend(loc="upper right", frameon=False, fontsize=12)
        plt.grid()

plt.show()
