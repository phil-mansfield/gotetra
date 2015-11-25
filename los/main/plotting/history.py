import numpy as np
import matplotlib.pyplot as plt
import deriv

#plot_type = "mass"
#plot_type = "radius"
plot_type = "gamma"
#plot_type = "dmdt"

H0 = 1.0 / 14e9

def H(a):
    return np.sqrt(0.3 / a**3 + 0.7) * H0

def dMdt(M, a):
    return deriv.vector_deriv(a, M) * H(a) * a

def is_valid(x):
    return (~np.isinf(x)) & (~np.isnan(x)) & (x > 0)

g_snaps, g_m_sp, g_r_sp, g_a, g_m_200m, g_m_200c, g_r_200m, g_r_200c = zip(*np.loadtxt("plot.dat", usecols=(1, 2, 3, 4, 5, 6, 7, 8)))
z0s = []
for (i, snap) in enumerate(g_snaps):
    if int(snap) == 100: z0s.append(i + 1)

starts = [0] + z0s[:-1]
ends = z0s

for (start, end) in zip(starts, ends):
    a = np.array(g_a[start:end])
    m_sp = np.array(g_m_sp[start:end])
    r_sp = np.array(g_r_sp[start:end])
    m_200m = np.array(g_m_200m[start:end])
    m_200c = np.array(g_m_200c[start:end])
    r_200m = np.array(g_r_200m[start:end])
    r_200c = np.array(g_r_200c[start:end])

    valid = is_valid(m_sp) & is_valid(r_sp) & is_valid(a)
    a = a[valid]
    m_sp = m_sp[valid]
    r_sp = r_sp[valid]
    m_200m = m_200m[valid]
    m_200c = m_200c[valid]
    r_200m = r_200m[valid]
    r_200c = r_200c[valid]
    
    r_sp = r_sp / a

    plt.figure()

    if plot_type == "mass":
        plt.plot(np.log10(a), np.log10(m_sp), lw=3, label=r"$M_{\rm sp}$")
        plt.plot(np.log10(a), np.log10(m_200m), lw=3, label=r"$M_{\rm 200m}$")
        plt.plot(np.log10(a), np.log10(m_200c), lw=3, label=r"$M_{\rm 200c}$")
        
        plt.xlabel(r"$\log_{10}a$", fontsize=18)
        plt.ylabel(r"$\log_{10}M$ [$M_\odot$]", fontsize=18)
        plt.legend(loc="lower left", frameon=False, fontsize=12)
        plt.grid()
    elif plot_type == "radius":
        plt.plot(np.log10(a), np.log10(r_sp), lw=3, label=r"$R_{\rm sp}$")
        plt.plot(np.log10(a), np.log10(r_200m), lw=3, label=r"$R_{\rm 200m}$")
        plt.plot(np.log10(a), np.log10(r_200c), lw=3, label=r"$R_{\rm 200c}$")
        
        plt.xlabel(r"$\log_{10}a$", fontsize=18)
        plt.ylabel(r"$\log_{10}R$ [Mpc/$h$]", fontsize=18)
        plt.legend(loc="lower left", frameon=False, fontsize=12)
        plt.grid()
    elif plot_type == "gamma":
        #la = np.log10(a)
        #lm_sp = np.log10(m_sp)
        #lm_200c = np.log10(m_200c)
        #lm_200m = np.log10(m_200m)
        la, lm_sp, lm_200m, lm_200c = np.log10(a), np.log10(m_sp), np.log10(m_200m), np.log10(m_200c)
        gam_sp = deriv.vector_deriv(la, lm_sp)
        gam_200m = deriv.vector_deriv(la, lm_200m)
        gam_200c = deriv.vector_deriv(la, lm_200c)

        plt.plot(la, gam_sp, lw=3, label=r"$R_{\rm sp}$")
        plt.plot(la, gam_200m, lw=3, label=r"$R_{\rm 200m}$")
        plt.plot(la, gam_200c, lw=3, label=r"$R_{\rm 200c}$")

        plt.xlabel(r"$\log_{10}a$", fontsize=18)
        plt.ylabel(r"$\Gamma$", fontsize=18)
        plt.yscale("symlog")
        plt.legend(loc="upper right", frameon=False, fontsize=12)
        plt.grid()
    elif plot_type == "dmdt":
        plt.plot(np.log10(a), dMdt(m_sp, a), lw=3, label=r"$M_{\rm sp}$")
        plt.plot(np.log10(a), dMdt(m_200m, a), lw=3, label=r"$M_{\rm 200m}$")
        plt.plot(np.log10(a), dMdt(m_200c, a), lw=3, label=r"$M_{\rm 200c}$")
        plt.xlabel(r"$\log_{10}a$", fontsize=18)
        plt.ylabel(r"$d M / dt$", fontsize=18)
        plt.yscale("symlog", linthreshy=10)
        plt.legend(loc="upper right", frameon=False, fontsize=12)
        plt.grid()

plt.show()
