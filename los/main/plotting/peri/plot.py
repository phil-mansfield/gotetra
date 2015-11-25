from __future__ import division

import colossus.Cosmology as cosmology
import matplotlib.colors as colors

import matplotlib.pyplot as plt
import numpy as np
import sys
import deriv

#dir = "data"
dir = "multi"

L = 62.5
#L = 125
#L = 250
#L = 500

"""
r_max_scale = 1.0
r_sp_scale = 0.8
r90_scale = 2.1
r200m_scale = 0.7
r200c_scale = 0.4
"""
r_max_scale = 1.0
r_sp_scale = 1.0
r90_scale = 1.0
r200m_scale = 1.0
r200c_scale = 1.0


SUBHALO_LIM = 40
PLOT_INDIV = False

if L == 62.5:
    m_ids = [522, 2991, 4250, 4302, 4565, 6862, 8092]
    sub_file, tree_file = "%s/h63_subs.dat" % dir, "%s/h63_tree.dat" % dir
    rad_file = "%s/h63_rad.dat" % dir
elif L == 125:
    m_ids = [276, 300, 367, 1268, 1568, 1673, 2019, 2189]
    sub_file, tree_file = "%s/h125_subs.dat" % dir, "%s/h125_tree.dat" % dir
    rad_file = "%s/h125_rad.dat" % dir
elif L == 250:
    m_ids = [201, 236, 296, 314, 339, 521, 726, 924, 1477]
    sub_file, tree_file = "%s/h250_subs.dat" % dir, "%s/h250_tree.dat" % dir
    rad_file = "%s/h250_rad.dat" % dir
elif L == 500:
    m_ids = [8, 51, 105, 241, 465, 562, 809, 902]
    sub_file, tree_file = "%s/h500_subs.dat" % dir, "%s/h500_tree.dat" % dir
    rad_file = "%s/h500_rad.dat" % dir
else:
    assert(0)
m_ids = None

G = 4.302e-9 * 0.7

params = {"flat":True, "H0":70, "Om0":0.27,
          "Ob0":0.0469, "sigma8":0.82, "ns":0.95}
cosmo = cosmology.setCosmology("meowCosmo", params)

sub_cols = map(np.array, zip(*np.loadtxt(sub_file)))
s_ids = np.array(sub_cols[0], dtype=int)
h_ids = np.array(sub_cols[2], dtype=int)
tree_cols = map(np.array, zip(*np.loadtxt(tree_file)))
tree_ids = np.array(tree_cols[0], dtype=int)
tree_snaps = np.array(tree_cols[1], dtype=int)
scales, xs, ys, zs, rs, ms = tree_cols[2:]
rad_cols = map(np.array, zip(*np.loadtxt(rad_file)))
rad_ids, _, m_sp, r_sp, r_min, r_max, r200m, m200c, gamma = rad_cols

print "text loaded"

rad_info = {}
for i in xrange(len(rad_ids)):
    rad_info[rad_ids[i]] = (m_sp[i], r_sp[i], r_min[i], r_max[i],
                            r200m[i], m200c[i], gamma[i])


prof_info = {}
def split_profs():
    rows = map(np.array, np.loadtxt(prof_file))
    n = (len(rows[0]) - 2) // 2
    vals = [(row[2:n+2], row[n+2:]) for row in rows]
    ids = np.array([row[0] for row in rows], dtype=int)
    snaps = np.array([row[1] for row in rows], dtype=int)

    start = 0
    for end in xrange(len(snaps)):
        if snaps[end] == 100:
            prof_info[ids[end]] = vals[start:end]
#split_profs()

class Halo(object):
    def __init__(self, start, end):
        self.id = tree_ids[end - 1]

        self.scales = scales[start:end]
        self.ts = cosmo.age(1/self.scales - 1)
        self.xs = xs[start:end]
        self.ys = ys[start:end]
        self.zs = zs[start:end]
        self.rs = rs[start:end]
        self.ms = ms[start:end]
        self.subs = []
    def add_sub(self, h):
        self.subs.append(h)
    def displace(self, host):
        for i in xrange(min(len(self.xs), len(host.xs))):
            self.xs[-i] -= host.xs[-i]
            if self.xs[-i] > L/2: self.xs[-i] -= L
            if self.xs[-i] < -L/2: self.xs[-i] += L
            self.ys[-i] -= host.ys[-i]
            if self.ys[-i] > L/2: self.ys[-i] -= L
            if self.ys[-i] < -L/2: self.ys[-i] += L
            self.zs[-i] -= host.zs[-i]
            if self.zs[-i] > L/2: self.zs[-i] -= L
            if self.zs[-i] < -L/2: self.zs[-i] += L
            
        self.ds = np.sqrt(self.xs**2 + self.ys**2 + self.zs**2)
        self.phis = np.arccos(self.zs / self.ds)
        self.ths = np.arctan2(self.ys, self.xs)
        self.host_rs = host.rs

        self.vxs = deriv.vector_deriv(self.ts, self.xs) * self.scales * 31.54
        self.vys = deriv.vector_deriv(self.ts, self.ys) * self.scales * 31.54
        self.vzs = deriv.vector_deriv(self.ts, self.zs) * self.scales * 31.54
        self.vs = np.sqrt(self.vxs**2 + self.vys**2 + self.zs**2)
        self.vrs = np.abs(deriv.vector_deriv(self.ts, self.ds) * self.scales * 31.54)

    def perihelion(self):
        l = min(len(self.ds), len(self.host_rs))
        i = np.argmin(self.ds[-l:])

        hi = i - l
        return (self.ds[hi], self.ds[hi]/self.host_rs[hi],
                self.phis[hi], self.ths[hi], self.ts[hi])
start = 0
hs = {}
for i in xrange(len(tree_ids)):
    if tree_ids[i] == -1:
        h = Halo(start, i)
        start = i + 1
        hs[h.id] = h
h = Halo(start, len(tree_ids))
hs[h.id] = h

print "created Halos."

hosts = []
for i, id in enumerate(s_ids):
    if s_ids[i] == h_ids[i]:
        hosts.append(hs[id])
        m_sp, r_sp, r_min, r_max, r200m, m200c, gamma = rad_info[id]
        #r_prof, rho_prof = zip(*prof_info[id])
        hs[id].m_sp = m_sp
        hs[id].r_sp = r_sp
        hs[id].r_min = r_min
        hs[id].r_max = r_max
        hs[id].r200m = r200m
        hs[id].gamma = gamma
        hs[id].m200c = m200c
        hs[id].r200c = (m200c/((cosmo.rho_c(0)*200)*1e9*4*np.pi/3))**0.333
        #hs[id].r_prof = r_prof
        #hs[id].rho_prof = rho_prof
        if m_ids is None:
            hs[id].m_id = None
        else:
            hs[id].m_id = m_ids[len(hosts) - 1]
    else:
        hs[h_ids[i]].add_sub(hs[id])
p_r200m, p_r_sp, p_r_max, p_r200c = [], [], [], []
b_r200m, b_r_sp, b_r_max, b_r200c, b_r90 = [], [], [], [], []
r85_r_sp, r90_r_sp, r95_r_sp = [], [], []
edges = None
sub_counts = []

ds_r200m, ds_r_sp, ds_r_max, xps_all = [], [], [], []

for (i, host) in enumerate(hosts):
    if i % 10 == 0: print i
    for sub in host.subs: sub.displace(host)

    ds, dps, xps, phis, ths, ts = [], [], [], [], [], []
    for sub in host.subs:
        dp, xp, _, _, t = sub.perihelion()
        ds.append(sub.ds[-1])
        dps.append(dp)
        xps.append(xp)
        phis.append(sub.phis[-1])
        ths.append(sub.ths[-1])
        ts.append(t)
    ds, dps, phis = np.array(ds), np.array(dps), np.array(phis)
    xps, ths, ts = np.array(xps), np.array(ths), np.array(ts)
    xs = ds / host.rs[-1]
    if PLOT_INDIV:
        plt.figure()
        plt.scatter(xs, xps, c=ts, s=70)
        plt.hot()
        plt.colorbar()

        xlo, xhi = plt.xlim()
        ylo, yhi = plt.ylim()
        yhi, xhi = 3, 3
        ylo, xlo = 0, 0
        hilim = max(xhi, yhi)
        lolim = min(xlo, ylo)
        plt.plot([lolim, hilim], [lolim, hilim], "k")
        plt.xlim(xlo, xhi)
        plt.ylim(ylo, yhi) 
    
        plt.plot([host.r_sp / host.r200m, host.r_sp / host.r200m], [0, yhi],
                 "r", lw=3, label=r"gotetra $R_{\rm sp}$")
        plt.plot([host.r_min / host.r200m, host.r_min / host.r200m], [0, yhi],
                 "r", lw=1)
        plt.plot([host.r_max / host.r200m, host.r_max / host.r200m], [0, yhi],
                 "r", lw=1, label=r"gotetra shell bounds")
        plt.plot([1, 1], [0, yhi], "--k",
                 lw=3, label=r"$R_{\rm 200m}$")

        if host.m_id is not None:
            plt.title("Halo %d" % host.m_id)
        else:
            plt.title(r"$\rm \log_{10}M_{\rm 200c}$ = %.1g $\Gamma$ = %.2f" % 
                      (host.m200c, host.gamma))
        plt.xlabel(r"$R(z=0)/R_{\rm 200m}(z=0)$")
        plt.ylabel(r"$R(z=z_{\rm peri})/R_{\rm 200m}(z=z_{\rm peri})$")
        plt.legend(loc="upper left")

    plt.figure(30)
    plt.title(r"$R_{\rm sp}$")

    def eps_eq(x, y):
        eps = 0.01
        return np.abs(x - y) < eps
    mask = (~eps_eq(xps, xs)) & (xps < 1) & (ds > 0)

    n = np.sum(mask)
    sub_counts.append(n)
    if n > SUBHALO_LIM:
        vals, edges, _ = plt.hist((ds / host.r_sp * r_sp_scale)[mask],
                                  histtype="step",
                                  cumulative=True, normed=True,
                                  bins=50, range=(0, 2), lw=1)
        plt.plot([1, 1], [0, 1.2], "--k")
        plt.ylim(0, 1.2)
        plt.ylabel(r"$N(<R/R_{\rm sp})/N_{\rm tot}$")
        plt.xlabel(r"$R/R_{\rm sp}$")

        b_r_sp.append(vals)

        plt.figure(31)
        plt.title(r"$R_{\rm max}$")
        vals, _, _ = plt.hist((ds / host.r_max * r_max_scale)[mask],
                              histtype="step",
                              cumulative=True, normed=True,
                              bins=50, range=(0, 2), lw=1)
        plt.plot([1, 1], [0, 1.2], "--k")
        plt.ylim(0, 1.2)
        plt.ylabel(r"$N(<R/R_{\rm max})/N_{\rm tot}$")
        plt.xlabel(r"$R/R_{\rm max}$")

        b_r_max.append(vals)

        plt.figure(32)
        plt.title(r"$R_{\rm 200m}$")
        vals, _, _ = plt.hist((ds / host.r200m * r200m_scale)[mask],
                              histtype="step",
                              cumulative=True, normed=True,
                              bins=50, range=(0, 2), lw=1)
        plt.plot([1, 1], [0, 1.2], "--k")
        plt.ylim(0, 1.2)
        plt.ylabel(r"$N(<R/R_{\rm 200m})/N_{\rm tot}$")
        plt.xlabel(r"$R/R_{\rm 200m}$")
        b_r200m.append(vals)

        plt.figure(33)
        plt.title(r"$R_{\rm 200c}$")
        vals, _, _ = plt.hist((ds / host.r200c * r200c_scale)[mask],
                              histtype="step",
                              cumulative=True, normed=True,
                              bins=50, range=(0, 2), lw=1)
        plt.plot([1, 1], [0, 1.2], "--k")
        plt.ylim(0, 1.2)
        plt.ylabel(r"$N(<R/R_{\rm 200c})/N_{\rm tot}$")
        plt.xlabel(r"$R/R_{\rm 200c}$")
        b_r200c.append(vals)

        plt.figure(34)
        r90 = np.percentile(ds[mask], 90)
        plt.title("$R_{90}$")
        vals, _, _ = plt.hist((ds / r90 * r90_scale)[mask], histtype="step",
                              cumulative=True, normed=True,
                              bins=50, range=(0, 2), lw=1)
        plt.plot([1, 1], [0, 1.2], "--k")
        plt.ylim(0, 1.2)
        plt.ylabel(r"$N(<R/R_{90})/N_{\rm tot}$")
        plt.xlabel(r"$R/R_{90}$")
        b_r90.append(vals)

        ds_r200m += list(ds[mask] / host.r200m)
        ds_r_sp += list(ds[mask] / host.r_sp)
        ds_r_max += list(ds[mask] / host.r_max)
        xps_all += list(xps[mask])

        p_r_sp.append(np.sum(mask & (ds < (host.r_sp*r_sp_scale))) / n)
        p_r_max.append(np.sum(mask & (ds < (host.r_max*r_max_scale))) / n)
        p_r200m.append(np.sum(mask & (ds < (host.r200m*r200m_scale))) / n)
        p_r200c.append(np.sum(mask & (ds < (host.r200c*r200c_scale))) / n)

        r85_r_sp.append((np.percentile(ds[mask], 85) - host.r_sp) / host.r_sp)
        r90_r_sp.append((np.percentile(ds[mask], 90) - host.r_sp) / host.r_sp)
        r95_r_sp.append((np.percentile(ds[mask], 95) - host.r_sp) / host.r_sp)

plt.figure(35)
plt.ylim(0, 1)
plt.ylabel("Cumulative Probability")
plt.xlabel(r"$N(<R^*)/N_{\rm tot}$")
plt.hist(p_r200c, normed=True, histtype="step", bins=100, range=(0,1),
         color="m", lw=3, label=r"$R^*=R_{\rm 200c}$", cumulative=True)
plt.hist(p_r200m, normed=True, histtype="step", bins=100, range=(0,1),
         color="r", lw=3, label=r"$R^*=R_{\rm 200m}$", cumulative=True)
plt.hist(p_r_sp, normed=True, histtype="step", bins=100, range=(0,1),
         color="g", lw=3, label=r"$R^*=R_{\rm sp}$", cumulative=True)
plt.hist(p_r_max, normed=True, histtype="step", bins=100, range=(0,1),
         color="b", lw=3, label=r"$R^*=R_{\rm max}$", cumulative=True)
plt.legend(loc="upper left")


xs = (edges[:-1] + edges[1:]) / 2

plt.figure(36)
m2s, m1s, med, p1s, p2s = np.percentile(
    b_r_sp, [50-95/2, 50-68/2, 50, 50+68/2, 50+95/2], axis=0,
)
plt.fill_between(xs, m2s, p2s, facecolor="green", alpha=0.3)
plt.fill_between(xs, m1s, p1s, facecolor="green", alpha=0.3)
plt.plot(xs, p2s, c="g", lw=1)
plt.plot(xs, m2s, c="g", lw=1)
plt.plot(xs, p1s, c="g", lw=1)
plt.plot(xs, m1s, c="g", lw=1)
plt.plot(xs, med, c="g", lw=3)
plt.ylim(0, 1)
plt.plot([1, 1], [0, 1], "--k")
plt.title(r"$R/R_{\rm sp}$")
plt.ylabel(r"$N(<R/R_{\rm sp})/N_{\rm tot}$")
if r_sp_scale == 1.0:
    plt.xlabel(r"$R/R_{\rm sp}$")
else:
    plt.xlabel(r"$R/(%.1f\ R_{\rm sp})$" % r_sp_scale)

plt.figure(37)
m2s, m1s, med, p1s, p2s = np.percentile(
    b_r_max, [50-95/2, 50-68/2, 50, 50+68/2, 50+95/2], axis=0,
)
plt.fill_between(xs, m2s, p2s, facecolor="blue", alpha=0.3)
plt.fill_between(xs, m1s, p1s, facecolor="blue", alpha=0.3)
plt.plot(xs, p2s, c="b", lw=1)
plt.plot(xs, m2s, c="b", lw=1)
plt.plot(xs, p1s, c="b", lw=1)
plt.plot(xs, m1s, c="b", lw=1)
plt.plot(xs, med, c="b", lw=3)
plt.ylim(0, 1)
plt.plot([1, 1], [0, 1], "--k")
plt.title(r"$R/R_{\rm max}$")
plt.ylabel(r"$N(<R/R_{\rm max})/N_{\rm tot}$")
if r_max_scale == 1.0:
    plt.xlabel(r"$R/R_{\rm max}$")
else:
    plt.xlabel(r"$R/(%.1f\ R_{\rm max})$" % r_max_scale)

plt.figure(38)
m2s, m1s, med, p1s, p2s = np.percentile(
    b_r200m, [50-95/2, 50-68/2, 50, 50+68/2, 50+95/2], axis=0,
)
plt.fill_between(xs, m2s, p2s, facecolor="red", alpha=0.3)
plt.fill_between(xs, m1s, p1s, facecolor="red", alpha=0.3)
plt.plot(xs, p2s, c="r", lw=1)
plt.plot(xs, m2s, c="r", lw=1)
plt.plot(xs, p1s, c="r", lw=1)
plt.plot(xs, m1s, c="r", lw=1)
plt.plot(xs, med, c="r", lw=3)
plt.ylim(0, 1)
plt.plot([1, 1], [0, 1], "--k")
plt.title(r"$R/R_{\rm 200m}$")
plt.ylabel(r"$N(<R/R_{\rm 200m})/N_{\rm tot}$")
if r200m_scale == 1.0:
    plt.xlabel(r"$R/R_{\rm 200m}$")
else:
    plt.xlabel(r"$R/(%.1f\ R_{\rm 200m})$" % r200m_scale)

plt.figure(39)
m2s, m1s, med, p1s, p2s = np.percentile(
    b_r200c, [50-95/2, 50-68/2, 50, 50+68/2, 50+95/2], axis=0,
)
plt.fill_between(xs, m2s, p2s, facecolor="magenta", alpha=0.3)
plt.fill_between(xs, m1s, p1s, facecolor="magenta", alpha=0.3)
plt.plot(xs, p2s, c="m", lw=1)
plt.plot(xs, m2s, c="m", lw=1)
plt.plot(xs, p1s, c="m", lw=1)
plt.plot(xs, m1s, c="m", lw=1)
plt.plot(xs, med, c="m", lw=3)
plt.ylim(0, 1)
plt.plot([1, 1], [0, 1], "--k")
plt.title(r"$R/R_{\rm 200c}$")
plt.ylabel(r"$N(<R/R_{\rm 200c})/N_{\rm tot}$")
if r200c_scale == 1.0:
    plt.xlabel(r"$R/R_{\rm 200c}$")
else:
    plt.xlabel(r"$R/(%.1f\ R_{\rm 200c})$" % r200c_scale)

plt.figure(40)
m2s, m1s, med, p1s, p2s = np.percentile(
    b_r90, [50-95/2, 50-68/2, 50, 50+68/2, 50+95/2], axis=0,
)
plt.fill_between(xs, m2s, p2s, facecolor="cyan", alpha=0.3)
plt.fill_between(xs, m1s, p1s, facecolor="cyan", alpha=0.3)
plt.plot(xs, p2s, c="c", lw=1)
plt.plot(xs, m2s, c="c", lw=1)
plt.plot(xs, p1s, c="c", lw=1)
plt.plot(xs, m1s, c="c", lw=1)
plt.plot(xs, med, c="c", lw=3)
plt.ylim(0, 1)
plt.plot([1, 1], [0, 1], "--k")
plt.title(r"$R/R_{90}$")
plt.ylabel(r"$N(<R/R_{90})/N_{\rm tot}$")
if r90_scale == 1.0:
    plt.xlabel(r"$R/R_{90}$")
else:
    plt.xlabel(r"$R/(%.1f\ R_{90})$" % r90_scale)


plt.figure(41)
plt.ylim(0, 6)
plt.hist(r85_r_sp, histtype="step", color="red", range=(-0.5, +0.5), bins=15,
         label=r"$R^*=R_{85}$", normed=True, lw=3)
plt.hist(r90_r_sp, histtype="step", color="blue", range=(-0.5, +0.5), bins=15,
         label=r"$R^*=R_{90}$", normed=True, lw=3)
plt.hist(r95_r_sp, histtype="step", color="green", range=(-0.5, +0.5), bins=15,
         label=r"$R^*=R_{95}$", normed=True, lw=3)
plt.legend(loc="upper left")
plt.xlabel(r"$R^*/R_{\rm s}$")
plt.ylabel(r"$\pi(R^*/R_{\rm sp})$")
lo, hi = plt.ylim()
plt.plot([0, 0], [lo, hi], "--k")


cmap = "afmhot"
plt.figure(42)
plt.hist2d(ds_r200m, xps_all, bins=40, range=[[0, 2], [0, 1]],
           cmap=cmap)
plt.title(r"$R_{\rm 200m}$")
plt.ylabel(r"$R(z=z_{\rm peri}) / R_{\rm 200m}(z=z_{\rm peri})$")
plt.xlabel(r"$R/R_{\rm 200m}$")
lo, hi = plt.ylim()
plt.plot([1, 1], [lo, hi], "w", lw=3)
plt.ylim(lo, hi)

plt.figure(43)
plt.hist2d(ds_r_sp, xps_all, bins=40, range=[[0, 2], [0, 1]],
           cmap=cmap)
plt.title(r"$R_{\rm sp}$")
plt.ylabel(r"$R(z=z_{\rm peri}) / R_{\rm 200m}(z=z_{\rm peri})$")
plt.xlabel(r"$R/R_{\rm sp}$")
plt.plot([1, 1], [lo, hi], "w", lw=3)
plt.ylim(lo, hi)

plt.figure(44)
plt.hist2d(ds_r_max, xps_all, bins=40, range=[[0, 2], [0, 1]],
           cmap=cmap)
plt.title(r"$R_{\rm max}$")
plt.ylabel(r"$R(z=z_{\rm peri}) / R_{\rm 200m}(z=z_{\rm peri})$")
plt.xlabel(r"$R/R_{\rm max}$")
plt.plot([1, 1], [lo, hi], "w", lw=3)
plt.ylim(lo, hi)

plt.show()
