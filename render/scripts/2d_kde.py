from __future__ import division

import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as intr
import healpy as hp
import healpy.pixelfunc as pf
import random
import os.path as path
import scipy.optimize as opt

def bin_by_thetas(rs, thetas, bins):
    pt_bins = [[] for _ in xrange(bins)]
    width = 2 * np.pi
    dtheta = width / bins
    for r, theta in zip(rs, thetas):
        idx = int((theta % width) / dtheta)
        pt_bins[idx].append(r)
    theta_bins = 2 * np.pi * (np.arange(bins) + 0.5) / bins
    return theta_bins, pt_bins

def gaussian_kde_spline(xs, h, range, bins=100):
    low, high = range
    eval_xs = np.linspace(low, high, bins)
    vals = np.zeros_like(eval_xs)
    for x in xs:
        vals += np.exp(-((x - eval_xs) / h)**2)
    vals /= np.sqrt(np.pi * 2) * len(xs) * h
    return intr.UnivariateSpline(eval_xs, vals, s=0)

def argmaxes(xs):
    is_max = (xs[1:-1] > xs[2:]) & (xs[1:-1] > xs[:-2])
    return np.arange(1, len(xs) - 1)[is_max]

class KDETree(object):
    _h_factor = 10
    _base_spline_r = False

    def __init__(self, rs, thetas, splits):
        max_r = np.max(rs)
        self._h = max_r / self._h_factor
        self._range = (0, max_r)

        self._sp_tree = [[
            gaussian_kde_spline(rs, self._h, self._range)
        ]]
        self._th_tree = [[np.pi]]

        self._grow_trees(rs, thetas, splits)
        self._sp_rs = np.linspace(0, max_r, 100)
        self._find_maxes()
        self._connect_maxes()

    def _grow_trees(self, rs, thetas, splits):
        for split in xrange(splits):
            bins = 2**(1+split)
            th_bins, r_bins = bin_by_thetas(rs, thetas, bins)

            sps = [None] * bins
            for (i, r_bin) in enumerate(r_bins):
                sps[i] = gaussian_kde_spline(r_bin, self._h, self._range)
            self._th_tree.append(th_bins)
            self._sp_tree.append(sps)
        self._th_tree = np.array(self._th_tree)
        self._sp_tree = np.array(self._sp_tree)

    def _find_maxes(self):
        self._maxes_tree = []
        for sps in self._sp_tree:
            maxes = [self._sp_rs[argmaxes(sp(self._sp_rs))] for sp in sps]
            self._maxes_tree.append(maxes)
        self._maxes_tree = np.array(self._maxes_tree)

    def _connect_maxes(self):
        self._conn_maxes = [[self._maxes_tree[0][0][0]]]
        for split, maxes in enumerate(self._maxes_tree[1:]):
            prev_maxes = self._conn_maxes[-1]
            curr_maxes = [None] * (len(prev_maxes) * 2)
            for node in xrange(len(maxes)):
                node_maxes = maxes[node]
                if len(maxes) == 0: continue
                node_prev_max = prev_maxes[node // 2]
                if node_prev_max is None:
                    conn_max = None
                else:
                    conn_idx = np.argmin(np.abs(node_maxes - node_prev_max))
                    conn_max = node_maxes[conn_idx]

                if (conn_max is not None and
                    (np.abs(conn_max - node_prev_max) < self._h or
                     split == 0)):
                    curr_maxes[node] = conn_max
                else:
                    sp = self.get_spline(split)
                    sp_r = sp(self._th_tree[split+1][node])
                    for max in node_maxes:
                        if np.abs(max - sp_r) < self._h:
                            curr_maxes[node] = max

            self._conn_maxes.append(curr_maxes)


    def plot(self, top_level_cs):
        plt.ylabel(r"r [Mpc/$h$]")
        plt.xlabel(r"$\theta$")
        plt.ylim(self._range)

        # Plot lines
        for row in xrange(len(self._th_tree)):
            for col in xrange((len(self._th_tree[row]))):
                th = self._th_tree[row][col]
                if row == 0: continue
                prev_th = self._th_tree[row-1][col//2]
                prev_r = self._conn_maxes[row-1][col//2]
                if prev_r is None: continue

                for r in self._maxes_tree[row][col]:
                    plt.plot([prev_th, th], [prev_r, r], "k--")

                conn_r = self._conn_maxes[row][col]
                if conn_r is None: continue
                plt.plot([prev_th, th], [prev_r, conn_r], "k", lw=2)

        # Plots points
        for row in xrange(len(self._th_tree)):
            for col in xrange((len(self._th_tree[row]))):
                th = self._th_tree[row][col]
                for r in self._maxes_tree[row][col]:
                    if row == len(self._th_tree) - 1:
                        plt.plot(th, r, "x", c=cs[col // 2])
                    elif row == len(self._th_tree) - 2:
                        plt.plot(th, r, "o", c=cs[col])
                    else:
                        plt.plot(th, r, "wo")


    def _get_finest_max(self, idx, level):
        for i, l in enumerate(xrange(level, -1, -1)):
            r = self._conn_maxes[l][idx // (2**i)]
            if r is not None: return r
        raise Exception("Impossible")

    def get_conn_maxes(self, level):
        ths = self._th_tree[level]
        maxes = self._conn_maxes[level]

        ret_maxes = [0] * len(maxes)
        for i, max in enumerate(maxes):
            ret_maxes[i] = self._get_finest_max(i, level)

        return np.array(ret_maxes), ths

    def get_spline(self, level):
        if self._base_spline_r:
            maxes, ths = self.get_conn_maxes(level)
            sp_ths = np.append(np.append(ths - 2*np.pi, ths), ths + 2*np.pi)
            sp_maxes = np.append(np.append(maxes, maxes), maxes)
            r_sp = intr.UnivariateSpline(sp_ths, sp_maxes, s=0)
        else:
            maxes, ths = self.get_conn_maxes(level)
            sp_ths = np.append(np.append(ths - 2*np.pi, ths), ths + 2*np.pi)
            sp_maxes = np.append(np.append(maxes, maxes), maxes)
            sp_xs, sp_ys = to_cartesian(sp_maxes, sp_ths)
            x_sp = intr.UnivariateSpline(sp_ths, sp_xs, s=0)
            y_sp = intr.UnivariateSpline(sp_ths, sp_ys, s=0)
            # Shhhhh... not really a spline...
            def r_sp(ths):
                xs, ys = x_sp(ths), y_sp(ths)
                rs, _ = to_radial(xs, ys)
                return rs
        return r_sp

    def filter_nearby(self, rs, thetas, level, dr=None):
        if dr is None: dr = self._h
        sp = self.get_spline(level)
        valid = np.abs(rs - sp(thetas)) < dr
        return rs[valid], thetas[valid], np.arange(len(valid))[valid]

def pretty_fig(n):
    plt.figure(n, figsize=(8, 8))
    plt.rc('text', usetex=True)
    plt.rc('font',size=19)
    plt.rc('xtick.major',pad=5); plt.rc('xtick.minor',pad=5)
    plt.rc('ytick.major',pad=5); plt.rc('ytick.minor',pad=5)

def to_radial(xs, ys):
    rs = np.sqrt(xs**2 + ys**2)
    thetas = np.arctan2(ys, xs)
    return rs, thetas

def to_cartesian(rs, thetas):
    ys, xs = np.sin(thetas) * rs, np.cos(thetas) * rs
    return xs, ys

def read_row(row):
    ts = [t for t in row.split(" ") if t != ""]
    return [float(ts[0]), float(ts[1]), float(ts[2]), float(ts[3])]

def read_pts(file_name):
    with open(file_name) as fp: s = fp.read()
    lines = s.split("\n")
    caustic_rows = map(read_row, lines)
    c_xs, c_ys, c_zs, c_rs = map(np.array, zip(*caustic_rows))

    if np.all(c_xs == 0):
        return c_ys * c_rs, c_zs * c_rs
    if np.all(c_ys == 0):
        return c_xs * c_rs, c_zs * c_rs
    if np.all(c_zs == 0):
        return c_xs * c_rs, c_ys * c_rs
    assert 0

def plot_func(f, axis, c="k", pts=200, lw=3):
    theta0s = np.linspace(0, np.pi * 2, pts)
    x0s, x1s = np.cos(theta0s), np.sin(theta0s)
    if axis == 0:
        xs, ys, zs = np.zeros(pts), x0s, x1s
    elif axis == 1:
        xs, ys, zs = x0s, np.zeros(pts), x1s
    else:
        xs, ys, zs = x0s, x1s, np.zeros(pts)
    phis = np.arctan2(ys, xs)
    thetas = np.arccos(zs, np.ones(pts))

    rs = f(phis, thetas)
    plt.plot(x0s*rs, x1s*rs, lw=lw, c=c)

def penna_coeffs(n, xs, ys, zs):
    phis = np.arctan2(ys, xs)
    rs = np.sqrt(xs**2 + ys**2 + zs**2)
    thetas = np.arccos(zs / rs)

    trig_coeffs = []
    for k in xrange(2):
        cos_k = np.cos(thetas)**k
        for j in xrange(0, n+1):
            sin_j = np.sin(phis)**j
            for i in xrange(0, n+1):
                cos_i = np.cos(phis)**i
                sin_ij = np.sin(thetas)**(i+j)
                trig_coeffs.append(sin_ij * cos_k * sin_j * cos_i)
    trig_coeffs = np.array(trig_coeffs)

    def fit_func(p, phis, thetas):
        sum = np.zeros(len(phis))
        for idx, coeffs in enumerate(trig_coeffs):
            sum += p[idx] * coeffs
        return sum

    def err_func(p, phis, thetas, rs): return rs - fit_func(p, phis, thetas)
    p_init = np.ones((n+1)*(n+1)*2)
    coeffs, _ = opt.leastsq(err_func, p_init, args=(phis, thetas, rs))
    return coeffs

def penna_func(c_ijks, N):
    def func(phis, thetas):
        trig_coeffs = []
        idx = 0
        for k in xrange(2):
            cos_k = np.cos(thetas)**k
            for j in xrange(0, N+1):
                sin_j = np.sin(phis)**j
                for i in xrange(0, N+1):
                    cos_i = np.cos(phis)**i
                    sin_ij = np.sin(thetas)**(i+j)
                    trig_coeffs.append(sin_ij * cos_k * sin_j * cos_i)
        trig_coeffs = np.array(trig_coeffs)

        sum = np.zeros(len(phis))
        for idx, coeffs in enumerate(trig_coeffs):
            sum += c_ijks[idx] * coeffs
        return sum
    return func

def read_row(row):
    ts = [t for t in row.split(" ") if t != ""]
    return [float(ts[0]), float(ts[1]), float(ts[2]), float(ts[3])]

def read_full(file_name):
    with open(file_name) as fp: s = fp.read()
    lines = s.split("\n")
    caustic_rows = map(read_row, lines)
    xs, ys, zs, rs = map(np.array, zip(*caustic_rows))
    return xs*rs, ys*rs, zs*rs

def read_plane(file_name):
    with open(file_name) as fp: s = fp.read()
    lines = s.split("\n")
    caustic_rows = map(read_row, lines)
    xs, ys, zs, rs = map(np.array, zip(*caustic_rows))
    if np.all(xs == 0):
        return ys*rs, zs*rs, 0
    elif np.all(ys == 0):
        return xs*rs, zs*rs, 1
    else:
        return xs*rs, ys*rs, 2

if __name__ == "__main__":
    file_dir = sys.argv[1]
    halo_id = int(sys.argv[2])

    x_plane_fname = path.join(file_dir, "%dh_x_caustic.txt" % halo_id)
    y_plane_fname = path.join(file_dir, "%dh_y_caustic.txt" % halo_id)
    z_plane_fname = path.join(file_dir, "%dh_z_caustic.txt" % halo_id)

    cs = ["Red", "Blue",
          "Green", "Black",
          "DarkOrange", "Magenta",
          "Cyan", "BlueViolet",
          "DarkGoldenRod", "DarkRed"]


    plane_files = [x_plane_fname, y_plane_fname, z_plane_fname]
    shell_idxs = []
    for n, plane_fname in enumerate(plane_files):
        pretty_fig(n)
        p_x0s, p_x1s, axis = read_plane(plane_fname)
        if axis == 0:
            plt.title("X Slice")
            plt.ylabel("Z [Mpc/$h$]")
            plt.xlabel("Y [Mpc/$h$]")
        elif axis == 1:
            plt.title("Y Slice")
            plt.ylabel("Z [Mpc/$h$]")
            plt.xlabel("X [Mpc/$h$]")
        else:
            plt.title("Z Slice")
            plt.ylabel("Y [Mpc/$h$]")
            plt.xlabel("X [Mpc/$h$]")
        plt.plot(p_x0s, p_x1s, "wo")

        plt.plot(p_x0s, p_x1s, "wo")
        plt.plot(0, 0, "kx")

        ylo, yhi = plt.ylim()
        xlo, xhi = plt.xlim()
        dist = max(max(yhi, xhi), -min(ylo, xlo))
        plt.xlim((-dist, dist))
        plt.ylim((dist, -dist))

        p_rs, p_thetas = to_radial(p_x0s, p_x1s)
        max_r = np.max(p_rs)
        h = max_r / 10

        tree = KDETree(p_rs, p_thetas, 4)

        # Plot nearby points
        dr = max_r / 20
        near_rs, near_ths, near_idxs = tree.filter_nearby(
            p_rs, p_thetas, 4, dr=dr,
        )
        shell_idxs.append(near_idxs)

        near_x0s, near_x1s = to_cartesian(near_rs, near_ths)
        plt.plot(near_x0s, near_x1s, "ro")

        # Plot tree leaves.
        x_rs, x_ths = tree.get_conn_maxes(4)
        x_x0s, x_x1s = to_cartesian(x_rs, x_ths)
        for i, (x_x0, x_x1) in enumerate(zip(x_x0s, x_x1s)):
            plt.plot(x_x0, x_x1, "x", c=cs[i // 2])

        # Plot lower tree nodes
        o_rs, o_ths = tree.get_conn_maxes(3)
        o_x0s, o_x1s = to_cartesian(o_rs, o_ths)
        for i, (o_x0, o_x1) in enumerate(zip(o_x0s, o_x1s)):
            plt.plot(o_x0, o_x1, "o", c=cs[i])

        """
        # Plot enclosing spline
        sp = tree.get_spline(4)
        sp_ths = np.linspace(0, 2*np.pi, 100)
        sp_rs = sp(sp_ths)
        sp_x0s, sp_x1s = to_cartesian(sp_rs, sp_ths)
        plt.plot(sp_x0s, sp_x1s, c="k")
        """

        pretty_fig(n + 3)
        if axis == 0: plt.title("X Slice KDE Tree")
        if axis == 1: plt.title("Y Slice KDE Tree")
        if axis == 2: plt.title("Z Slice KDE Tree")
        tree.plot(cs)

    # Plot Penna fit.
    p_xs, p_ys, p_zs = [], [], []
    for n, plane_file in enumerate(plane_files):
        xs, ys, zs = read_full(plane_file)
        valid = shell_idxs[n]
        p_xs = np.append(p_xs, xs[valid])
        p_ys = np.append(p_ys, ys[valid])
        p_zs = np.append(p_zs, zs[valid])

    N = 3
    c_ijks = penna_coeffs(N, p_xs, p_ys, p_zs)
    fit = penna_func(c_ijks, N)
    for n in xrange(3):
        plt.figure(n)
        plot_func(fit, n)

    plt.show()
