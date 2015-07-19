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

def gaussian_kde_spline(xs, h, range, bins=100):
    low, high = range
    eval_xs = np.linspace(low, high, bins)
    vals = np.zeros_like(eval_xs)
    for x in xs:
        vals += np.exp(-((x - eval_xs) / h)**2)
    vals /= np.sqrt(np.pi * 2) * len(xs) * h
    return intr.UnivariateSpline(eval_xs, vals, s=0)

def plot_kde(ipix, pt_bin, h, range, c="k"):
    low, high = range
    sp = gaussian_kde_spline(pt_bin, h, range)

    curve_rs = np.linspace(low, high, 100)
    plt.plot(curve_rs, sp(curve_rs), lw=3, c=c)
    plt.plot(pt_bin, sp(pt_bin), "o", c=c, label=r"ipix = %d" % ipix)

def kde_peak_heights(ipix, r_bins, h, range):
    low, high = range
    max_rs, max_heights = [], []
    for r_bin in r_bins:
        sp = gaussian_kde_spline(r_bin, h, range)
        curve_rs = np.linspace(low, high, 100)
        maxi = np.argmax(sp(curve_rs))
        max_rs.append(curve_rs[maxi])
        max_heights.append(sp(curve_rs)[maxi])
    return np.array(max_rs), np.array(max_heights)

def read_row(row):
    ts = [t for t in row.split(" ") if t != ""]
    return [float(ts[0]), float(ts[1]), float(ts[2]), float(ts[3])]

def read_pts(file_name):
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


def plot_func(f, axis, c="k", pts=200, lw=3, label=None):
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
    if label is None:
        plt.plot(x0s*rs, x1s*rs, lw=lw, c=c)
    else:
        plt.plot(x0s*rs, x1s*rs, lw=lw, c=c, label=label)

def penna_coeff(n, xs, ys, zs):
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

def healpix_bins(nside, rs, phis, thetas):
    npix = pf.nside2npix(nside)
    ipix = pf.ang2pix(nside, thetas, phis)
    r_bins = [[] for _ in xrange(npix)]
    for (r, i) in zip(rs, ipix):
        r_bins[i].append(r)
    return r_bins

def mode_radius(rs):
    vals, edges = np.histogram(rs, bins=2*(len(rs))**0.33)
    bin_rs = (edges[1:] + edges[:-1]) / 2
    return bin_rs[np.argmax(vals)]

def chi2(vals, fit_vals, nu):
    return np.sum((vals - fit_vals)**2 /
                  (len(vals) - 1 - nu))

def multi_fit(xs, ys, zs, r_lims, N=2):
    rs = np.sqrt(xs**2 + ys**2 + zs**2)
    chi2s, coeffs = [], []
    for lim in r_lims:
        valid = rs < lim
        lim_xs, lim_ys, lim_zs = xs[valid], ys[valid], zs[valid]
        lim_rs = rs[valid]
        nu = (N+1)**2 * 2
        if len(lim_rs) <= (N+1)**2 * 2: continue
        c_ijks = penna_coeff(N, lim_xs, lim_ys, lim_zs)
        func = penna_func(c_ijks, N)

        lim_thetas, lim_phis = pf.vec2ang(
            np.array(zip(lim_xs, lim_ys, lim_zs))
        )

        chi2s.append(chi2(lim_rs, func(lim_phis, lim_thetas), nu))
        coeffs.append(c_ijks)
    return coeffs, chi2s

if __name__ == "__main__":
    file_dir = sys.argv[1]
    halo_id = int(sys.argv[2])

    rs_fname = path.join(file_dir, "%d_caustic.txt" % halo_id)
    x_plane_fname = path.join(file_dir, "%dh_x_caustic.txt" % halo_id)
    y_plane_fname = path.join(file_dir, "%dh_y_caustic.txt" % halo_id)
    z_plane_fname = path.join(file_dir, "%dh_z_caustic.txt" % halo_id)

    xs, ys, zs = read_pts(rs_fname)
    rs = np.sqrt(xs**2 + ys**2 + zs**2)
    vecs = np.array(zip(xs, ys, zs))
    thetas, phis = pf.vec2ang(vecs)

    nside = 2
    npix = pf.nside2npix(nside)

    r_bins = healpix_bins(nside, rs, phis, thetas)
    h = mode_radius(rs) / 5
    max_r = np.max(rs)

    cs = ["Red", "Blue",
          "Green", "Black",
          "DarkOrange", "Magenta",
          "Cyan", "BlueViolet",
          "DarkGoldenRod", "DarkRed"]

    plot_bins = random.sample(zip(np.arange(len(r_bins)), r_bins), 8)

    pretty_fig(3)
    for i, (ipix, bin) in enumerate(plot_bins):
        plot_kde(ipix, bin, h, (0, max_r), c=cs[i])
    plt.xlabel(r"$r$ [Mpc/$h$]")
    plt.ylabel(r"KDE")
    plt.title(r"KDE Profiles")
    plt.legend(fontsize=12)

    pretty_fig(4)
    peak_rs, peak_heights = kde_peak_heights(ipix, r_bins, h, (0, max_r))
    plt.plot(peak_rs, peak_heights, "wo")
    r_cut = 1.3 * mode_radius(rs)
    plt.plot(peak_rs[peak_rs < r_cut], peak_heights[peak_rs < r_cut], "ro")
    plt.xlabel(r"$r_{\rm KDE}$ [Mpc/$h$]")
    plt.ylabel(r"$h_{\rm KDE}$")
    plt.title(r"KDE Peak Locations")

    valid_peaks = peak_heights > 0

    N = 2
    ipix = np.arange(npix)[valid_peaks]
    peak_rs = peak_rs[valid_peaks]
    kde_xs, kde_ys, kde_zs = pf.pix2vec(nside, ipix) * peak_rs
    c_ijks = penna_coeff(N, kde_xs, kde_ys, kde_zs)

    cut_ipix = ipix[peak_rs < r_cut]
    cut_peak_rs = peak_rs[peak_rs < r_cut]
    cut_kde_xs,cut_kde_ys,cut_kde_zs = pf.pix2vec(nside,cut_ipix)*cut_peak_rs
    cut_cijks = penna_coeff(N, cut_kde_xs, cut_kde_ys, cut_kde_zs)

    cut_raw_cijks = penna_coeff(N, xs[rs<r_cut], ys[rs<r_cut], zs[rs<r_cut])
    raw_cijks = penna_coeff(N, xs, ys, zs)


    for n, plane_fname in enumerate([x_plane_fname,
                                     y_plane_fname,
                                     z_plane_fname]):
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

        p_rs = np.sqrt(p_x0s**2 + p_x1s**2)
        plt.plot(p_x0s, p_x1s, "wo")
        plt.plot(0, 0, "kx")

        plot_func(penna_func(raw_cijks, N), axis, c="k", lw=1,
                  label="Full")
        plot_func(penna_func(cut_raw_cijks, N), axis, c="r", lw=1,
                  label=r"Full + $\eta$")
        plot_func(penna_func(c_ijks, N), axis, c="k",
                  label="KDE")
        plot_func(penna_func(cut_cijks, N), axis, c="r",
                  label=r"KDE + $\eta$")

        ylo, yhi = plt.ylim()
        xlo, xhi = plt.xlim()
        dist = max(max(yhi, xhi), -min(ylo, xlo))
        plt.xlim((-dist, dist))
        plt.ylim((dist, -dist))
        plt.legend(fontsize=12)

    plt.show()
