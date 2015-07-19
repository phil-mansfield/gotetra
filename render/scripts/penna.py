from __future__ import division

import abc
import sys
import matplotlib.pyplot as plt
import scipy.optimize as opt
from matplotlib import patches
import numpy as np
import numpy.random as rand
import numpy.linalg as linalg


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

def fast_penna_coeffs(n, xs, ys, zs):
    N = 2*n + 1
    rs = np.sqrt(xs**2 + ys**2 + zs**2)
    rNs = rs**(N+1)
    M = []
    for k in xrange(2):
        zks = zs**k
        for j in xrange(n+1):
            yjs = ys**j
            for i in xrange(n+1):
                xis = xs**i
                M.append(rs**(N-i-j-k)*xis*yjs*zks)
    M = np.matrix(M)
    return np.array(rNs * linalg.pinv(M))[0]

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

def find_max_r(rs):
    vals, edges = np.histogram(rs, bins=2*(len(rs))**0.33)
    bin_rs = (edges[1:] + edges[:-1]) / 2

    max = bin_rs[np.argmax(vals)]
    dr = edges[1] - edges[0]
    for (i, r) in enumerate(bin_rs):
        if r > max * 1.3:
            low_edges, high_edges = edges[:i], edges[i:-1]
            low_vals, high_vals = vals[:i], vals[i:]
            break
    pretty_fig(7)
    plt.bar(low_edges, low_vals, width=dr, color="r")
    plt.bar(high_edges, high_vals, width=dr, color="w")
    plt.xlabel("$r$ [Mpc/$h$]")
    plt.ylabel("Counts")

    return bin_rs[np.argmax(vals)]

def filter_around_max(xs, ys, zs, eta=0.3):
    rs = np.sqrt(xs**2 + ys**2 + zs**2)
    max_r = find_max_r(rs)
    valid = (rs < max_r*(1+eta)) & (rs > max_r*(1-eta))
    return max_r, xs[valid], ys[valid], zs[valid]

def filter_around(vals, filter_vals, target, eta=0.3):
    valid = (filter_vals < target*(1+eta)) & (filter_vals > target*(1-eta))
    return vals[valid]

def chi2(func, phis, thetas, rs, n):
    deltas = func(phis, thetas) - rs
    return np.sum(deltas**2) / (len(rs) - n - 1) / np.mean(rs)**2

def multi_fit(xs, ys, zs, r_lims, n):
    rs = np.sqrt(xs**2 + ys**2 + zs**2)
    for lim in r_lims:
        lim_xs, lim_ys, lim_zs = xs[rs < lim], ys[rs < lim], zs[rs < lim]
        lim_rs = rs[rs < lim]
        print "R limit = %g (%d points)" % (lim, len(lim_xs))
        if (n+1)*(n+1)*2 >= len(lim_xs): continue
        print len(lim_xs)
        c_ijks = penna_coeff(n, lim_xs, lim_ys, lim_zs)
        phis = np.arctan2(lim_ys, lim_xs)
        thetas = np.arccos(lim_zs / lim_rs)
        print "Chi^2 = %g" % chi2(penna_func(c_ijks, n), phis,
                                  thetas, lim_rs, n)

if __name__ == "__main__":
    rs_file_name = sys.argv[1]
    plane_file_name = sys.argv[2]

    eta = 0.3

    raw_xs, raw_ys, raw_zs = read_pts(rs_file_name)
    max_r, xs, ys, zs = filter_around_max(raw_xs, raw_ys, raw_zs, eta=eta)

    raw_p_x0s, raw_p_x1s, axis = read_plane(plane_file_name)
    p_rs = np.sqrt(raw_p_x0s**2 + raw_p_x1s**2)
    p_x0s = filter_around(raw_p_x0s, p_rs, max_r, eta=eta)
    p_x1s = filter_around(raw_p_x1s, p_rs, max_r, eta=eta)

    pretty_fig(0)
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

    plot_func(lambda a, b: max_r*(1+eta), axis, c="r", lw=1)
    plot_func(lambda a, b: max_r*(1-eta), axis, c="r", lw=1)
    plt.plot(raw_p_x0s, raw_p_x1s, "wo")
    plt.plot(p_x0s, p_x1s, "ro")

    furthest_r = np.max(np.sqrt(xs**2 + ys**2 + zs**2))
    multi_fit(xs, ys, zs, np.linspace(0, furthest_r, 20), 4)

    N = 4
    c_ijks = penna_coeff(N, xs, ys, zs)
    plot_func(penna_func(c_ijks, N), axis)

    ylo, yhi = plt.ylim()
    xlo, xhi = plt.xlim()
    dist = max(max(yhi, xhi), -min(ylo, xlo))
    plt.xlim((-dist, dist))
    plt.ylim((dist, -dist))
    plt.show()
