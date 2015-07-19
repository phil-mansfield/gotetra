from __future__ import division

import abc
import sys
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
import numpy.random as rand
import scipy.interpolate as intr

class AbstractMetropolis(object):
    __metaclass__ = abc.ABCMeta

    def sample(self, steps, T, sample_rate=1, burn_in=0):
        if self.tracks_stats(): self.reset_stats()

        snapshots, accepts = [], 0
        for step in xrange(steps):
            updates = self.update_candidates()
            for update in updates:
                dE = self.dE(update)
                if dE < 0 or rand.random() < np.exp(-dE/T):
                    self.change(update)
                    accepts += 1
            if self.tracks_stats and burn_in <= step:
                self.update_stats()
            if sample_rate is not None and step % sample_rate == 0:
                snapshots.append(self.state())
        return np.array(snapshots), accepts

    @abc.abstractmethod
    def update_candidates(self): pass
    @abc.abstractmethod
    def dE(self, update): pass
    @abc.abstractmethod
    def change(self, update): pass
    @abc.abstractmethod
    def state(self): pass
    @abc.abstractmethod
    def tracks_stats(self): pass

class Potential1D(AbstractMetropolis):
    def __init__(self, f, x0, dx):
        self.f, self.x, self.dx = f, x0, dx
    def update_candidates(self):
        dx = self.dx * 2 * rand.random() - self.dx
        return [self.x + dx]
    def dE(self, update):
        return self.f(update) - self.f(self.x)
    def change(self, update):
        self.x = update
    def state(self):
        return self.x
    def tracks_stats(self):
        return False

def bin_by_thetas(rs, thetas, bins):
    pt_bins = [[] for _ in xrange(bins)]
    width = 2 * np.pi
    dtheta = width / bins
    for r, theta in zip(rs, thetas):
        idx = int((theta % width) / dtheta)
        pt_bins[idx].append(r)
    theta_bins = 2 * np.pi * (np.arange(bins) + 0.5) / bins
    return theta_bins, pt_bins

class Membrane(AbstractMetropolis):
    def __init__(self, pt_rs, pt_thetas, bins, dr=0.1):
        self.bins = bins
        self.dr = dr
        _, self.pt_bins = bin_by_thetas(pt_rs, pt_thetas, bins)

        for (i, pt_bin) in enumerate(self.pt_bins):
            self.pt_bins[i] = sorted(pt_bin)

        self.rs = np.ones(bins) * np.mean(pt_rs)
        for i, pts in enumerate(self.pt_bins):
            if len(pts) != 0: self.rs[i] = pts[0]
        self.weights = np.sqrt(map(len, self.pt_bins))
        self.weights /= max(self.weights)

        self.angle_edges = np.linspace(0, 2*np.pi, bins+1)
        self.reset_stats()

    def tracks_stats(self): return True

    def reset_stats(self):
        self.rs_sum = np.zeros(self.bins)
        self.rs_sqr = np.zeros(self.bins)
        self.trials = 0

    def update_stats(self):
        self.trials += 1
        self.rs_sum += self.rs
        self.rs_sqr += self.rs*self.rs

    def get_mean(self):
        return self.rs_sum / self.trials

    def get_std(self):
        return np.sqrt(self.rs_sqr - self.rs_mean**2)

    def update_candidates(self):
        idx = rand.randint(self.bins)
        dr = (rand.rand()*2 - 1) * self.dr
        return [(idx, dr + self.rs[idx])]

    def min_E(self, r_min, r, w):
        #return -np.exp(-(r-r_min)**2 / 4) * 4 * w
        return (r - r_min)**2

    def dE(self, update):
        i, r1 = update
        n = self.bins
        r0, rl, rr = self.rs[i], self.rs[(i+n-1)%n], self.rs[(i+1)%n]
        if len(self.pt_bins[i]) == 0:
            E0 = (r0 - rl)**2 + (r0 - rr)**2
            E1 = (r1 - rl)**2 + (r1 - rr)**2
        else:
            r_pt = self.pt_bins[i][0]
            w = self.weights[i]
            E0 = (r0 - rl)**2 + (r0 - rr)**2 + self.min_E(r0, r_pt, w)
            E1 = (r1 - rl)**2 + (r1 - rr)**2 + self.min_E(r1, r_pt, w)
        return E1 - E0

    def change(self, update):
        idx, r = update
        self.rs[idx] = r

    def state(self):
        return self.rs, self.angle_edges

def pretty_fig(n):
    plt.figure(n, figsize=(8, 8))
    plt.rc('text', usetex=True)
    plt.rc('font',size=19)
    plt.rc('xtick.major',pad=5); plt.rc('xtick.minor',pad=5)
    plt.rc('ytick.major',pad=5); plt.rc('ytick.minor',pad=5)

def plot_circle_segments(
    ax, rs, angle_edges, c="k", label=None, center=(0,0),
):
    angle_edges *= 360 / (2 * np.pi)
    for r, a_lo, a_hi in zip(rs, angle_edges[:-1], angle_edges[1:]):
        arc = patches.Arc(center, 2*r, 2*r, theta1=a_lo, theta2=a_hi,
                          color=c, linewidth=2)
        ax.add_patch(arc)

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

def plot_kde(theta, pt_bin, h, range, c="k"):
    low, high = range
    sp = gaussian_kde_spline(pt_bin, h, range)

    curve_rs = np.linspace(0, 1, 100)
    plt.plot(curve_rs, sp(curve_rs), lw=3, c=c)
    plt.plot(pt_bin, sp(pt_bin), "o", c=c, label=r"$\theta =%f.1$" % theta)

def plot_kde_max(fig, theta, pt_bin, h, range, c="w"):
    low, high = range
    sp = gaussian_kde_spline(pt_bin, h, range)
    curve_rs = np.linspace(0, 1, 100)
    max_r = curve_rs[np.argmax(sp(curve_rs))]
    x, y = max_r * np.cos(theta), max_r * np.sin(theta)
    if sp(max_r) > 4 and len(pt_bin) > 5:
        fig.plot(x, y, "o", c=c)

if __name__ == "__main__":
    file_name = sys.argv[1]

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 8, forward=True)
    pt_xs, pt_ys = read_pts(file_name)
    pt_rs, pt_thetas = to_radial(pt_xs, pt_ys)

    """
    T = 0.5
    for c in ["Orange", "r", "g", "b", "k"]:
        membrane = Membrane(pt_rs, pt_thetas, 25, dr=1)
        _, accepts = membrane.sample(200 * 1000, T, sample_rate=None)
        print "%d acceptances" % accepts
        rs, edges = membrane.state()
        rs_mean = membrane.get_mean()
        plot_circle_segments(ax, rs_mean, edges, c=c)
    """
    plt.plot(pt_xs, pt_ys, "ow")
    ylo, yhi = plt.ylim()
    xlo, xhi = plt.xlim()
    dist = max(max(yhi, xhi), -min(ylo, xlo))
    plt.xlim((-dist, dist))
    plt.ylim((dist, -dist))

    theta_bins, pt_bins = bin_by_thetas(pt_rs, pt_thetas, 25)

    cs = ["Red", "Blue",
          "Chartreuse", "Cyan",
          "DarkOrange", "Magenta",
          "Black", "BlueViolet",
          "DarkGoldenRod", "DarkRed"]

    pretty_fig(2)
    h = 0.05
    for theta, pt_bin, c in zip(theta_bins, pt_bins, cs):
        plot_kde_max(ax, theta, pt_bin, h, (0, 1), c=c)
        plot_kde(theta, pt_bin, h, (0, 1), c=c)
    plt.legend(fontsize=12, loc="upper left")

    for theta, pt_bin in zip(theta_bins, pt_bins)[len(cs):]:
        plot_kde_max(ax, theta, pt_bin, h, (0, 1), c="k")

    plt.show()
