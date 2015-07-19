from __future__ import division

import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import scipy.stats as stats
import scipy.interpolate as intr
import deriv

class SplinedBinConverter(object):
    def __init__(self, low, high, bins):
        self.low, self.high, self.bins = low, high, bins
        dx = (self.high - self.low) / bins
        self.edges = np.linspace(low, high, bins + 1)
    
    def dx(self):
        return (self.high - self.low) / self.bins

    def convert(self, xs, ys):
        ys = ys[(self.edges[0] <= xs) & (xs <= self.edges[-1])]
        xs = xs[(self.edges[0] <= xs) & (xs <= self.edges[-1])]
        if len(xs) <= 3 or len(ys) <= 3: return None, None

        low_edge_idx = np.searchsorted(self.edges, xs[0])
        if low_edge_idx == 0: low_edge_idx = 1
        high_edge_idx = np.searchsorted(self.edges, xs[-1])

        sp = intr.UnivariateSpline(xs, ys, s=0)

        if high_edge_idx == low_edge_idx:
            return (np.array([xs[0], xs[-1]]),
                    np.array([sp.integral(xs[0], xs[-1]) / (xs[-1] - xs[0])]))

        edges = self.edges[low_edge_idx - 1: high_edge_idx + 1]

        first = self._first_bin(edges, xs, sp)
        mid = self._mid_bins(edges, xs, sp)
        last = self._last_bin(edges, xs, sp)
        
        return edges, append(first, mid, last)
        
    def _first_bin(self, edges, xs, sp):
        if xs[0] == edges[1]: return []
        return [sp.integral(xs[0], edges[1]) / (edges[1] - xs[0])]

    def _mid_bins(self, edges, xs, sp):
        vals = np.zeros(len(edges) - 3)
        for i in xrange(len(vals)):
            start_edge, end_edge = edges[i + 1], edges[i + 2]
            vals[i] = (sp.integral(start_edge, end_edge) /
                       (end_edge - start_edge))
        return vals

    def _last_bin(self, edges, xs, sp):
        if xs[-1] == edges[-2]: return []
        return [sp.integral(edges[-2], xs[-1]) / (xs[-1] - edges[-2])]

def append(*arrays):
    out, idx = np.zeros(sum(map(len, arrays))), 0
    for array in arrays:
        for i in xrange(len(array)):
            out[idx] = array[i]
            idx += 1
    return out

def pretty_fig(n):
    """ pretty_fig(n) is equivalent to plt.figure(n), except it also
    sets a number of options which make the resulting plot look nicer.
    """
    plt.figure(n)
    plt.rc('text', usetex=True)
    plt.rc('font',size=19)
    plt.rc('xtick.major',pad=5); plt.rc('xtick.minor',pad=5)
    plt.rc('ytick.major',pad=5); plt.rc('ytick.minor',pad=5)

def nan_split(rs, rhos):
    """ nan_split(rs, rhos) splits up rs and rhos into lists of contiguous
    non-NaN seqeunces of values.
    """
    rs_group, rhos_group = [], []
    start_i = 0
    prev_nan = False

    for i in xrange(len(rs)):
        if np.isnan(rhos[i]):
            if not prev_nan:
                if i != start_i:
                    rs_group.append(np.array(rs[start_i: i]))
                    rhos_group.append(np.array(rhos[start_i: i]))
            prev_nan = True
        else:
            if prev_nan:
                start_i = i
            prev_nan = False
    if not prev_nan:
        rs_group.append(np.array(rs[start_i:len(rs)]))
        rhos_group.append(np.array(rhos[start_i:len(rhos)]))
    return rs_group, rhos_group

def r_sp(rs, rhos, derivs, lim=-5):
    curr_min = rhos <= np.minimum.accumulate(rhos)
    idxs = signal.argrelextrema(derivs, np.less)[0]
    idxs = np.array([idx for idx in idxs if idx != 0 and idx != len(rs) - 1])
    if len(idxs) == 0: return np.nan
    idxs = idxs[curr_min[idxs]]
    if len(idxs) == 0: return np.nan
    idxs = idxs[derivs[idxs] < lim]
    if len(idxs) == 0: return np.nan
    min_idx = idxs[np.argmin(derivs[idxs])]    
    return rs[min_idx]

# rs must be evenly spaced.
def find_splashback(rs, rhos, sg_window):
    sm_rhos = signal.savgol_filter(rhos, sg_window, 4)
    dr = (rs[-1] - rs[0]) / (len(rs) - 1)
    sm_derivs = signal.savgol_filter(rhos, sg_window, 4, deriv=1, delta=dr)
    return r_sp(rs, sm_rhos, sm_derivs, lim=-5)

def splashback_range(rs, rhos, range_max=101):
    if rs[0] <= 0: rs, rhos = rs[1:], rhos[1:]
    r_sps = []
    lrs, lrhos = np.log10(rs), np.log10(rhos)
    for sg_window in xrange(2, range_max//2):
        sg_window = sg_window * 2 + 1
        if sg_window >= len(lrs): break
        r_sp = find_splashback(lrs, lrhos, sg_window)
        if np.isnan(r_sp): break
        r_sps.append(r_sp)
    return 10**np.array(r_sps)

cs = ["k", "r", "b", "g", "pink", "orange",
      "brown", "y", "DarkMagenta", "LightSlateGray"]

if __name__ ==  "__main__":
    rows = np.loadtxt(sys.argv[1])
    vec_xs, vec_ys, vec_zs = rows[:3]
    idxs = np.arange(len(vec_xs) // 2) * 2
    vec_xs, vec_ys, vec_zs = vec_xs[idxs], vec_ys[idxs], vec_zs[idxs]
    rows = rows[3:]
    cols = map(lambda *x: x, *rows)
    profile_count = len(cols) // 2

    bins = 200
    range_max = 100
    log_low, log_high = -2, 0
    log_converter = SplinedBinConverter(log_low, log_high, bins)
    
    max_cutoff, std_cutoff = 0.1, 0.1

    n, m = 0, 0

    maxes, stds, valid_rs = [], [], []
    for i in xrange(profile_count):
        rs_group, rhos_group = nan_split(cols[2*i], cols[2*i + 1])
        for rs, rhos in zip(rs_group, rhos_group):
            R = np.nan
            if rs[0] <= 0: rs, rhos = rs[1:], rhos[1:]
            edges, vals = log_converter.convert(np.log10(rs), np.log10(rhos))
            if edges is None: continue 
            rs, rhos = 10**((edges[:-1] + edges[1:]) / 2), 10**vals

            if len(rs) <= 21: continue
            rs_range = splashback_range(rs, rhos, range_max)
            if len(rs_range) * 2 <= 21: continue

            r_mean, r_std = np.mean(rs_range), np.std(rs_range)
            drs = np.abs(rs_range[1:] - rs_range[:-1])
            dr_max, dr_max_mean = np.max(drs), np.mean(drs)
            dr_sign_max_mean = np.mean(rs_range[1:] - rs_range[:-1])
            rs_range_diff = np.max(rs_range) - np.min(rs_range)

            # Howwwwwww??
            if dr_max == 0 or r_std == 0: continue
            m += 1

            maxes.append(dr_max / r_mean)
            stds.append(r_std / r_mean)

            # Figure 0 
            plt.figure(0)
            c = cs[i] if i < len(cs) else "w"
            plt.plot(dr_max / r_mean, r_std / r_mean, "o")

            if i < 10:
                plt.figure(4)
                windows = np.arange(2, len(rs_range) + 2)*2 + 1
                plt.plot(windows, rs_range, c=c, lw=3)
                is_good = dr_max/r_mean<max_cutoff and r_std/r_mean<std_cutoff
                plt.plot(windows, rs_range, "o",  c=c, lw=3,
                         label="%.3f %.3f %.3f %.5f %s" % 
                         (dr_max / r_mean, r_std / r_mean,
                          rs_range_diff / r_mean,
                          np.abs(dr_sign_max_mean / r_mean),
                          "*" if is_good else ""))

            if dr_max / r_mean < max_cutoff and r_std / r_mean < std_cutoff:
                n += 1
                R = rs_range[0]
                valid_rs.append(R)
                print ("%9.4g %9.4g %9.4g %9.4g %9d" % 
                       (vec_xs[i], vec_ys[i], vec_zs[i], R, i))

    sys.stderr.write("%d LoS's found, %d are good\n" % (m, n))

    # Figure 1
    plt.figure(1)
    plt.hist(np.log10(maxes), bins=40, range=(-2, 0))
    
    # Figure 2
    plt.figure(2)
    plt.hist(np.log10(stds), bins=40, range=(-2, 0))
    
    # Figure 3
    plt.figure(3)
    plt.hist(valid_rs, bins=40, range=(0, 1))

########################
# Do all the plotting. #
########################

# Figure 0
pretty_fig(0)
plt.ylabel(r"stdev($R_{\rm sp}$) [Mpc/$h$]")
plt.ylim((1e-2, 1))
plt.yscale("log")
plt.xlabel(r"max($R_{\rm sp}$) [Mpc/$h$]")
plt.xlim((1e-2, 1))
plt.xscale("log")

plt.grid()
plt.savefig("%s_std_max1.png" % sys.argv[2])

#Figure 1
pretty_fig(1)
plt.xlabel(r"$\log_{10}\ $max($R_{\rm sp}$) [Mpc/$h$]")
plt.savefig("%s_max_hist.png" % sys.argv[2])

# Figure 2
pretty_fig(2)
plt.xlabel(r"$\log_{10}\ $stdev($\Delta R_{\rm sp}$) [Mpc/$h$]")
plt.savefig("%s_std_hist.png" % sys.argv[2])

#Figure 3
pretty_fig(3)
plt.xlabel(r"$R_{\rm sp}$ [Mpc/$h$]")
plt.savefig("%s_r_hist.png" % sys.argv[2])

# Figure 4
pretty_fig(4)
plt.xlim((0, 125))
plt.ylabel(r"$R_{\rm sp}$ [Mpc/$h$]")
plt.ylim((0, 1))
plt.legend(fontsize=12)
plt.savefig("%s_r_window2.png" % sys.argv[2])
