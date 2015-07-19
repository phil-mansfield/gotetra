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

def naive_caustic_min(rs, rhos, derivs):
    min_deriv = np.argmin(derivs)
    if min_deriv == 0 or min_deriv == len(rs) - 1: return [], [], []
    return [rs[min_deriv]], [derivs[min_deriv]], [rhos[min_deriv]]

def local_caustic_min(rs, rhos, derivs):
    idxs = signal.argrelextrema(derivs, np.less)[0]
    idxs = np.array([idx for idx in idxs if idx != 0 and idx != len(rs) - 1])
    if len(idxs) == 0: return [], [], []
    min_idx = idxs[np.argmin(derivs[idxs])]
    return [rs[min_idx]], [derivs[min_idx]], [rhos[min_idx]]

def caustic_min(rs, rhos, derivs, lim=-5):
    curr_min = rhos <= np.minimum.accumulate(rhos)
    idxs = signal.argrelextrema(derivs, np.less)[0]
    idxs = np.array([idx for idx in idxs if idx != 0 and idx != len(rs) - 1])
    if len(idxs) == 0: return [], [], []
    idxs = idxs[curr_min[idxs]]
    if len(idxs) == 0: return [], [], []
    idxs = idxs[derivs[idxs] < lim]
    if len(idxs) == 0: return [], [], []
    min_idx = idxs[np.argmin(derivs[idxs])]    
    return [rs[min_idx]], [derivs[min_idx]], [rhos[min_idx]]

# rs must be evenly spaced.
def find_splashback(rs, rhos, sg_window):
    sm_rhos = signal.savgol_filter(rhos, sg_window, 4)
    dr = (rs[-1] - rs[0]) / (len(rs) - 1)
    sm_derivs = signal.savgol_filter(rhos, sg_window, 4, deriv=1, delta=dr)
    r_sps, deriv_sps, _ = caustic_min(rs, sm_rhos, sm_derivs, lim=0)
    if len(r_sps) > 0:
        return r_sps[0]
    else:
        return np.nan

def splashback_range(rs, rhos, range_max=101):
    r_sps, deriv_sps = [], []
    for sg_window in xrange(2, range_max//2):
        sg_window = sg_window * 2 + 1
        if sg_window >= len(rs): break
        r_sp, deriv_sp = find_splashback(rs, rhos, sg_window)
        if np.isnan(r_sp): break
        r_sps.append(r_sp)
    return np.array(r_sps)

if len(sys.argv) != 4:
    print "Proper usage: %s input_file plot_name window_size" % sys.argv[0]
    exit(1)

rows = np.loadtxt(sys.argv[1])
vec_xs, vec_ys, vec_zs = rows[:3]
idxs = np.arange(len(vec_xs) // 2) * 2
vec_xs, vec_ys, vec_zs = vec_xs[idxs], vec_ys[idxs], vec_zs[idxs]
rows = rows[3:]
cols = map(lambda *x: x, *rows)
profile_count = len(cols) // 2

cs = ["k", "r", "b", "g", "pink", "orange",
      "brown", "y", "DarkMagenta", "LightSlateGray"]

log_low, log_high = -2, 0
bins = 200
sg_window = int(sys.argv[3])
prof_min, prof_max = 1e-2, 1e6
deriv_min, deriv_max = -20, 10
caustic_type = "low_local"
plot_count = 10

log_converter = SplinedBinConverter(log_low, log_high, bins)

min_list, caustic_mins, deriv_list, rho_list = [], [], [], []
fig_7_list = []
print "# %9s %9s %9s %9s" % ("x comp", "y comp", "z comp", "R_caustic")
for i in range(profile_count):
    c = cs[i % len(cs)]
    rs_group, rhos_group = nan_split(cols[2*i], cols[2*i + 1])

    for (rs, rhos) in zip(rs_group, rhos_group):
        rhos = rhos[rs > 0]
        rs = rs[rs > 0]
        if len(rs) <= 1: continue
        pretty_fig(0)
        if i < plot_count:
            plt.plot(rs, rhos, c=c, lw=3)
        pretty_fig(8)
        if i < plot_count:
            plt.plot(rs, rhos, c=c, lw=3)


        pretty_fig(1)
        edges, vals = log_converter.convert(np.log10(rs), np.log10(rhos))
        if edges is None: continue 
        mids = (edges[:-1] + edges[1:]) / 2

        if len(mids) <= sg_window: continue
        smoothed_rhos = 10**signal.savgol_filter(vals, sg_window, 4)
        smoothed_rs = 10**mids
        if i < plot_count:
            plt.plot(smoothed_rs, smoothed_rhos, c=c, lw=3)

        pretty_fig(2)
        smoothed_derivs = signal.savgol_filter(
            vals, sg_window, 4, delta=log_converter.dx(), deriv=1,
        )

        if i < plot_count:
            plt.plot(smoothed_rs, smoothed_derivs, c=c, lw=3)

        if i < plot_count:
            if len(rs) <= sg_window: continue
            pretty_fig(9)
            lin_s_rhos = 10**signal.savgol_filter(
                np.log10(rhos), sg_window, 4,
            )
            plt.plot(rs, lin_s_rhos, c=c, lw=3)

            pretty_fig(10)
            dr = (rs[-1] - rs[0]) / (len(rs) - 1)
            lin_s_derivs = signal.savgol_filter(
                np.log10(rhos), sg_window, 4, delta=dr, deriv=1,
            )
            plt.plot(rs, lin_s_derivs, c=c, lw=3)

        if caustic_type == "naive":
            caustic_mins, derivs, rhos = naive_caustic_min(
                smoothed_rs, smoothed_rhos, smoothed_derivs
            )
        elif caustic_type == "local":
            caustic_mins, derivs, rhos = local_caustic_min(
                smoothed_rs, smoothed_rhos, smoothed_derivs
            )
        elif caustic_type == "low_local":
            caustic_mins, derivs, rhos = low_local_caustic_min(
                smoothed_rs, smoothed_rhos, smoothed_derivs
            )

        if i < plot_count:
            for cm in caustic_mins:
                pretty_fig(0)
                plt.plot([cm, cm], [prof_min, prof_max], "-", c=c)
                pretty_fig(1)
                plt.plot([cm, cm], [prof_min, prof_max], "-", c=c)
                pretty_fig(2)
                plt.plot([cm, cm], [deriv_min, deriv_max], "-", c=c)
                fig_7_list.append((cm, rhos[0], c))
                pretty_fig(8)
                plt.plot([cm, cm], [prof_min, prof_max], "-", c=c)
                pretty_fig(9)
                plt.plot([cm, cm], [prof_min, prof_max], "-", c=c)
                pretty_fig(10)
                plt.plot([cm, cm], [-50, +50], "-", c=c)
                
        min_list += caustic_mins
        deriv_list += derivs
        rho_list += rhos

        if i < plot_count:
            lr_range, deriv_range = splashback_range(mids, vals)
            r_range = 10**lr_range
            windows = np.arange(2, len(r_range) + 2) * 2 + 1
            # I've stopped even caring at this point.
            pretty_fig(11)
            plt.plot(windows, r_range, c=c, lw=3)
            m_dr = np.max(np.abs(r_range[1:] - r_range[:-1]))
            mean, stddev = np.mean(r_range), np.std(r_range)
            plt.plot(windows, r_range, "o", c=c, 
                     label="%.3f %.3f" % (m_dr/mean, stddev/mean))
            
            pretty_fig(12)
            plt.plot(1 / windows, deriv_range, c=c, lw=3)
            plt.plot(1 / windows, deriv_range, "o", c=c)

    r_c = np.nan if len(caustic_mins) == 0 else caustic_mins[-1]
    print ("  %9.4g %9.4g %9.4g %9.4g" % 
           (vec_xs[i], vec_ys[i], vec_zs[i], r_c))

def fig_str(n):
    # In retrospect, an array would've been better...
    if n == 0:
        return "prof"
    elif n == 1:
        return "smooth"
    elif n == 2:
        return "deriv"
    elif n == 3:
        return "hist"
    elif n == 4:
        return "depth"
    elif n == 5:
        return "lin_hist"
    elif n == 6:
        return "rho"
    elif n == 7:
        return "r_rho"
    elif n == 8:
        return "lin_prof"
    elif n == 9: 
        return "lin_smooth"
    elif n == 10:
        return "lin_deriv"
    elif n == 11:
        return "r_window"
    elif n == 12:
        return "deriv_window"

if profile_count > 1:        
    bin_num = 40

    pretty_fig(3)
    plt.hist(np.log10(min_list), range=(-2, 0), bins=bin_num)

    pretty_fig(4)
    deriv_means, edges, _ = stats.binned_statistic(
        np.log10(min_list), deriv_list,
        statistic="mean", bins=bin_num, range=(-2, 0),
    )
    deriv_sqrs, edges, _ = stats.binned_statistic(
        np.log10(min_list), np.array(deriv_list)**2,
        statistic="mean", bins=bin_num, range=(-2, 0),
    )
    deriv_stds = np.sqrt(deriv_sqrs - deriv_means**2)

    counts, _ = np.histogram(np.log10(min_list), range=(-2, 0), bins=bin_num)
    deriv_errs = deriv_means / np.sqrt(counts)
    rs = (edges[1:] + edges[:-1]) / 2
    plt.errorbar(10**rs, deriv_means, fmt="o", c="k", yerr=deriv_stds)
    plt.errorbar(10**rs, deriv_means, fmt= "o", c="b", yerr=deriv_errs)

    pretty_fig(5)
    plt.hist(min_list, range = (0, 1), bins=bin_num)

    pretty_fig(6)
    plt.hist(np.log10(rho_list), bins=25, range=(-2, 2))

    pretty_fig(7)
    plt.plot(min_list, rho_list, "wo")
    for r, rho, c in fig_7_list:
        plt.plot(r, rho, "o", c=c)
    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel(r"$\rho$")
    plt.xlabel(r"$R_{\rm sp}$ [Mpc/$h$]")
    plt.xlim(10**log_low, 10**log_high)
    plt.savefig("%s_%s.png" % (sys.argv[2], fig_str(7)))

def save_profile(n, log_r):
    plt.figure(n)
    plt.xlim(10**log_low, 10**log_high)
    plt.ylim(prof_min, prof_max)
    plt.yscale("log")
    if log_r: plt.xscale("log")
    plt.ylabel(r"$\rho$ [$\rho_{\rm m}$]")
    plt.xlabel(r"$r$ [Mpc/$h$]")
    plt.savefig("%s_%s.png" % (sys.argv[2], fig_str(n)))

def save_deriv(n, log_r):
    plt.figure(n)
    plt.xlim(10**log_low, 10**log_high)
    if log_r:
        plt.ylim(deriv_min, deriv_max)
        plt.xscale("log")
        plt.ylabel(r"$d \ln\rho / d \ln r$")
    else:
        plt.ylabel(r"$d \ln\rho / d r$")
        plt.ylim(-50, +50)
    plt.xlabel(r"$r$ [Mpc/$h$]")
    plt.savefig("%s_%s.png" % (sys.argv[2], fig_str(n)))

def save_hist(n, is_r):
    plt.figure(n)
    plt.ylabel(r"$N$")
    if is_r:
        plt.xlabel(r"$\log_{10}\ R_{\rm sp}$ [Mpc/$h$]")
        plt.xlim(log_low, log_high)
    else:
        plt.xlabel(r"$\log_{10}\ \rho_{\rm sm}$ [$\rho_{\rm m}$]")
        plt.xlim(-2, 2)
    plt.savefig("%s_%s.png" % (sys.argv[2], fig_str(n)))

def save_window(n, is_r):
    plt.figure(n)
    if is_r:
        plt.xlabel(r"Window size")
        plt.ylabel(r"$r$ [Mpc/$h$]")
        plt.ylim((0, 1))
        plt.legend(loc="upper right", fontsize=12)
    else:
        plt.xlabel(r"Window size")
        plt.ylabel(r"$d \ln \rho / d \ln r$")
    plt.savefig("%s_%s.png" % (sys.argv[2], fig_str(n)))

def save_lin_hist(n):
    plt.figure(n)
    plt.ylabel(r"$N$")
    plt.xlabel(r"$R_{\rm sp}$ [Mpc/$h$]")
    plt.xlim(10**log_low, 10**log_high)
    plt.savefig("%s_%s.png" % (sys.argv[2], fig_str(n)))

# lol
save_profile(0, True)
save_profile(1, True)
save_profile(8, False)
save_profile(9, False)
save_deriv(2, True)
save_deriv(10, False)
if profile_count > 1:
    save_hist(3, True)
    save_deriv(4, True)
    save_lin_hist(5)
    save_hist(6, False)

save_window(11, True)
save_window(12, False)
