from __future__ import division

import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import scipy.stats as stats
import scipy.interpolate as intr
import deriv
import profile

# figure numbers
fig_profile = 0
fig_valid_profile = 1
fig_mean_profile = 2
fig_mean_edge = 3

cs = ["Red", "Blue",
      "Chartreuse", "Cyan",
      "DarkOrange", "Magenta",
      "Black", "BlueViolet",
      "DarkGoldenRod", "DarkRed"]

def color(i): return cs[i] if i < len(cs) else "w"

def pretty_fig(n):
    """ pretty_fig(n) is equivalent to plt.figure(n), except it also
    sets a number of options which make the resulting plot look nicer.
    """
    plt.figure(n, figsize=(8, 8))
    plt.rc('text', usetex=True)
    plt.rc('font',size=19)
    plt.rc('xtick.major',pad=5); plt.rc('xtick.minor',pad=5)
    plt.rc('ytick.major',pad=5); plt.rc('ytick.minor',pad=5)

def read_profiles(fname):
    rows = np.loadtxt(sys.argv[1])
    vec_xs, vec_ys, vec_zs = rows[:3]
    idxs = np.arange(len(vec_xs) // 2) * 2
    vec_xs, vec_ys, vec_zs = vec_xs[idxs], vec_ys[idxs], vec_zs[idxs]
    rows = rows[3:]
    cols = np.array(map(lambda *x: x, *rows))
    idxs = np.arange(len(cols) // 2)
    rs, rhos = cols[idxs*2], cols[idxs*2  + 1]

    n = len(rs)
    return vec_xs[:n], vec_ys[:n], vec_zs[:n], rs[:n], rhos[:n]

def write_caustic(profs, sm_profs, fname):
    lines, idx = [], 0
    for sm_prof in sm_profs:
        for i in xrange(sm_prof.segments()):
            if sm_prof.r_sp(i) is not None:
                lines.append(sm_prof.string(i, color(idx)))
        if sm_prof.contains_caustic(): idx += 1
    with open(fname, "w+") as fp:
        fp.write("\n".join(lines))
    """
    lines = []
    for sm_prof in sm_profs:
        for i in xrange(sm_prof.segments()):
            r = sm_prof.r_sp(i)
            if r is not None:
                vec = r * sm_prof.unit()
                lines.append("%g %g %g" % (vec[0], vec[1], vec[2]))
    with open(fname, "w+") as fp:
        fp.write("\n".join(lines))
    """

def plot_profiles(profs, sm_profs, out_prefix):
    pretty_fig(fig_profile)
    for i, sm_prof in enumerate(sm_profs[:len(cs)]):
        for j in xrange(sm_prof.segments()):
            plt.plot(sm_prof.rs(j), sm_prof.rhos(j), c=color(i), lw=3)
            r = sm_prof.r_sp(j)
            if r is not None: plt.plot([r, r], [1e-2, 1e5], c=color(i))
        
    plt.xlabel(r"$R$ [Mpc/$h$]")
    plt.xscale("log")
    plt.xlim((0.01, 1))

    plt.ylabel(r"$\rho / \rho_{\rm m}$")
    plt.yscale("log")
    plt.ylim((1e-2, 1e5))

    plt.savefig("%s_profiles.png" % out_prefix)
        
def plot_valid_profiles(profs, sm_profs, out_prefix):
    pretty_fig(fig_valid_profile)
    idx, i = 0, 0
    while idx < len(cs) and i < len(sm_profs):
        sm_prof = sm_profs[i]
        if not sm_prof.contains_caustic(): 
            i += 1
            continue
        for j in xrange(sm_prof.segments()):
            plt.plot(sm_prof.rs(j), sm_prof.rhos(j), c=color(idx), lw=3)
            r = sm_prof.r_sp(j)
            if r is not None: plt.plot([r, r], [1e-2, 1e5], c=color(idx))

        idx += 1
        i += 1
        

    plt.xlabel(r"$R$ [Mpc/$h$]")
    plt.xscale("log")
    plt.xlim((0.01, 1))

    plt.ylabel(r"$\rho / \rho_{\rm m}$")
    plt.yscale("log")
    plt.ylim((1e-2, 1e5))

    plt.savefig("%s_valid_profiles.png" % out_prefix)

def savgol_nonuniform(xs, ys, pts, window, wrap=True):
    if wrap:
        xs = np.append(xs - 2*np.pi, [xs, xs + 2*np.pi])
        ys = np.append(ys, [ys, ys])
    sp = intr.UnivariateSpline(xs, ys, s=0)
    new_xs = np.linspace(np.min(xs), np.max(xs), pts * 3)
    new_ys = sp(new_xs)
    smoothed_ys = signal.savgol_filter(new_ys, window, 4)
    if wrap:
        return new_xs[pts: 2*pts], smoothed_ys[pts: 2*pts]
    else:
        return new_xs, smoothed_ys 

def plot_mean_edge(profs, sm_profs, out_prefix):
    pretty_fig(fig_mean_edge)

    sg_window = sm_profs[0].sg_window()
    sm_group = profile.ProfileGroup2D(sm_profs)
    sm_rs, sm_angles = sm_group.rs(), sm_group.angles()
    sm_xs, sm_ys, sm_zs = sm_group.coords(sm_rs, sm_angles)

    def select2d(x, y, z):
        if sm_group.axis() == 0:
            return y, z
        elif sm_group.axis() == 1:
            return x, z
        else:
            return x, y

    sm_x0s, sm_x1s = select2d(sm_xs, sm_ys, sm_zs)
    x0_str, x1_str = select2d(r"X [Mpc/$h$]", r"Y [Mpc/$h$]", r"Z [Mpc/$h$]")

    plt.plot(sm_x0s, sm_x1s, ".", c="k", lw=3, label="Smoothed points")
    plt.xlabel(x0_str)
    plt.ylabel(x1_str)

    ylo, yhi = plt.ylim()
    plt.ylim((yhi, ylo))

    plt.legend(fontsize=12)
    plt.savefig("%s_mean_edge.png" % out_prefix)
    return

    for (max_window, c) in zip([21, 61, 81], ["DarkBlue", "b", "DeepSkyBlue"]):
        msm_profs = []
        for prof in profs:
            msm_profs.append(
                prof.multi_smoothed_profile(sg_window, max_window)
            )

        group = profile.ProfileGroup2D(msm_profs)
        rs, angles = group.rs(), group.angles()
        xs, ys, zs = group.coords(rs, angles)
        x0s, x1s = select2d(xs, ys, zs)
        plt.plot(
            x0s, x1s, ".", c=c, lw=3,
            label="Multi-smoothed points (%d)" % max_window,
        )
    """
    def my_min(xs):
        if len(xs) == 0: return np.nan
        return np.min(xs)
    min_rs, min_edges, counts = stats.binned_statistic(
        np.array(angles), np.array(rs), statistic=my_min,
        range=(-np.pi, np.pi), bins=70,
    )

    #######################################################
    # Check if there are more counts than average or not. #
    #######################################################

    min_angles = (min_edges[1:] + min_edges[:-1]) / 2

    def wrap(vals): return np.append(vals, vals[0])
    def rm_nan(xs, ys):
        idxs = np.logical_not(np.isnan(xs) | np.isnan(ys))
        return xs[idxs], ys[idxs]

    min_angles, min_rs = rm_nan(min_angles, min_rs)
    sg_angles, sg_rs = savgol_nonuniform(min_angles, min_rs, 100, 19)
    min_angles, min_rs = wrap(min_angles), wrap(min_rs)
    sg_angles, sg_rs = wrap(sg_angles), wrap(sg_rs)

    min_xs, min_ys, min_zs = group.coords(min_rs, min_angles)
    sg_xs, sg_ys, sg_zs = group.coords(sg_rs, sg_angles)
    """

    plt.xlabel(x0Str)
    plt.ylabel(x1Str)

    ylo, yhi = plt.ylim()
    plt.ylim((yhi, ylo))

    plt.legend(fontsize=12)
    plt.savefig("%s_mean_edge.png" % out_prefix)

if __name__ == "__main__":
    vec_xs, vec_ys, vec_zs, r_profs, rho_profs = read_profiles(sys.argv[1])
    units = np.array(zip(vec_xs, vec_ys, vec_zs))
    out_prefix = sys.argv[2]
    sg_window = int(sys.argv[3])
    caustic_file = sys.argv[4]

    profs, sm_profs = [], []
    for (i, (raw_rs, raw_rhos)) in enumerate(zip(r_profs, rho_profs)):
        c = cs[i] if i < len(cs) else "w"
        prof = profile.RawProfile(raw_rs, raw_rhos, units[i])
        sm_prof = prof.smoothed_profile(sg_window)

        profs.append(prof)
        sm_profs.append(sm_prof)
        
    write_caustic(profs, sm_profs, caustic_file)
    plot_profiles(profs, sm_profs, out_prefix)
    plot_valid_profiles(profs, sm_profs, out_prefix)

    if np.all(vec_xs == 0) or np.all(vec_ys == 0) or np.all(vec_zs == 0):
        plot_mean_edge(profs, sm_profs, out_prefix)
