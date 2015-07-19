from __future__ import division

import sys
import numpy as np
import scipy.signal as signal
import scipy.stats as stats
import scipy.interpolate as intr
import deriv
import abc

class BinConverter(object):
    """ `BinConverter` converts data binned according to one scheme into
    a different binning scheme.

    Conversion is handled through integrating splines as opposed to simple
    rebinning specifically so data can be rebinned from linear ot logarithmic
    scales without losing significant amounts of information.
    """
    def __init__(self, low, high, bins):
        """ `__init__` creates a new `BinConverter` instance corresponding
        which will correspond to `bins` bins with an upper bound at `high`
        and a lower bound at `low`.
        """
        self.low, self.high, self.bins = low, high, bins
        dx = (self.high - self.low) / bins
        self.edges = np.linspace(low, high, bins + 1)
    
    def dx(self):
        """ `dx` returns the width of an individual output bin. """
        return (self.high - self.low) / self.bins

    def convert(self, xs, ys):
        """ `convert` takes an input curve defined by `xs` and `ys` and
        converts it to the binning scheme used by this `BinConverter`
        instance.
        """
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

        if len(first) == 0: edges = edges[1:]
        if len(last) == 0: edges = edges[:-1]

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
    """ append appends two numpy arrays together np.append is an abominaiton.
    """
    out, idx = np.zeros(sum(map(len, arrays))), 0
    for array in arrays:
        for i in xrange(len(array)):
            out[idx] = array[i]
            idx += 1
    return out

def nan_split(xs, ys):
    """ `nan_split` splits the input data into lists which correspond to
    contiguous non-NaN sequences in `ys` and `xs`.
    """
    xs_group, ys_group = [], []
    start_i = 0
    prev_nan = False

    for i in xrange(len(xs)):
        if np.isnan(ys[i]) or np.isnan(xs[i]):
            if not prev_nan:
                if i != start_i:
                    xs_group.append(np.array(xs[start_i: i]))
                    ys_group.append(np.array(ys[start_i: i]))
            prev_nan = True
        else:
            if prev_nan:
                start_i = i
            prev_nan = False
    if not prev_nan:
        xs_group.append(np.array(xs[start_i:len(xs)]))
        ys_group.append(np.array(ys[start_i:len(ys)]))
    return np.array(xs_group), np.array(ys_group)

class AbstractProfile(object):
    """ `AbstractProfile` is an abstract base class representing segmented
    density profiles.
    """
    __metaclass__ = abc.ABCMeta
    lderiv_lim = -5
    @abc.abstractmethod
    def __init__(self, rs_group, rhos_group, unit, lderivs_group=None):
        """ `__init__` creates a segmented density profile pointing in the
        direction pointed by the unit vector `unit`. The number of 

        `lderivs_group` gives d ln(rho) / d ln(r) for the profile if already
        computed. If not, it will be computed through fourth order finite
        differences.
        """
        self._unit = unit
        self._rs_group, self._rhos_group = rs_group, rhos_group
        self._lrs_group, self._lrhos_group = [], []
        for rs, rhos in zip(rs_group, rhos_group):
            valid = rhos > 0
            rs, rhos = rs[valid], rhos[valid]
            self._lrs_group.append(np.log10(rs))
            self._lrhos_group.append(np.log10(rhos))
                
        if lderivs_group is None:
            size_check = lambda xs: len(xs) > 5
            self._rs_group = np.array(filter(size_check, self._rs_group))
            self._lrs_group = np.array(filter(size_check, self._lrs_group))
            self._rhos_group = np.array(filter(size_check, self._rhos_group))
            self._lrhos_group = np.array(filter(size_check, self._lrhos_group))

            self._lderivs_group = []
            for i in xrange(self.segments()):
                self._lderivs_group.append(deriv.vector_deriv(
                    self._lrs_group[i], self._lrhos_group[i], order=4,
                ))
        else:
            self._lderivs_group = lderivs_group

        self._r_sp_idxs = [None]*self.segments()
        self._r_p_idxs = [None]*self.segments()

    def segments(self):
        """ `subprofiles` returns the number of segments in this profile.
        """
        return len(self._rs_group)

    def rs(self, i):
        """ `rs` returns the r-values of the `i`th segment of this profile.
        """
        return self._rs_group[i]

    def rhos(self, i):
        """ `rhos` returns the rho-values of the `i`the segment of this profile.
        """
        return self._rhos_group[i]

    def lderivs(self, i): 
        """ `lderivs` returns the d ln(rho) / d ln(r) values of the `i`th
        segment of this profile.
        """
        return self._lderivs_group[i]

    def r_sp(self, i):
        """ `r_sp` returns the splashback radius of `i`th segment of this
        profile and `None` if there is no splashback radius for this segment.

        Results of previous calls to `r_sp` are saved, so repeated calls
        are fast.
        """
        if self._r_sp_idxs[i] is None:
            self._r_sp_idxs[i] = self.r_sp_idx(i)
        if self._r_sp_idxs[i] == -1:
            return None
        return self._rs_group[i][self._r_sp_idxs[i]]

    def r_sp_idx(self, i):
        """ `r_sp_idx` returns the index of the splashback radius of the
        `i`th segment of this profile and `-1` if there is no splashback radius 
        for this segment.

        Results of previous calls to `r_sp_idx` are _not_ saved.
        """
        rs, rhos = self._rs_group[i], self._rhos_group[i]
        lderivs = self._lderivs_group[i]

        curr_min = rhos <= np.minimum.accumulate(rhos)
        idxs = signal.argrelextrema(lderivs, np.less)[0]
        idxs = np.array([idx for idx in idxs if idx != 0 and idx != len(rs)-1])
        if len(idxs) == 0: return -1
        idxs = idxs[curr_min[idxs]]
        if len(idxs) == 0: return -1
        idxs = idxs[lderivs[idxs] < self.lderiv_lim]
        if len(idxs) == 0: return -1
        return idxs[np.argmin(lderivs[idxs])]

    def contains_caustic(self):
        for i in xrange(self.segments()):
            if self.r_sp(i) is not None: return True
        return False

    def plateau(self, i):
        """ `plateau` returns the starting and ending radii of the mid-profile
        plateau in the `i`th segment of this profile and `(None, None)` if
        there is no plateau for this segment.

        Results of previous calls to `r_sp` are saved, so repeated calls
        are fast.
        """
        if self._r_p_idxs[i] is None:
            self._r_p_idxs[i] = self.plateau_idxs(i)
        if self._r_p_idxs[i][0] == -1:
            return None, None
        low_idx, high_idx = self._p_idxs
        return self._rs_group[i][low_idx], self._rs_group[i][high_idx]

    def plateau_idxs(self, i):
        """ `plateau` returns the indexes of the starting and ending radii of
        the mid-profile plateau in the `i`th segment of this profile and `-1`
        if there is no plateau for this segment.

        Results of previous calls to `r_sp` are __not__ saved.
        """
        if self._r_sp_idxs[i] is None:
            self._r_sp_idxs[i] = self.r_sp_idx(i)
        sp_idx = self._r_sp_idxs[i]
        if sp_idx == -1: return -1, -1
        
        p_end_idx = sp_idx
        for i in xrange(sp_idx-1, -1, -1):
            if rhos[i + 1] > rhos[i]:
                p_end_idx = i + 1
                break
        for p_start_idx in xrange(p_end_idx + 1):
            if rhos[p_start_idx] < rhos[p_end_idx]: break

        if p_start_idx == p_end_idx: return -1, -1
        return p_start_idx, p_end_idx

    def string(self, i, c):
        """ `string` converts this profile into human-readable summary
        information neccessary for other parts of my plotting pipline.
        """
        r, x, y, z = self.r_sp(i), self._unit[0], self._unit[1], self._unit[2]
        r = r if r is not None else np.nan
        return "%12.4g %12.4g %12.4g %12.4g %15s" % (x, y, z, r, c)

    def unit(self):
        """ Unit returns this profile's unit vector."""
        return self._unit

class RawProfile(AbstractProfile):
    """ `RawProfile` is a segmented density profile which can be used to create
    smoothed density profiles.
    """
    def __init__(self, rs, rhos, unit, conv=BinConverter(-2, 1, 300)):
        if rs[0] == 0: rs, rhos = rs[1:], rhos[1:]
        raw_rs_group, raw_rhos_group = nan_split(rs, rhos)
        raw_lrs_group = map(np.log10, raw_rs_group)
        raw_lrhos_group = map(np.log10, raw_rhos_group)
        rs_group, rhos_group = [], []
        for raw_lrs, raw_lrhos in zip(raw_lrs_group, raw_lrhos_group):
            lr_edges, lrhos = conv.convert(raw_lrs, raw_lrhos)
            if lr_edges is None: continue
            rs_group.append(10**((lr_edges[1:] + lr_edges[:-1])/2))
            rhos_group.append(10**lrhos)
        AbstractProfile.__init__(self, rs_group, rhos_group, unit)

    def smoothed_profile(self, sg_window):
        """ `smoothed_profile` returns a `SmoothedProfile` instance
        corresponding to this profile smoothed by a Savizky-Golay filter with
        the specified window size.
        """
        return SmoothedProfile(
            self._lrs_group, self._lrhos_group, self._unit, sg_window,
        )

    def multi_smoothed_profile(self, sg_window, max_window=81):
        return MultiSmoothedProfile(
            self._lrs_group, self._lrhos_group, self._unit,
            sg_window, max_window=max_window,
        )

class SmoothedProfile(AbstractProfile):
    """ `SmoothedProfile` is a segmented density profile which has been smoothed
    by a Savitzky-Golay filter.
    """
    def __init__(self, lrs_group, lrhos_group, unit, sg_window):
        sm_rhos_group, sm_lderivs_group, sm_rs_group = [], [], []
        self._sg_window = sg_window
        for lrs, lrhos in zip(lrs_group, lrhos_group):
            if len(lrs) <= sg_window: continue
            sm_rhos = signal.savgol_filter(lrhos, sg_window, 4)
            dlr = (lrs[-1] - lrs[0]) / len(lrs)
            sm_lderivs = signal.savgol_filter(
                lrhos, sg_window, 4, deriv=1, delta=dlr
            )
            sm_rs_group.append(10**lrs)
            sm_rhos_group.append(10**sm_rhos)
            sm_lderivs_group.append(sm_lderivs)

        AbstractProfile.__init__(
            self, sm_rs_group, sm_rhos_group,
            unit, sm_lderivs_group,
        )

    def sg_window(self): return self._sg_window

class MultiSmoothedProfile(AbstractProfile):
    _max_dr_r = 0.1
    _max_std_r = 0.1
    _stop_early = True

    def __init__(self, lrs_group, lrhos_group, unit, sg_window, max_window=81):
        self._unit = unit
        self._sg_window = sg_window
        self._prof_idx = (sg_window - 5) // 2
        self._max_window = max_window

        self._profs_group, self._valid_profs_group = [], []
        for lrs, lrhos in zip(lrs_group, lrhos_group):
            profs, ok = self._multi_smooth(lrs, lrhos, unit)
            self._profs_group.append(profs)
            if ok: self._valid_profs_group.append(profs)

    def _multi_smooth(self, lrs, lrhos, unit):
        profs = []
        rs = []
        prev_r = None
        for window in xrange(5, self._max_window, 2):
            prof = SmoothedProfile([lrs], [lrhos], unit, window)
            if prof.segments() == 0: break
            r = prof.r_sp(0)
            if r is None or prof.segments() == 0: break
            if (self._stop_early and prev_r is not None and
                abs(r - prev_r) / r > self._max_dr_r):
                return profs, False
            prev_r = r
            profs.append(prof)
            rs.append(r)

        if len(profs) <= self._prof_idx: return profs, False

        rs = np.array(rs)
        r_mean, r_std = np.mean(rs), np.std(rs)
        dr_max = np.max(rs[1:] - rs[:-1])
        dr_r, std_r = dr_max / r_mean, r_std / r_mean
        return profs, True

    def segments(self):
        return len(self._valid_profs_group)

    def r_sp(self, i):
        return self._valid_profs_group[i][self._prof_idx].r_sp(0)

    def r_sp_idx(self, i):
        return self._valid_profs_group[i][self._prof_idx].r_sp_idx(0)

    def plateau(self, i):
        return self._valid_profs_group[i][self._prof_idx].plateau(0)

    def plateau_idxs(self, i):
        return self._valid_profs_group[i][self._prof_idx].plateau_idxs(0)

    def string(self, i, c):
        return self._valid_profs_group[i][self._prof_idx].string(i, c)

    def unit(self): return self._unit

    def contains_caustic(self, i): return True

    #########################################
    # MultiSmoothedProfile specific methods #
    #########################################
    def full_segments(self): return len(self._profs_group)

    def smoothed_full_profiles(self, i):
        profs = self._profs_group[i]
        return np.arange(len(profs)) * 2 + 5, profs

class ProfileGroup2D(object):
    def __init__(self, profs):
        self._rs, units = [], []
        for prof in profs:
            for i in xrange(prof.segments()):
                r = prof.r_sp(i)
                if r is not None:
                    self._rs.append(r)
                    units.append(prof.unit())
        comps = map(np.array, zip(*units))
        print len(profs), len(units)

        if np.all(comps[0] == 0):
            self._axis = 0
        elif np.all(comps[1] == 0):
            self._axis = 1
        elif np.all(comps[2] == 0):
            self._axis = 2
        else:
            raise ValueError("Profiles are not constrained to plane.")
        
        self._j0 = 1 if self._axis == 0 else 0
        self._j1 = 1 if self._axis == 2 else 2
        self._angles = [0] * len(units)
        for (i, unit) in enumerate(units):
            self._angles[i] = np.arctan2(unit[self._j1], unit[self._j0])

    def rs(self): return self._rs

    def angles(self): return self._angles
            
    def binned_mean(self, bins):
        mean, edges = stats.binned_statistic(
            self._angles, self._rs, statistic="mean",
            bins=bins, range=(0, 2 * np.pi),
        )
        return (edges[1:] + edges[:-1]) / 2, mean
        
    def binned_std(self, bins):
        mean, edges = stats.binned_statistic(
            self._angles, self._rs, statistic="mean",
            bins=bins, range=(0, 2 * np.pi),
        )
        sqr, _ = stats.binned_statistic(
            self._angles, self._rs, statistic="mean",
            bins=bins, range=(0, 2 * np.pi),
        )
        return (edges[1:] + edges[:-1]) / 2, np.sqrt(sqr - mean*mean)

    def spline(self):
        idxs = np.argsort(self._angles)
        s_angles = self._angles[idxs]
        angles = append(s_angles - 2*np.pi, s_angles, s_angles + 2*np.pi)
        rs = append(self._rs, self._rs, self._rs)
        return intr.UnivariateSpline(angles, rs, s=0)

    def savgol_filter(self, window, pts):
        sp = self.spline()
        angles = np.linspace(0, 2*np.pi, pts)
        raw_rs = sp(angles)
        sm_rs = signal.savgol_filter(raw_rs, window, 4)
        return angles, sm_rs

    def coords(self, rs, angles):
        x0s = np.zeros(len(rs))
        x1s, x2s = rs * np.cos(angles), rs * np.sin(angles)
        if self._axis == 0:
            return x0s, x1s, x2s
        elif self._axis == 1:
            return x1s, x0s, x2s
        else:
            return x1s, x2s, x0s

    def axis(self): return self._axis
