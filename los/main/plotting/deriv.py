import math
import scipy.interpolate as intr
import numpy as np

def vector_deriv(xs, ys, order=2):
    """ Function vector_deriv takes two np.arrays representing the x and y
    values of a function and returns a numpy array of the derivative of that
    function evaluated at those x values.
    
    Keyword arguments:
    order -- indicates the order of the returned derivative in terms of the
             local x-spacing. Only orders 2 and 4 are implemented. (default=2)
    """

    if len(xs) != len(ys):
        raise ValueError("len(xs) = %d and len(ys) = %d." % (len(xs), len(ys)))

    if order == 2:
        dy_dxs = ((np.roll(ys, +1) - np.roll(ys, -1)) /
                 (np.roll(xs, +1) - np.roll(xs, -1)))
        dy_dxs[0] =(-3.0*ys[0] + 4.0*ys[1] - ys[2]) / (xs[2] - xs[0])
        dy_dxs[-1] =-(-3.0*ys[-1] + 4.0*ys[-2] - ys[-3]) / (xs[-1] - xs[-3])
    elif order == 4:
        dy_dxs = ((-np.roll(ys, +2) + 8.0*np.roll(ys, +1) -
                   8.0*np.roll(ys, -1) + np.roll(ys, -2)) /
                  (3.0*(np.roll(xs, +2) - np.roll(xs, -2))))
        dy_dxs[0] = ((-3.0*ys[4] + 16.0*ys[3] - 36.0*ys[2] + 48.0*ys[1] -
                     25.0*ys[0]) / (3.0*(xs[4] - xs[0])))
        dy_dxs[-2]= ((-3.0*ys[-1] - 10.0*ys[-2] + 18.0*ys[-3] - 6.0*ys[-4] +
                     ys[-5]) / (3.0*(xs[-5] - xs[-1])))
        dy_dxs[1]= ((-3.0*ys[0] - 10.0*ys[1] + 18.0*ys[2] - 6.0*ys[3] +
                     ys[4]) / (3.0*(xs[4] - xs[0])))
        dy_dxs[-1] = ((-3.0*ys[-5] + 16.0*ys[-4] - 36.0*ys[-3] + 48.0*ys[-2] -
                      25.0*ys[-1]) / (3.0*(xs[-5] - xs[-1])))
    else:
        raise ValueError("order %s unrecognized or unimplemented" % str(order))

    return dy_dxs

def spline_deriv(xs, ys, s):
    return intr.UnivariateSpline(xs, ys, s=s).derivative()(xs)
    

def log_to_linear(log_xs, log_ys, dlogy_dlogxs, base=None):
    """ Function log_to_linear converts a logarithmic derivative to a linear
    derivative.

    Keyword arguments:
    base -- The log base. None indicates the natural base (default = None)
    """
    if base is None:
        return dlogx_dlogys * np.exp(log_ys - log_xs)
    else:
        return dlogx_dlogys * base**(log_ys - log_xs)

def linear_to_log(xs, ys, dy_dxs, base=None):
    """ Function linear_to_log converts a linear derivative to a logarithmic
    derivative.

    Keyword arguments:
    base -- The log base. None indicates the natural base (default = None)
    """
    return dy_dxs * xs / ys

def linear_to_log_base(xs, dy_dxs, base=None):
    """ Function linear_to_log_base converts a linear-base derivative to a
    logarithmic-base derivative.

    Keyword arguments:
    base -- The log base. None indicates the natural base (default = None)
    """
    if base is None:
        return dy_dxs * xs
    else:
        return dy_dxs * (xs * math.log(base))

def log_to_linear_base(log_xs, dlogy_dlogxs, base=None):
    """ Function linear_to_log_base converts a logarithmic-base derivative to
    a logarithmic-base derivative.

    Keyword arguments:
    base -- The log base. None indicates the natural base (default = None)
    """
    if base is None:
        return dlogy_dlogxs / np.exp(log_xs)
    else:
        return dlogy_dlogxs / (base**log_xs * math.log(base))


class ChurazovInterpolant(object):
    """ Class ChurazovInterpolant implements the interpolation scheme
    described in Appendix B of:

    Comparison of approximately isothermal gravitational potentials
    of elliptical galaxies based on X-ray and optical data
    - E. Churazov et al, MNRAS, 2010

    This method is designed to be particularly adept at finding logarithmic
    derivatives.
    """

    ###################
    # Special Methods #
    ###################

    def __init__(self, xs, ys, delta_r=0.3):
        """ Function __init__ creates a new ChurazovInterpolant based on the
        given x and y points.

        Optional Parameters:
        delta_r -- Parameter controling the width of the weighting function.
                   (default = 0.3)
        """
        self._fs = np.log(xs)
        self._ss = np.log(ys)
        self._delta_r = delta_r

    ##################
    # Public Methods #
    ##################

    def curve(self, x0s):
        """ Function curve returns the value of the fitted curve evaluated
        at each of the points in x0s.
        """
        a_params, b_params = zip(*map(self._local_params, x0s))
        return np.exp(a_params * np.log(x0s) + b_params)

    def log_deriv(self, x0s):
        """ Function log_deriv returns dln(y) / dln(x) evaluated at each
        of the points in x0s.
        """
        return np.array([self._local_params(x0)[0] for x0 in x0s])

    def deriv(self, x0s):
        """ Function log_deriv returns dln(y) / dln(x) evaluated at each
        of the points in x0s.
        """
        a_params, b_params = zip(*map(self._local_params, x0s))
        return np.exp(a_params * np.log(x0s) + b_params) * a_params / x0s

    ###################
    # Private Methods #
    ###################

    def _local_params(self, x0):
        ws = self._weights(x0)
        fs = self._fs
        ss = self._ss

        fws_sum = sum(fs * ws * ss)
        w_sum = sum(ws)
        ws_sum = sum(ws * ss)
        fw_sum = sum(fs * ws)
        ffw_sum = sum(fs * fs * ws)

        a = ((fws_sum * w_sum - ws_sum * fw_sum) /
             (ffw_sum * w_sum - fw_sum**2))
        b = (ws_sum - a * fw_sum) / w_sum
        return a, b

    def _weights(self, x0):
        f0s = math.log(x0) * np.ones(len(self._fs))
        return np.exp(- (f0s - self._fs)**2 / (2 * self._delta_r**2))
