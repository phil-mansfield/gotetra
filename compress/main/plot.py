from __future__ import division

import collections
import sys

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import matplotlib.cm as cm
import matplotlib.colors as colors

s=40
L = 62.5

class MultiStd(object):
    def __init__(self, x, y):
        self._x = np.array(x)
        self._y = np.array(y)
        self._sum = np.cumsum(self._y)
        self._sqr = np.cumsum(self._y*self._y)

    def std(self, i, w):
        sum = self._sum[i + w//2 + 1] - self._sum[i - w//2]
        sqr = self._sqr[i + w//2 + 1] - self._sqr[i - w//2]
        return np.sqrt(sqr / w - (sum / w)**2)

    def std_line(self, w):
        idx = np.arange(w//2, len(self._x) - w//2 - 1)
        stds = [self.std(i, w) for i in idx]
        xs = self._x[idx]
        return xs, np.array(stds)

class MultiMin(object):
    def __init__(self, y):
        self._y = np.array(y)

    def min_line(self, w):
        dq = collections.deque(maxlen=w)
        m = np.zeros(len(self._y))
        w2 = w//2
        for i, y in enumerate(self._y):
            if len(dq) != 0:
                hd_idx, _ = dq[0]
                if len(dq) != 0 and hd_idx + w <= i:
                    dq.popleft()

            while True:
                if len(dq) == 0: break
                _, tl_val = dq[-1]
                if tl_val > y:
                    dq.pop()
                else: break

            dq.append((i, y))
            if i < w2: continue
            m[i - w2] = dq[0][1]
        for i in xrange(w2): m[i] = m[w2]
        return m

class MultiMax(object):
    def __init__(self, y):
        self._y = np.array(y)

    def max_line(self, w):
        dq = collections.deque(maxlen=w)
        m = np.zeros(len(self._y))
        w2 = w//2
        for i, y in enumerate(self._y):
            if len(dq) != 0:
                hd_idx, _ = dq[0]
                if len(dq) != 0 and hd_idx + w <= i:
                    dq.popleft()

            while True:
                if len(dq) == 0: break
                _, tl_val = dq[-1]
                if tl_val < y:
                    dq.pop()
                else: break

            dq.append((i, y))
            if i < w2: continue
            m[i - w2] = dq[0][1]
        for i in xrange(w2): m[i] = m[w2]
        return m

class MultiRange(object):
    def __init__(self, x, y):
        self._x = np.array(x)
        self._y = np.array(y)
        self._multi_min = MultiMin(self._y)
        self._multi_max = MultiMax(self._y)

    def range_line(self, w):
        mins = self._multi_min.min_line(w)
        maxes = self._multi_max.max_line(w)
        r = maxes - mins
        return self._x[w//2: -w//2 - 1], r[w//2: -w//2 - 1]

def naive_window_min(w, x):
    m = np.zeros(len(x))
    for i in xrange(len(x)):
        start, end = i - w//2, i + w//2 + 1
        m[i] = np.min(x[max(0, start): end])
    return m
    
def pad(xs, w):
    left = np.zeros(w//2, dtype=bool)
    right = np.zeros(w//2 + 1, dtype=bool)
    return np.append(left, xs)

def smooth(bs, w):
    nbs = np.zeros(len(bs), dtype=bool)
    for i in range(len(bs)):
        nbs[i] = (bs[max(0, i - w//2): i + w//2 + 1]).any()
    return nbs

x0, xs, ys, zs = zip(*np.loadtxt(sys.argv[1]))
n = len(xs)
x0 = np.array(x0)
xs = np.array(xs)
ys = np.array(ys)
zs = np.array(zs)
x0, xs, ys, zs = x0[:n-2], xs[:n-2], ys[:n-2], zs[:n-2]

def cost(n):
    return int(np.ceil(np.log(np.abs(n) + 1) / np.log(2))) + 1

def quant(xs):
    return np.abs(map(int, (xs - np.mean(xs)) * 500))
    
def range_cost(qs):
    lo = np.min(qs)
    ns = np.append(qs - lo, [len(qs), lo])
    return np.sum(map(cost, ns)), ns

def dp(xs):
    qs = quant(xs)
    costs = np.zeros(len(qs) + 1, dtype=int)
    dp_ns = [None] * (len(qs) + 1)
    dp_ns[0] = np.array([], dtype=int)

    for i in xrange(1, len(costs)):
        min, min_ns = range_cost(qs[0: i])
        mini = 0
        lo = np.min(qs[0: i])
        for j in xrange(1, i):
            c, ns = range_cost(qs[j:i])
            if c + costs[j] < min:
                min = c + costs[j]
                mini = j
                min_ns = ns
                lo = np.min(qs[j:i])
        costs[i] = min
        dp_ns[i] = np.append(dp_ns[mini], min_ns)
    return costs, dp_ns

def halo_to_ints(xs):
    costs, ns = dp(xs)
    return costs[-1] + 23 + len(xs) / float(2), ns[-1]

plt.figure(0)

plt.plot(x0, xs, "k", lw=1)
plt.scatter(x0, xs, c="k", linewidths=0, s=s)
plt.xlabel("$X(t=0)$", fontsize=18)
plt.ylabel("$X(t)$", fontsize=18)
plt.xlim(0, 62.5)

plt.figure(1)
plt.plot(ys, xs, "k", lw=1)
plt.scatter(ys, xs, c="k", linewidths=0, s=s)
plt.xlabel(r"$Y(t)$", fontsize=18)
plt.ylabel("$X(t)$", fontsize=18)
plt.xlim(0, 62.5)

diff = xs[1:] - xs[:-1]
is_slope = np.zeros(len(diff) + 1, dtype=bool)
is_pos = True
start = 0

n_slope = 3
diff_slope = 0.0625

class BaseHalo(object):
    def __init__(self, x0, xs):
        self.xs, self.x0 = xs, x0
        if len(xs) == 0:
            self.nil_halo = True
            return
        else:
            self.nil_halo = False

        self.min = np.min(self.xs)
        self.max = np.max(self.xs)
        self.med = np.median(self.xs)
        self.s = np.sign(self.xs - self.med)

    def is_halo(self): return True

    def contains(self, h):
        return self.min < h.med < self.max

class Halo(BaseHalo):
    def __init__(self, range, x0, xs):
        lo, hi = range
        BaseHalo.__init__(self, x0[lo:hi], xs[lo:hi])

class BaseVoid(object):
    def __init__(self, range, x0, xs):
        lo, hi = range
        self.xs, self.x0 = xs[lo:hi], x0[lo:hi]

    def is_halo(self): return False

class Void(object):
    def __init__(self, range, x0, xs):
        lo, hi = range
        self.xs, self.x0 = xs[lo:hi], x0[lo:hi]

    def is_halo(self): return False

class LagrangianLine(object):
    def __init__(self, void_ranges, x0, xs):
        # Create halo ranges.
        if len(void_ranges) == 0:
            halo_ranges = [(0, len(xs))]
        else:
            void_start, void_end = void_ranges[0][0], void_ranges[-1][1]
            if void_ranges[0][0] != 0:
                halo_ranges = [(0, void_start)]
            else:
                halo_ranges = []

            for i in xrange(len(void_ranges) - 1):
                halo_ranges.append((void_ranges[i][1], void_ranges[i+1][0]))

            if void_end != len(xs):
                halo_ranges.append((void_end, len(xs)))
                
        # Convert ranges to objects.
        if len(void_ranges) == 0:
            self.objs = [Halo(halo_ranges[0], x0, xs)]
        elif len(halo_ranges) == 0:
            self.objs = [Void(void_ranges[0], x0, xs)]
        else:
            self.objs = []
            if halo_ranges[0][0] < void_ranges[0][0]:
                # halos are first.
                for i in xrange(len(void_ranges)):
                    self.objs.append(Halo(halo_ranges[i], x0, xs))
                    self.objs.append(Void(void_ranges[i], x0, xs))
                if len(halo_ranges) > len(void_ranges):
                    self.objs.append(Halo(halo_ranges[-1], x0, xs))
            else:
                for i in xrange(len(halo_ranges)):
                    self.objs.append(Void(void_ranges[i], x0, xs))
                    self.objs.append(Halo(halo_ranges[i], x0, xs))
                if len(void_ranges) > len(halo_ranges):
                    self.objs.append(Void(void_ranges[-1], x0, xs))

        self._merge()
   
    def _merge(self):
        if len(self.objs) <= 1: return

        while True:
            for i in xrange(1, len(self.objs)):
                if i == len(self.objs) - 1: break
                h = self.objs[i]
                if h.is_halo() and h.nil_halo:
                    self._combine_halo(i-1, i, i+1)
                    break
            if i >= len(self.objs) - 1: break

        while True:
            for i in xrange(1, len(self.objs)):
                if i == len(self.objs) - 1: break
                v = self.objs[i]
                if not v.is_halo():
                    h_l, h_r = self.objs[i-1], self.objs[i+1]
                    if h_l.contains(h_r) or h_r.contains(h_l):
                        self._combine_halo(i-1, i, i+1)
                        break
            if i >= len(self.objs) - 1: break
                

    def _combine_halo(self, *idxs):
        base = np.array([])
        
        f_x0 = flatten([self.objs[i].x0 for i in idxs])
        f_xs = flatten([self.objs[i].xs for i in idxs])
        h = BaseHalo(f_x0, f_xs)

        idxs = sorted(idxs)
        self.objs[idxs[-1]] = h
        for i in (idxs[:-1])[::-1]: self.objs.pop(i)
    
    def _combine_void(self, *idxs):
        base = np.array([])
        
        f_x0 = flatten([self.objs[i].x0 for i in idxs])
        f_xs = flatten([self.objs[i].xs for i in idxs])
        v = BaseVoid(f_x0, f_xs)
        
        idxs = sorted(idxs)
        self.objs[idxs[-1]] = v
        for i in (idxs[:-1])[::-1]: self.objs.pop(i)

    def plot(self):
        for i in xrange(len(self.objs)):
            if self.objs[i].is_halo():
                self._plot_halo(i)
            else:
                self._plot_void(i)

    def _plot_halo(self, i):
        h = self.objs[i]
        if h.nil_halo: return
        plt.scatter(h.x0, h.xs, linewidth=0, c="red")
        plt.plot([h.x0[0], h.x0[-1]], [h.med, h.med], "red")

    def _plot_void(self, i):
        plt.scatter(self.objs[i].x0, self.objs[i].xs, linewidth=0, c="blue")

    def cost(self):
        cost = 0.0
        ns = np.array([], dtype=int)
        for obj in self.objs:
            if obj.is_halo():
                c, n = halo_to_ints(obj.xs)
                cost += c
                ns = np.append(ns, n)
        return cost, ns

    def halo_ps(self):
        n = 0
        for obj in self.objs:
            if obj.is_halo(): n += len(obj.xs)
        return n

    def void_ps(self):
        n = 0
        for obj in self.objs:
            if not obj.is_halo(): n += len(obj.xs)
        return n

def flatten(arrays):
    base = arrays[0]
    for i in xrange(1, len(arrays)):
        base = np.append(base, arrays[i])
    return base

ranges = []
for i in xrange(len(diff)):
    if is_pos:
        if diff[i] < diff_slope:
            is_pos = False
            if i - start >= n_slope:
                # The +1 usage is not a mistake.
                ranges.append((start+1, i))
                start = i
    else:
        if diff[i] > diff_slope:
            is_pos = True
            start = i

line = LagrangianLine(ranges, x0, xs)

plt.figure(0)
line.plot()

c, ns = line.cost()
h_ps = line.halo_ps()
v_ps = line.void_ps()

print "Bits per halo float: %g" % (c / h_ps)
print "Min bits per float: %g" % (c / (h_ps + v_ps))
print "Max bits per float: %g" % ((c  + (23 * v_ps)) / (h_ps + v_ps))

plt.figure()
plt.hist(ns, bins=100)
plt.figure()

plt.figure()
plt.hist(map(cost, ns))

plt.show()
