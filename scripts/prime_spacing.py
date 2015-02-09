import numpy as np
import matplotlib.pyplot as plt

PRIME_COUNT = 1000
WIDTH = 1<<10

vol = WIDTH * WIDTH * WIDTH

ps = [2, 3]

def append_prime(ps):
    p = ps[-1]

    while True:
        p+=2
        for pp in ps:
            if pp*pp > p:
                ps.append(p)
                return
            if p % pp == 0:
                break

for i in xrange(PRIME_COUNT):
    append_prime(ps)

axis = np.arange(0, WIDTH, dtype=int)
xs, ys, zs = np.meshgrid(axis, axis, axis)

xs = np.reshape(xs, vol)
ys = np.reshape(ys, vol)
zs = np.reshape(zs, vol)

dists = []

for p in ps[100:]:
    print p
    p_idxs = np.arange(0, vol, p, dtype=int)[1:]
    p_xs, p_ys, p_zs = xs[p_idxs], ys[p_idxs], zs[p_idxs]

    x_leg = np.minimum(p_xs**2, (WIDTH - p_xs)**2)
    y_leg = np.minimum(p_ys**2, (WIDTH - p_ys)**2)
    z_leg = np.minimum(p_zs**2, (WIDTH - p_zs)**2)

    dist = min(x_leg**2 + y_leg**2 + z_leg**2)
    dists.append(dist)

print "MEOW?"

dists = np.sqrt(dists)
plt.title("%d" % WIDTH)
plt.plot(ps[100:], dists)
plt.show()
