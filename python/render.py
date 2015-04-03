import gotetra
import sys
import matplotlib.pyplot as plt
import numpy as np

IMGS = 1

fname = sys.argv[1]

if len(sys.argv) > 2:
    dim = sys.argv[2]
else:
    dim = "X"
axis = ["Z", "Y", "X"].index(dim) # (Why did I write it like this?)
axis_x = 2 - axis

hd = gotetra.read_header(fname)
if hd.axis != -1 and axis_x != hd.axis:
    print ("specified projection axis '%s' did not match gtet " + 
           "rendering axis, '%s'") % (dim, ["X", "Y", "Z"][hd.axis])
    exit(1)

if len(sys.argv) > 3:
    depth = int(sys.argv[3])
else:
    depth = hd.render.min_projection_depth

if hd.axis != -1:
    depth = hd.loc.span[axis_x] / hd.pw

if len(sys.argv) > 4:
    cmap_name = sys.argv[4]
else:
    cmap_name = "afmhot"

if len(sys.argv) > 5:
    savefile = sys.argv[5]
else:
    savefile = None

cmaps = {
    "afmhot": ("afmhot", -2, 4),
    "ocean": ("ocean", -4, 3),
    "RdBu_r": ("RdBu_r", -3, 3),
    "cubehelix": ("cubehelix", -2, 4),
    "test": ("RdGy_r", -3, 3)
}
cmap, vmin, vmax = cmaps[cmap_name]

rhos = gotetra.read_grid(fname)
rhos = np.clip(rhos, 10**vmin, 10**vmax)

if hd.axis == -1:
    cs = [map(int, np.array(rhos.shape) * float(i + 1) / (IMGS + 1))
          for i in xrange(IMGS)]
else:
    cs = [None]

for c in cs:
    # TODO: generalize this so that it doesn't require the ugly "X", "Y", "Z"
    # conditionals and copy/pasted code.

    if hd.axis == -1:
        mid = hd.loc.origin[axis_x] + c[axis_x] * hd.pw
    else:
        mid = hd.loc.origin[axis_x] + hd.loc.span[axis_x] / 2

    if dim == "Z":
        if hd.axis == -1:
            slice = rhos[c[axis] - depth / 2: c[axis] + depth / 2 + 1, :, :]
        else:
            slice = rhos
        extent = [hd.loc.origin[1], hd.loc.origin[1] + hd.loc.span[1],
                  hd.loc.origin[0], hd.loc.origin[0] + hd.loc.span[0]] 
        xlabel = "$x$ [Mpc/$h$]"
        ylabel = "$y$ [Mpc/$h$]"
        title = (r"$z$ = $%.4g$ Mpc/$h$, $\Delta z$ = $%.4g$ Mpc/$h$" %
                 (mid, hd.pw * depth))
        if hd.axis == -1:
            px, py = rhos.shape[1], rhos.shape[2]
        else:
            px, py = rhos.shape[0], rhos.shape[1]

    elif dim == "Y":
        if hd.axis == -1:
            slice = rhos[:, c[axis] - depth / 2: c[axis] + depth / 2 + 1, :]
        else:
            slice = rhos
        extent = [hd.loc.origin[2], hd.loc.origin[2] + hd.loc.span[2],
                  hd.loc.origin[0], hd.loc.origin[0] + hd.loc.span[0]] 
        xlabel = "$x$ [Mpc/$h$]"
        ylabel = "$z$ [Mpc/$h$]"
        title = (r"$y$ = $%.4g$ Mpc/$h$, $\Delta y$ = $%.4g$ Mpc/$h$" %
                 (mid, hd.pw * depth))
        if hd.axis == -1:
            px, py = rhos.shape[0], rhos.shape[2]
        else:
            px, py = rhos.shape[0], rhos.shape[1]

    else:
        if hd.axis == -1:
            slice = rhos[:, :, c[axis] - depth / 2: c[axis] + depth / 2 + 1]
        else:
            slice = rhos
        extent = [hd.loc.origin[2], hd.loc.origin[2] + hd.loc.span[2],
                  hd.loc.origin[1], hd.loc.origin[1] + hd.loc.span[1]] 
        xlabel = "$y$ [Mpc/$h$]"
        ylabel = "$z$ [Mpc/$h$]"
        title = (r"$x$ = $%.4g$ Mpc/$h$, $\Delta x$ = $%.4g$ Mpc/$h$" %
                 (mid, hd.pw * depth))
        if hd.axis == -1:
            px, py = rhos.shape[0], rhos.shape[1]
        else:
            px, py = rhos.shape[0], rhos.shape[1]

    slice = np.log10(slice)
    if hd.axis == -1:
        img = np.mean(slice, axis=axis)
    else:
        img = slice

    print px, py
    # I hate this.
    extent[0], extent[1], extent[2], extent[3] = extent[2], extent[3], extent[1], extent[0]

    if savefile is not None:
        print px, py
        fig = plt.figure(frameon=False, figsize=(px/300.0, py/300.0))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis("off")
        ax.imshow(img, vmin=vmin, vmax=vmax, cmap=cmap,
                  interpolation="nearest", extent=extent)
        fig.savefig(savefile, dpi=300)
    else:
        plt.figure()
        plt.imshow(img, vmin=vmin, vmax=vmax, cmap=cmap,
                   interpolation="nearest", extent=extent)
        plt.title(title)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.colorbar()

if savefile is None:
    plt.show()
