import os.path as path
import os

Ls = [2000, 500, 250, 125, 63]
fs = ["%s_id.dat", "%s_mean.dat", "%s_med.dat", "%s_plot_l4.dat"]

for i in xrange(len(fs)):
    Lfs = [path.join("L%d" % L, fs[i] % ("L%d" % L)) for L in Ls]
    out = fs[i] % "multi"
    os.system("cat %s > %s" % (" ".join(Lfs), out))
