import os
import sys

N = sys.argv[1]
out = sys.argv[2]
fnames = sys.argv[3:]

for name in fnames:
    print name
    os.system("python density_comp.py %s %s %s" % (name, N, out))
