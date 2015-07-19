import numpy as np
import sys

width = int(sys.argv[1])
out = sys.argv[2]
inputs = sys.argv[3:]

grid = np.zeros(width * width * width)
for fname in inputs:
    grid += np.fromfile(fname)

grid.tofile(out)
