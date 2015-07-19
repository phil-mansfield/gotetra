from __future__ import division
import numpy as np
import random
import sys
import os.path as path

pixels = 400
k = 5

# Remember: 7 arguents.
shell_script="""#!/bin/sh
#SBATCH --job-name=tetra.Render_%s
#SBATCH --output=%s/Render_%s.out
#SBATCH --error=%s/Render_%s.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --time=2:00:00
#SBATCH --mem=24GB
#SBATCH --account=pi-kravtsov
module load go/1.3
go build -gcflags=-B github.com/phil-mansfield/gotetra && time ~/code/go/src/github.com/phil-mansfield/gotetra/main/main -Render %s %s"""

# Remember: 3 arguments
render_script="""[Render]
Quantity = Density
Input  = %%s
Output = %%s
ImagePixels = %d
Particles = %%d""" % pixels

# Remember: 5 arguments
ball_script="""[Ball "%s"]
X = %g
Y = %g
Z = %g
Radius = %g
"""
def particles(L, l, pixels):
    return max(500 * (L / (3 * l))**3 * (pixels / 500)**3, 20)

print particles(63, 2 * 1.1648, 400)
assert 0



def gen_scripts(name, script_dir, in_dir, out_dir, row):
    _, R_200c, x, y, z = row
    R = R_200c * 7

    shell_name = path.join(script_dir, "%s_shell.sh" % name)
    render_name = path.join(script_dir, "%s_render.txt" % name)
    ball_name = path.join(script_dir, "%s_halo.txt" % name)
    
    shell_text = shell_script % (name, script_dir, name, script_dir,
                                 name, render_name, ball_name)
    render_text = render_script % (in_dir, out_dir, particles(63, 2*R, pixels))
    ball_text = ball_script % (name, x, y, z, R)

    with open(shell_name, "w+") as fp: fp.write(shell_text)
    with open(render_name, "w+") as fp: fp.write(render_text)
    with open(ball_name, "w+") as fp: fp.write(ball_text)

if __name__ == "__main__":
    fname = sys.argv[1]
    script_dir, in_dir, out_dir = sys.argv[2], sys.argv[3], sys.argv[4]
    script_dir = path.abspath(script_dir)
    in_dir = path.abspath(in_dir)
    out_dir = path.abspath(out_dir)

    rows = np.loadtxt(fname, usecols=(1,2,3,4,5))

    # Milky Way mass range.
    lims = [1.5e12, 0.8e12]
    def is_lim(prev, curr):
        for lim in lims:
            if prev >= lim > curr: return True
        return False

    idxs = []
    for (i, row) in enumerate(rows):
        M200c, R200c, x, y, z = row
        curr = M200c
        if i == 0:
            prev = M200c
            continue
        if is_lim(prev, curr): idxs.append(i)
        prev = curr

    rows12 = rows[idxs[0]: idxs[1]+1]
    sam12 = random.sample(rows12, k)

    for row in sam12:
        name = "m_%g" % row[0]
        gen_scripts(name, script_dir, in_dir, out_dir, row)
