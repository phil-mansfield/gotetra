#!/usr/bin/env python
from __future__ import division

import time

import sys
import os
import subprocess
import os.path as path

r_defs = ["vir", "200c", "500c", "200m"]
str_r_defs = "[" + " | ".join(r_defs) + "]"

usage_str = (
    """Correct usage:
$ %s halo_id r_def filter_size
  halo_id     - n | 1 <= n
  r_def       - %s
  filter_size - 2 n + 1 | 0 <= n""" % (sys.argv[0], str_r_defs)
)

def check(val):
    if not val:
        print usage_str
        exit(1)
def is_even(n):
    return (n // 2) * 2 == n

check(len(sys.argv) == 4)

try:
    if sys.argv[1] == "None":
        halo_id = None
    else:
        halo_id = int(sys.argv[1])
    filter_size = int(sys.argv[3])
except:
    check(False)
r_def = sys.argv[2]

check(r_def in r_defs)
check(filter_size > 0 and not is_even(filter_size))

gtet_dir = "/home/mansfield/code/go/src/github.com/phil-mansfield/gotetra/main/small_halos"
sh_dir = gtet_dir

hlist_file = "/project/surph/diemer/Box_L0063_N1024_CBol/Rockstar/hlists/hlist_1.00000.list"
snap_file = "/project/surph/diemer/Box_L0063_N1024_CBol/Snaps/snapdir_100/snapshot_100.0"
gtet_snap_dir = "/project/surph/mansfield/data/sheet_segments/Box_L0063_N1024_G0008_CBol/snapdir_100"


opt_str = "R%s_F%d" % (r_def, filter_size)
halo_str = "%dh" % halo_id
out_dir = path.join("output", opt_str, halo_str)
plot_dir = path.join("plots", opt_str, halo_str)
x_plot_dir = path.join(plot_dir, "x_proj")
y_plot_dir = path.join(plot_dir, "y_proj")
z_plot_dir = path.join(plot_dir, "z_proj")

if not path.exists(out_dir):
    subprocess.check_call("mkdir -p %s" % out_dir, shell=True)
if not path.exists(plot_dir):
    subprocess.check_call("mkdir -p %s" % plot_dir, shell=True)
if not path.exists(x_plot_dir):
    subprocess.check_call("mkdir -p %s" % x_plot_dir, shell=True)
if not path.exists(y_plot_dir):
    subprocess.check_call("mkdir -p %s" % y_plot_dir, shell=True)
if not path.exists(z_plot_dir):
    subprocess.check_call("mkdir -p %s" % z_plot_dir, shell=True)

gtet_file = path.join(gtet_dir, "%dh.gtet" % halo_id)
sh_file = path.join(gtet_dir, "%dh_R%s_sh.txt" % (halo_id, r_def))

out_prefix = path.join(out_dir, halo_str)

lines_file = path.join(out_dir, "%s_lines.txt" % halo_str)
x_lines_file = path.join(out_dir, "%s_x_lines.txt" % halo_str)
y_lines_file = path.join(out_dir, "%s_y_lines.txt" % halo_str)
z_lines_file = path.join(out_dir, "%s_z_lines.txt" % halo_str)

caustic_file = path.join(out_dir, "%s_caustic.txt" % halo_id)
x_caustic_file = path.join(out_dir, "%s_x_caustic.txt" % halo_str)
y_caustic_file = path.join(out_dir, "%s_y_caustic.txt" % halo_str)
z_caustic_file = path.join(out_dir, "%s_z_caustic.txt" % halo_str)

plot_prefix = path.join(plot_dir, halo_str)
x_plot_prefix = path.join(x_plot_dir, halo_str + "_x")
y_plot_prefix = path.join(y_plot_dir, halo_str + "_y")
z_plot_prefix = path.join(z_plot_dir, halo_str + "_z")


subprocess.check_call("go build read_halo.go", shell=True)
t0 = time.time()
subprocess.check_call(
    "./read_halo %s %s %s %s %s %s" % (
        hlist_file, snap_file, gtet_snap_dir, gtet_dir, gtet_dir, r_def,
    ), shell=True,
)
t1 = time.time()
subprocess.check_call("go build profile.go", shell=True)
t2 = time.time()
subprocess.check_call("./profile %s %s %s" %
                      (gtet_file, sh_file, out_prefix),
                      shell=True)
t3 = time.time()
print t1 - t0, t3 - t2
print "Radial Profile:"

subprocess.check_call("python plot_filament_remover.py %s %s %d %s" %
                      (lines_file, plot_prefix, filter_size, caustic_file),
                      shell=True)

print "X Profile:"
subprocess.check_call("python plot_filament_remover.py %s %s %d %s" %
                      (x_lines_file, x_plot_prefix, filter_size,x_caustic_file),
                      shell=True)

print "Y Profile:"
subprocess.check_call("python plot_filament_remover.py %s %s %d %s" %
                      (y_lines_file, y_plot_prefix, filter_size,y_caustic_file),
                      shell=True)
print "Z Profile:"
subprocess.check_call("python plot_filament_remover.py %s %s %d %s" %
                      (z_lines_file, z_plot_prefix, filter_size,z_caustic_file),
                      shell=True)
