import sys
import gotetra as gtet
import time

def basic_usage(filename):
    """ The output of gotetra is a custom binary format signified by the
    file extention '.gtet'. This format consists of pretty big header (your
    version uses 336 bytes) followed by either one or three arrays of 32-bit
    floats.

    You can read these files by calling gtet.read_header(filename) and
    gtet.read_grid(filename).
    """

    hd = gtet.read_header(filename)
    grid = gtet.read_grid(filename)

    """  To help you keep track of the different fields in the header and to
    ingnore the ones you don't care about, it's broken up into four subfields:
    type   - System information about the file: sometimes useful.
    loc    - Dimensional information about the grid: often useful.
    render - Information about the parameters used to render the grid: rarely
             useful.
    cosmo  - Information about the underlying cosmology used in the simulation:
             often useful.
    """

    print "Each square pixel has a width of %g Mpc/h." % hd.loc.pixel_width

    """ The most important fields are given abbreviated names:
    """
    
    print "Each square pixel (still) has a width of %g Mpc/h." % hd.pw
    print "The dimensions of the grid are %s" % hd.dim
    if hd.axis == -1:
        print "The gird is not a projection."
    else:
        print "The grid is projected along the %s-axis" % ("xyz"[hd.axis])

def grid_type(filename):
    """ read_grid(filename) will return a single grid if the .gtet file
    represents a scalar quantity and a tuple of three grids if the .gtet file
    represents a vector quantity.

    If you don't know which one you're looking at, test for it like this:
    """
    hd = gtet.read_header(filename)
    if hd.type.is_vector_grid:
        xs, ys, zs = gtet.read_grid(filename)
        print "%s is a vector grid" % filename
    else:
        vals = gtet.read_grid(filename)
        print "%s is a scalar grid" % filename

    """ Additionally, if you don't know what physical quantity is represented
    by this file, you can test for it like this:
    """
    if hd.type.grid_type == gtet.DENSITY:
        # Units: rho_mean
        print "%s is a density grid" % filename
    elif hd.type.grid_type == gtet.DENSITY_GRADIENT:
        # Units: rho_mean / (Mpc/h)
        print "%s is a gradient grid" % filename
    elif hd.type.grid_type == gtet.VELOCITY:
        # Units: km/s
        print "%s is a velocity grid" % filename
    elif hd.type.grid_type == gtet.VELOCITY_DIVERGENCE:
        # Units: km/s / (Mpc/h)
        print "%s is a divergence grid" % filename
    elif hd.type.grid_type == gtet.VELOCITY_CURL:
        # Units: km/s / (Mpc/h)
        print "%s is a curl grid" % filename
    
def indexing(filename):
    """ I always get confused by who uses what indexing convention, so I'll
    be explicit about it.

    .gtet files use C-like array ordering, so the x axis is the last index
    """
    hd = gtet.read_header(filename)

    if hd.type.is_vector_grid:
        vals, _, _ = gtet.read_grid(filename)
        print "the grid at X = 200, Y = 200, Z = 20 is %g" % vals[20, 200, 200]
    else:
        vals = gtet.read_grid(filename)
        print "the grid at X = 200, Y = 200, Z = 20 is %g" % vals[20, 200, 200]

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("%s requires a target file" % (sys.argv[0]))
        exit(1)    

    filename = sys.argv[1]
    basic_usage(filename)
    grid_type(filename)
    indexing(filename)
