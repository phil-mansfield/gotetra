from __future__ import division, print_function

""" gotetra.py provides functions for interacting with gotetra output files.
When run as a script it will print the header information of the target
gotetra output file to stdout.

Supported functions:

read_header(filename) -> gotetra.Header
read_grid(filename) -> numpy.ndarray

Supported classes:

Header
  TypeInfo
  CosmoInfo
  RenderInfo
  LocationInfo
  VelocityInfo
"""

import array
import struct
import sys

import numpy as np

DENSITY             = 0
DENSITY_GRADIENT    = 1
VELOCITY            = 2
VELOCITY_DIVERGENCE = 3
VELOCITY_CURL       = 4

class Sizes(object):
    def __init__(self, ver):
        if ver == 1:
            self.header = 232
            self.type_info = 24
            self.cosmo_info = 64
            self.render_info = 40
            self.location_info = 104
        elif ver == 2:
            self.header = 336
            self.type_info = 24
            self.cosmo_info = 64
            self.render_info = 40
            self.location_info = 104
            self.velocity_info = 104
        else:
            print("Unrecognized gotetra output version, %d." % ver)
            exit(1)

def _read_endianness_version(s):
    flag = struct.unpack("q", s)
    ver = flag[0] & 0xffffffff
    end = -1 if ver >> 31 != 0 else 0
    ver = 1 + (~ver & 0x00000000ffffffff) if end else ver
    return end, ver

def _read_endianness_flag(s):
    return _read_endianness_version(s)[0]
def _read_version(s):
    return _read_endianness_version(s)[1]

def read_header(filename):
    """ read_header returns the header information at the top of a gotetra
    output file as a Header object.
    """
    with open(filename, "r") as fp:
        flag_s = fp.read(8)
        ver = _read_version(flag_s)
        sizes = Sizes(ver)
        s = fp.read(sizes.header)
    return Header(flag_s + s, sizes)

def read_grid(filename):
    """ read_grid returns the grid data stored in a gotetra output file as a
    3D numpy array if the file represents a scalar field or as a 3-tuple of
    3D numpy arrays if the file represents a vector field.

    The numpy arrays use a C-like element order, meaning the last index
    corresponds to the x-coordinate.
    """

    hd = read_header(filename)
    
    def maybe_swap(xs):
        endianness = sys.byteorder
        if endianness == "little" and hd.type.endianness_flag == -1:
            return
        elif endianness == "big" and hd.type.endianness_flag == 0:
            return
        xs.byteswap()

    n = 1
    for i in range(3):
        if i != hd.axis: n *= hd.dim[i]

    if hd.axis == 0: j, k = 1, 2
    if hd.axis == 1: j, k = 0, 2
    if hd.axis == 2: j, k = 0, 1

    if hd.type.is_vector_grid:
        xs, ys, zs = array.array("f"), array.array("f"), array.array("f")
        with open(filename, "rb") as fp:
            fp.read(hd.sizes.header + 8)
            xs.fromfile(fp, n)
            ys.fromfile(fp, n)
            zs.fromfile(fp, n)

        maybe_swap(xs)
        maybe_swap(ys)
        maybe_swap(zs)
        
        if hd.axis == -1:
            xs = np.reshape(xs, (hd.dim[2], hd.dim[1], hd.dim[0]))
            ys = np.reshape(ys, (hd.dim[2], hd.dim[1], hd.dim[0]))
            zs = np.reshape(zs, (hd.dim[2], hd.dim[1], hd.dim[0]))
        else:
            xs = np.reshape(xs, (hd.dim[k], hd.dim[j]))
            ys = np.reshape(ys, (hd.dim[k], hd.dim[j]))
            zs = np.reshape(zs, (hd.dim[k], hd.dim[j]))

        return np.array([xs, ys, zs])
    else:
        xs = array.array("f")
        with open(filename, "rb") as fp:
            fp.read(hd.sizes.header + 8)
            xs.fromfile(fp, n)

        maybe_swap(xs)
        if hd.axis == -1:
            xs = np.reshape(xs, (hd.dim[2], hd.dim[1], hd.dim[0]))
        else:
            xs = np.reshape(xs, (hd.dim[k], hd.dim[j]))
        return xs

class Header(object):
    """ Header contains header information from a gotetra header file. It
    contains the fields:
        type   : TypeInfo
        cosmo  : CosmoInfo
        render : RenderInfo
        loc    : LocationInfo
        vel    : VelocityInfo - Only included in phase-space grids.

    These contain many fields, but the three most useful are reproduced as
    top-level fields:
        dim  : np.array - The dimensions of the grid in pixels. Equivalent to
                          Header.loc.pixel_span.
        pw   : float    - The width of a single pixel. Equivalent to
                          Header.loc.pixel_width.
        axis : int      - The axis over which the image is projected over. If
                          no projection was performed and the array is 3D,
                          this will be set to -1.
    
    (These are the only fields which are truly neccessary to form images. The
    others can be learned as needed.)
    """
    def __init__(self, s, sizes):
        end = _read_endianness_flag(s[0:8])
        self.version = _read_version(s[0:8])
        ver = self.version
        self.sizes = sizes

        assert len(s) == self.sizes.header + 8
        type_start = 8
        type_end = self.sizes.type_info + type_start
        cosmo_start = type_end
        cosmo_end = cosmo_start + self.sizes.cosmo_info
        render_start = cosmo_end
        render_end = render_start + self.sizes.render_info
        loc_start = render_end
        loc_end = render_end + self.sizes.location_info

        if ver >= 2:
            vel_start = loc_end
            vel_end = vel_start + self.sizes.velocity_info

        if ver == 1:
            assert(loc_end == self.sizes.header + 8)
        elif ver >= 2:
            assert(vel_end == self.sizes.header + 8)


        self.type = TypeInfo(s[type_start: type_end], end)
        if self.type.header_size != 8 + sizes.header:
            print(
                "Header size %d is different from expected, %d." % 
                (self.type.header_size, sizes.header + 8)
            )
        self.cosmo = CosmoInfo(s[cosmo_start: cosmo_end], end)
        self.render = RenderInfo(s[render_start: render_end], end)
        self.loc = LocationInfo(s[loc_start: loc_end], end)
        if ver >= 2:
            pass

        self.dim = self.loc.pixel_span
        self.pw = self.loc.pixel_width
        self.axis = self.render.projection_axis

    def __str__(self):
        return "\n".join([
            "TypeInfo:", str(self.type), "",
            "CosmoInfo:", str(self.cosmo), "",
            "RenderInfo:", str(self.render), "",
            "LocationInfo:", str(self.loc),
        ])

def little_endian(end): return end == -1
def endian_unpack(fmt, s, end):
    fmt = "<" + fmt if little_endian(end) else ">" + fmt
    return struct.unpack(fmt, s)

class TypeInfo(object):
    """ TypeInfo contains system information about the file. This is primarily
        useful for systems purposes. Its fields are:
            endianness_flag : int - Flag indicating the endianness of the input
                                    file's endianness.
            header_size     : int - Number of bytes in the input file's header.
            grid_type       : int - A flag indicating the type of the
                                    information stored in the file.
            is_vector_grid  : bool - Flag indicating whether the file is a
                                     vector field or a scalar field.
    """
    def __init__(self, s, end):
        self.endianness_flag = end

        fmt = "qqq"
        data = endian_unpack(fmt, s, self.endianness_flag)

        self.header_size = data[0]
        self.grid_type = data[1]
        self.is_vector_grid = not (self.grid_type == DENSITY or 
                                   self.grid_type == VELOCITY_DIVERGENCE)

    def __str__(self):
        return "\n".join([
            "    endianness_flag = %s" % self.endianness_str(),
            "    header_size     = %d" % self.header_size,
            "    grid_type       = %s" % self.grid_type_str(),
            "    is_vector_grid  = %r" % self.is_vector_grid,
        ])


    def endianness_str(self):
        if little_endian(self.endianness_flag):
            return "Little Endian"
        else:
            return "Big Endian"

    def grid_type_str(self):
        if self.grid_type == DENSITY:
            return "Density"
        elif self.grid_type == DENSITY_GRADIENT:
            return "Density Gradient"
        elif self.grid_type == VELOCITY:
            return "Velocity"
        elif self.grid_type == VELOCITY_DIVERGENCE:
            return "Velocity Divergence"
        elif self.grid_type == VELOCITY_CURL:
            return "Velocity Curl"


class CosmoInfo(object):
    """ CosmoInfo contains information about the data file's underlying
    cosmology, as well as other useful physical information. Its fields are:
        redshift     : float
        scale_factor : float
        omega_m      : float
        omega_l      : float
        h0           : float - units are (km / s) / Mpc
        rho_mean     : float - units are (M_sun / h) / (Mpc / h)^3
        rho_critical : float - units are (M_sun / h) / (Mpc / h)^3
        box_width    : float - units are Mpc / h
    """
    def __init__(self, s, end):
        fmt = "d" * 8
        data = endian_unpack(fmt, s, end)
        self.redshift = data[0]
        self.scale_factor = data[1]
        self.omega_m = data[2]
        self.omega_l = data[3]
        self.h0 = data[4]
        self.rho_mean = data[5]
        self.rho_critical = data[6]
        self.box_width = data[7]

    def __str__(self):
        return "\n".join([
            "    redshift     = %.4g" % self.redshift,
            "    scale_factor = %.4g" % self.scale_factor,
            "    omega_m      = %.4g" % self.omega_m,
            "    omega_l      = %.4g" % self.omega_l,
            "    h0           = %.4g" % self.h0,
            "    rho_mean     = %.4g" % self.rho_mean,
            "    rho_critical = %.4g" % self.rho_critical,
            "    box_width    = %.4g" % self.box_width,
        ])

class RenderInfo(object):
    """ RenderInfo contains information about the free parameters used to fine
    tune gotetra when rendering the density distribution. Its fields are
        particles            : int - Number of particles per tetrahedron.
        total_pixels         : int - Number of pixels required to render to a
                                     single side of the sim box. Equal to
                                     box_width / pixel_width.
        subsample_length     : int - Level of subsampling used. 1 indicates no
                                     subsampling.
        min_projection_depth : int - Suggested minimum layers to use when
                                     forming a projected image. WARNING: This
                                     feature is totally broken.
        projection_axis      : int - If the image is projected along an axis,
                                     the index of that axis. Otherwise this
                                     will be set to -1.
    """

    def __init__(self, s, end):
        fmt = "qqqqq"
        data = endian_unpack(fmt, s, end)

        self.particles = data[0]
        self.total_pixels = data[1]
        self.subsample_length = data[2]
        self.min_projection_depth = data[3]
        self.projection_axis = data[4]

    def __str__(self):
        return "\n".join([
            "    particles            = %d" % self.particles,
            "    total_pixels         = %d" % self.total_pixels,
            "    subsample_length     = %d" % self.subsample_length,
            "    min_projection_depth = %d" % self.min_projection_depth,
            "    projection_axis      = %s" % self._projection_axis_str(),
        ])

    def _projection_axis_str(self):
        if self.projection_axis == -1:
            return "None"
        elif self.projection_axis == 0:
            return "X"
        elif self.projection_axis == 1:
            return "Y"
        elif self.projection_axis == 2:
            return "Z"
        else:
            assert(0)


class LocationInfo(object):
    """ LocationInfo contains information on the dimensions and location of the
    grid. Its fields are:
        origin : float numpy.array     - Vector to the bottommost cornder of
                                         the grid. Units are Mpc / h.
        span : float numpy.array       - Physical dimensions of box. Units are
                                         Mpc / h.
        pixel_origin : int numpy.array - Vector to the bottommost corner of
                                         the grid. Units are pixels.
        pixel_span : int numpy.array   - Dimensions of box. Units are pixels.
        pixel_width : float            - The width of a single pixel. Units are
                                         Mpc / h.
    """
    def __init__(self, s, end):
        fmt = ("d" * 6) + ("q" * 6) + "d"
        data = endian_unpack(fmt, s, end)

        self.origin = np.array([data[0], data[1], data[2]])
        self.span = np.array([data[3], data[4], data[5]])
        self.pixel_origin = np.array([data[6], data[7], data[8]])
        self.pixel_span = np.array([data[9], data[10], data[11]])
        self.pixel_width = data[12]

    def __str__(self):
        return "\n".join([
            ("    origin = [%.4g, %.4g, %.4g]" % 
             (self.origin[0], self.origin[1], self.origin[2])),
            ("    span = [%.4g, %.4g, %.4g]" %
             (self.span[0], self.span[1], self.span[2])),
            ("    pixel_origin = [%d, %d, %d]" %
             (self.pixel_origin[0],self.pixel_origin[1],self.pixel_origin[2])),
            ("    pixel_span = [%d, %d, %d]" %
             (self.pixel_span[0], self.pixel_span[1], self.pixel_span[2])),
            "    pixel_width = %.4g" % self.pixel_width,
        ])

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("%s requires a target file" % (sys.argv[0]))
        exit(1)

    print(read_header(sys.argv[1]))
