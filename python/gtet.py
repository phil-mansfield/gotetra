import array
import struct
import numpy as np

class File(object):
    def __init__(self, file_name):
        """ Reads a gotetra grid file located at file_name. Returns an object
        with various header fields along with a grid of positions ordered
        according to their positions in the ICs.

        Data grids:
            x - A shape (n, n, n, 3) numpy array containing the positions of
                every particle in the file in comoving Mpc/h. The 1st [2nd, 3rd]
                index corresponds to the x [y, z] dimension in the ICs, and the
                4th index to the coordinates of the position vector. Positions
                in this array are unique and do not appear in the x array of
                other files.
           xn - The same as x, except with shape (n+1, n+1, n+1, 3). This is the
                "neighbor" array and contains all the particles which are one
                index to the x [y, z] + 1 of the position vectors in x,
                including a buffer region which is shared with the x array of
                adjacent files.
        v, vn - The same as above, except velocities in physical km/s.

        Header fields:
                  z - Redshift.
            omega_m - (rho_m/rho_c)(z=0).
               h100 - H0 / 100 (km/s/Mpx).
              count - Number of paritlces in the box.
        count_width - count^(1/3).
            x_width - The width of x on one side.
           xn_width - The width of xn on one side.
           file_idx - Unique index of the file.
         file_cells - Width of the file grid on one side. file_cells*x_width =
                      count_width.
           x_origin - Position of the "lowest" corner of the bounding box around
                      x positions in the file in comoving Mpc/h. Will always be
                      within the simulation domain.
             x_span - Displacement vector from the "lowest" to "highest" corner
                      of the position bounding box in comoving Mpc/h.
           v_origin - Same as x_origin but for velocities in physical km/s.
             v_span - Same as x_origin but for velocities in physical km/s.
        """
        fp = open(file_name, "rb")

        _ = fp.read(8) 
        cosmo_info = struct.unpack("dddd", fp.read(4*8))
        index_info = struct.unpack("qqqqqqq", fp.read(7*8))
        misc_info = struct.unpack("dd", fp.read(2*8))
        position_bounds = struct.unpack("ffffff", fp.read(6*4))
        velocity_bounds = struct.unpack("ffffff", fp.read(6*4))
        
        self.z, self.scale, self.omega_m, self.h100 = cosmo_info
        (self.count, self.count_width,
         self.segment_width, self.grid_width, self.grid_count,
         self.file_idx, self.file_cells) = index_info
        _, self.box_width = misc_info

        self.x_width, self.xn_width = self.segment_width, self.grid_width
        
        self.x_origin = np.array(position_bounds[:3])
        self.x_span = np.array(position_bounds[3:])
        self.v_origin = np.array(velocity_bounds[:3])
        self.v_span = np.array(velocity_bounds[3:])

        n = int(self.grid_width**3)
        
        x = array.array("f")
        x.fromfile(fp, 3*n)
        x = np.array(x)
        
        v = array.array("f")
        v.fromfile(fp, 3*n)
        v = np.array(v)
        
        gw = self.grid_width
        self.xn = x.reshape((gw, gw, gw, 3))
        self.vn = v.reshape((gw, gw, gw, 3))
        self.x = self.xn[:-1,:-1,:-1,:]
        self.v = self.vn[:-1,:-1,:-1,:]

        fp.close()

if __name__ == "__main__":
    file_name = "/home/phil/code/src/github.com/phil-mansfield/tetra_deformation/data/L500/sheet000_snap100.dat"

    f = File(file_name)
    #print(f.xn[0,0,:,:])
    #print()
    #print(f.xn[0,:,0,:])
    #print()
    #print(f.xn[:,0,0,:])
