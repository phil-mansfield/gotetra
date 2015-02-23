## Output format

Output files will each consist of a binary header followed by a `float32` array (now that we aren't constrained to use a numpy-readable format we can save ourselves some disk space and `scp` time). I'm going to err on the side of including too much in the header instead of too little for now. The header will have the following format:

    type GridHeader struct {
        Type TypeInfo
        Cosmo CosmoInfo
        Render RenderInfo
        Loc LocationInfo
    }

    type TypeInfo struct {
        Endianness int64
        HeaderSize int64
        GridType int64
        IsVectorGrid int64
    }

    type CosmoInfo struct {
        Redshift, ScaleFactor float64
        OmegaM, OmegaL, Hubble float64
        RhoMean, RhoCritical float64
        BoxWidth float64
    }

    type RenderInfo struct {
        Particles int64
        TotalPixels int64
        SubsampleLength int64
        MinProjectionDepth int64
    }

    type LocationInfo struct {
        Origin, Span Vector
        PixelOrigin, PixelSpan IntVector
        PixelWidth float64
    }

    type Vector [3]float64
    type IntVector [3]int64
    
####Field descriptions
    
* `Endianness` - A flag indicating the endianness of the data in this file. -1 indicates "little endian" byte ordering and 0 indicates "big endian" byte ordering. (I choose these two values because they're the only two int values which are the same for both endiannesses). This way files can written and read on two different architectures without ambiguity.
* `HeaderSize` - The size of the header in bytes. Can be checked to ensure consistency.
* `Cosmo` - Basic information about the cosmology of the simulation that generated this grid. Contains some redundant information for ease of use.
* `Render` - Basic information about the rendering configuration variables.
* `GridType` - An integer flag indicating the type of information stored in the grid. Can be set in the range 1 - 4, corresponding to density, density gradient, divergence, and curl, respectively. 
* `IsVec` - An integer set to 0 if one grid is stored in the file and 1 if there are three (x, y, and z). Redundant with `GridType`.
* `Origin` - Vector pointing to the lowermost x, y, and z coordinates of the grid.
* `Span` - The x-width, y-width, and z-width of the grid.
* `PixelOrigin` - The (x, y, z) coordinates of the lowermost pixel in the grid.
* `PixelSpan` - The x-width, y-width, and z-width of the grid in pixels.
* `PixelWidth` - The width of a single pixel in Mpc/h.
* `Len` - The length of the grid. Redundant with `PixelSpan`, but easier to use.
   
As this header is fairly large, I will provide implementations for reading it into popular languages (read: Python and C). The interfaces will look like the following (note that I am changing my variable naming scheme to be consistent with language standards. Bye bye PascalCase ;____; )

#### Python:

See `gotetra.py` in the current directory. Note that for ease of use, some long but commonly used fields have been reproduced with shorter names at the top namespace of the class.

#### C:

    typedef struct { ... } gotetra_header_t
    void gotetra_read_header(char *filename, gotetra_header_t *hd)
    void gotetra_read_grid(char *filename, float *grid)
    void gotetra_read_vec(char *filename, float *x_grid, float *y_grid, float *z_grid)
    int gotetra_grid_length(char *filename) // For convenience
    static const int DENSITY, GRADIENT, DIVERGENCE, CURL 