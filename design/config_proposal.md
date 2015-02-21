###Config File Syntax

Full example config files are provided in `config_example.txt` and `bounds_example.txt`.

Config files are organized by possibly named sections and all variable declarations occur within these sections. For example, suppose that you wanted to make a config file specifying the cats that I owned:

    [CatOwner]
    Name = "Phil"
    Cats = 2
    
    [Cat "Ed"]
    Color = "White"
    Age = 5
    
    [Cat "Ned"]
    Color = "Black"
    Age = 9
    # Optional parameter:
    FangSize = "Very long"

Think of each section as a struct and the declarations within it as struct fields.

###Config File Sections

Since different `gotetra_cmd` modes have different config files, they will consist of exactly one section specifying what type of config file they are. This will be one of `[ConvertCatalog | Density | RenderDensity | FullDensity]` and will prevent the command being called on the wrong header type. Both `RenderDensity` and `Density` config files are used with the `Density` flag.

Bounds files can contain any number of sections, but must only be one of two types: `[Ball | Box]`. Both will render onto a rectangular grid, but will require different information.  All sections should have some sort of user-defined name "big_void", "h13151612", etc. That name will be used as the output file name.

### Config File Variables

##### ConvertCatalog

Standard parameters:
* `Input` - Path to the directory containing simulation catalogs to be converted.
* `Output` - Path to directory where gotetra catalogs will be written to.
* `Simulation` - ["GadgetX.Y.Z" | "ArtX.Y.Z" | ...]. Right now the only format that's supported is whatever Gadget version I'm currently reading i.
* `Cells` - int specifying how many cells to break the 

Optional parameters:
* `OutputName` - Format string for a file output name which accepts three integers, x, y, and z. Default is `density_%d%d%d.dat`.
* `IteratedInput` - Replaces `Input`. A format string which accepts a single integer and allows for multiple catalogs to be converted through a single job. i.e. instead of `Input = snapshots/snapdir_032/`, it would be `IteratedInput = snapshots/snapdir_%03d`.
* `IteratedOutput` - Replaces `Output` in the same way.
* `IterationStart` - First number iterated over for the formatted input and output directories. Default is 0.
* `IterationEnd` - Last number iterated over for the formatted input and output directories. Default is the first number which doesn't correspond to any files.

##### Density

Standard parameters:

* `Input` - Input directory.
* `Output` - Output directory.
* `Pixels` - The number of pixels that would be required to render one side length of the full box. 
* `Particles` - The number of particles per tetrahedron.

Optional parameters:

* `SubsampleLength` - Scale which the phase sheet is subsampled on. Default is 1.

##### RenderDensity

A variant of `Density` which is specifically tuned towards creating uniform images.

* `Input` - Input directory.
* `Output` - Output directory.
* `MinProjectionDepth` - The minimum number of cells you are intending to project over. Controls the number of particles per tetrahadron.
* `PixelLimit` - The number of pixels required to render the longest side length of each `Box` and each `Ball` specified in the Bounds files.

Optional parameters:

* `SubsampleLength` - Scale which the phase sheet is subsampled on. Default is 1.
* `AppendName` - List of parameters to append to fine name before writing, akin to the old `h15252115_XW315_YW314_ZW315_NP20000_G5000_S2.dat`system. Default is no parameters appended.

##### FullDensity

Standard parameters:

* `Input` - Input directory.
* `Output` - Output directory.
* `Pixels` - The number of pixels that would be required to render one side length of the full box. 
* `Particles` - The number of particles per tetrahedron.

Optional parameters:

* `SubsampleLength` - Scale which the phase sheet is subsampled on. Default is 1.
* `AppendName` - List of parameters to append to fine name before writing, akin to the old `h15252115_XW315_YW314_ZW315_NP20000_G5000_S2.dat`system. Default is no parameters appended.

##### Box

Standard Parameters:

* `X`, `Y`, `Z` - coordinates of the lowermost corner of the bounding box.
* `XWidth`, `YWidth`, `ZWidth` - width of the bounding box

Optional Parameters:

* `AlignToGrid` - Move the `X`, `Y`, and `Z` to the nearest grid vertex prior to all other calculations. Default is true for `RenderDensity` and false for `Density`

##### Ball

Renders a bounding box which completely encloses a sphere. Useful for rendering halos.

* `X`, `Y`, `Z` - Coordinates of ball center.
* `Radius` - Radius of the ball.

Optional Parameters:

* `AlignToGrid` - Move the `X`, `Y`, and `Z` to the nearest grid vertex prior to all other calculations. Default is true for `RenderDensity` and false for `Density`
* `RadiusMultiplier` - Factor to multiply the ball radius by prior to rendering. Redundant with `Radius`, but convenient. Default value is 1.