# User Interface Proposal

I'll be revamping the command line interface for `gotetra` to that it is easier to use and easier to maintain. To emphasize the transformation from a hacky script to a fully-formed aspect of `gotetra`, the functionality previously provided by `main.go` will now be provided by the module and binary called `gotetra_cmd`.

##Command Line Interface

`gotetra_cmd` can be run in a number of modes. Each mode is configured by one or more configuration files which are passed to it at the command line. The modes are `DefaultConfig`, `ConvertCatalog`, `Density`, and `FullDensity`. Other flags will be added when I get this code to deal with vector quantities (giv, grad, curl).

##### DefaultConfig

    $ ./gotetra_cmd -DefaultConfig [ConvertCatalog | Density | FullDensity ]

`DefaultConfig` prints a configuration file to `stdout` corresponding to the default values for the specified mode. This can of course be converted to a usable config file via `$ ./gotetra_cmd ... > MyConfig.gcfg` This means that configuration files will not need to be commited with the repository and that users will never be put in a position where they have no example files to work off of.

##### ConvertCatalog

    $ ./gotetra_cmd -ConvertCatalog Config.txt

`ConvertCatalog` converts raw simulation particle catalogs and converts them to `gotetra`-readable catalogs. The config file specifies directory locations, initial catalog format, output file size, etc.

##### Density

    $ ./gotetra_cmd -Density Config.txt Bounds1.txt Bounds2.txt ...

`Density` replaces the old `BoundedDensity` command which computes the density field on some set of sub-boxes. The first config file specifies things like file loactions and maximum projection depth. The remaining files are each configuration files specifying boxes that will be rendered.

The density grids for all the specified boxes will be written to output files with names specified within the Bounds files. An optional paramters in the main config file will allow the user to append information like particle count, 

##### FullDensity

    $ ./gotetra_cmd -Density Config.txt

`FullDensity` is a special case of `Density` which corresponds to an entire box being rendered. There are a host of special optimizations which can (and should) be made in this case, meaning that it entails its own special mode.