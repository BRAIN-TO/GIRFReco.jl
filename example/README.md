# Quick Guide for Running GIRFReco Example

You may find code details of the example script [here](https://brain-to.github.io/GIRFReco.jl).

## System Configuration

Though other lower versions may still be compatible, our testing environment is Julia with version 1.9.3, and we recommend to use the same or higher version of Julia.

We also recommend to leave at least 30 GB disk space for the data of demonstration and all dependent Julia packages.

The reconstructed Cartesian and Spiral images in NIfTI format will be in the folder `joss_data_zenodo/results`. We recommend to use [FSLeyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes) as the viewer to open them.

## Steps of Running the GIRFReco Example

### (1) Cloning the package

Download or clone the repo to your local by:

```
git clone git@github.com:BRAIN-TO/GIRFReco.jl.git
```

Then enter the `example` folder:

```
cd GIRFReco.jl/example
```

Alternatively, for the one who are using Visual Studio Code, simply open the `example` sub-folder from the Menu option `File -> Open Folder...`.


### (2) Run in Julia REPL

We recommend to use Visual Studio Code with Julia extension (steps of installation can be found [here](https://code.visualstudio.com/docs/languages/julia)) to avoid possible image displaying issue, especially for those using SSH and X11 forwarding.

The Julia REPL can be launched from the Command Palette (For Windows and Linux: Press `Shift + Ctrl + P`; For Mac: Press `Shift + Command + P`) by searching the command `Julia: Start REPL`.

In the launched REPL, simpley execute the following command after launching REPL to run the whole demonstration script (including data download):

```julia
include("run_example.jl")
```

If you only want to download the demonstration dataset, run:

```julia
include("download_data.jl")
```