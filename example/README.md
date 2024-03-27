# Quick Guide for Running GIRFReco Example

You may find code details of the example script [here](https://brain-to.github.io/GIRFReco.jl).

## System Configuration

Though other lower versions may still be compatible, our testing environment is Julia with version 1.9.3, and we recommend to use the same or higher version of Julia.

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

Now you have two options to run the example script in Julia command-line or REPL.

### (2a) Run in Julia Command-line

You can run the example script by simply executing the following command in the `example` folder:

```
julia --project run_example.jl
```

If you only want to run data downloading part, run

```
julia --project download_data.jl
```

### (2b) Run in Julia REPL

Simpley execute the following command after launching REPL:

```julia
include("joss_demo.jl")
```

If you only want to run data downloading part, run

```julia
include("download_data.jl")
```