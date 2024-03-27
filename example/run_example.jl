using Pkg

Pkg.activate(".")
Pkg.instantiate()

include("download_data.jl")

if isdir("./joss_data_zenodo") # need to download # Run the Demo
    include("joss_demo.jl")
end