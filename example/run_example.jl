using Pkg

Pkg.activate(".")

include("download_data.jl")

if isdir("./joss_data_zenodo") # need to download # Run the Demo
    include("joss_demo.jl")
end