using Pkg

Pkg.activate(".") #Maybe something else is necessary here
Pkg.add("ZipFile")

# Unzip data
using ZipFile
unzip("./data.zip", exdir="./data")

# GIRFReco Setup of Data environment
using GIRFReco
download_example()

# Run the Demo
include("joss_demo.jl")
