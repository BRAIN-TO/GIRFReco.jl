using GIRFReco
using MRISimulation
using MRIReco
using MRIFiles
using ImageUtils
using Test
using Scratch
using Statistics
using JLD2

global const tmpdir = @get_scratch!("temp")
@info "For the output of tests, please refer to the directory $tmpdir"

include("testio.jl")
include("testutils.jl")
include("testintegration.jl")