using GIRFReco
using MRISimulation
using MRIReco
using ImageUtils
using Test
using Scratch

global const tmpdir = @get_scratch!("temp")
@info "For the output of tests, please refer to the directory $tmpdir"

include("testutils.jl")
include("testintegration.jl")