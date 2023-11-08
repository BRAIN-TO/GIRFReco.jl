## Recon loop manipulator dimension selector of what has to be reconstructed
# Probably the spiral_recon-Wrapper should be a function to take selector as Argument
# but then running as a script for debugging might be less convenient
@info "Starting in $(pwd())"

global selector = Dict{Symbol,Any}()
global is_called_from_global_recon::Bool = true
for avg = 1:4 #1:4
    for dif = 0:30 #1:10 # 0:30 # 0:6 # 0 is b=0 image
        selector[:avg] = avg
        selector[:seg] = 1
        selector[:dif] = dif
        @info "\n\n\n\tReconstructing avg=$(selector[:avg]), seg=$(selector[:seg]), dif=$(selector[:dif])\n\n"
        include("spiral_recon_single_interleave.jl")
    end
end
is_called_from_global_recon = false;