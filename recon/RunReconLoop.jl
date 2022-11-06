## Recon loop manipulator dimension selector of what has to be reconstructed
# Probably the SpiralRecon-Wrapper should be a function to take selector as Argument
# but then running as a script for debugging might be less convenient

global selector = Dict{Symbol,Any}()
global isCalledFromReconLoopGlobal::Bool = true
for avg = 1 #1:4
    for dif = 0:30 # 0:6 # 0 is b=0 image
        selector[:avg] = avg;
        selector[:seg] = 1;
        selector[:dif] = dif;
        @info "\n\tReconstructing avg=$(selector[:avg]), seg=$(selector[:seg]), dif=$(selector[:dif])\n\n"
        include("SpiralRecon_Cleaned_Mar2022_Human_SingleIntlv.jl")
    end
end
isCalledFromReconLoopGlobal = false;