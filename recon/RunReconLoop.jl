global selector = Dict{Symbol,Any}()
 
for avg = 1:4
    for dif = 1:6
        selector[:avg] = avg;
        selector[:seg] = 1;
        selector[:dif] = dif;
        @info "Reconstructing avg=$(selector[:avg]), seg=$(selector[:seg]), dif=$(selector[:dif])"
        include("SpiralRecon_Cleaned_Mar2022_Human_SingleIntlv.jl")
    end
end