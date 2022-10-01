using Flux

function ml_cost(x::Matrix{T},y::Matrix{Complex{T}},z::Matrix{Complex{T}}, β) where T

    Ψ = (sum(abs.(conj.(y) .* z) .* (1 .- cos.(angle.(z) .- angle.(y) .- x)),dims=[1,2]) + β * R(x))[1]

end

function R(x::Matrix{T}) where T

    return T.(0.5) * sum(abs2,diff(x,dims=1),dims=[1,2]) + T.(0.5) * sum(abs2,diff(x,dims=2),dims=[1,2])

end

function pcg_ml_est_fieldmap(y,z,β) 

    d = 1
    κ = conj.(y) .* z
    x = angle.(κ)
    m = abs.(κ)

    trust_step = 0.005 ./ (m .+ 4*β)

    while d < 100

        gs = gradient(Flux.params(x)) do
            ml_cost(x, y,z,β) # need to interpolate the y z and β to have better performance 
        end

        x̄ = gs[x]
        x .-= trust_step .* x̄ 

        if mod(d,10) == 0
            println(ml_cost(x, y, z, β))
        end

        d += 1 # TODO add early stopping criteria

    end

    return x

end


function estimateB0Maps(imData,slices, TE1,TE2,β,isrotated)

    b0Maps = Complex.(zeros(size(imData)[1:3]))
    for slice in slices
        b0Maps[:,:,slice] = pcg_ml_est_fieldmap(imData[:,:,slice,1,1],imData[:,:,slice,2,1],β) ./ ((TE2 - TE1)/1000)
    end

    b0Maps = mapslices(isrotated ? x->x : x-> rotl90(x),b0Maps,dims=(1,2)) 

    return real.(b0Maps)

end