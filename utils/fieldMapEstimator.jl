
function ml_cost(x::Matrix{T},y::Matrix{Complex{T}},z::Matrix{Complex{T}}, β) where T

    Ψ = (sum(abs.(conj.(y) .* z) .* (1 .- cos.(angle.(z) .- angle.(y) .- x)),dims=[1,2]) + β * R(x))[1]

end

function R(x::Matrix{T}) where T

    return T.(0.5) * sum(abs2,diff(x,dims=1),dims=[1,2]) + T.(0.5) * sum(abs2,diff(x,dims=2),dims=[1,2])

end

function pcg_ml_est_fieldmap(y,z)

    d = 1
    x = angle.(conj.(y) .* z)

    while d < 10000

        gs = gradient(Flux.params(x)) do
            ml_cost(x, y, z,0.01)
        end
        x̄ = gs[x]
        x .-= 0.01 .* x̄

        if mod(d,100) == 0
            println(ml_cost(x, y, z,0.01))
        end

        d += 1

    end

    return x

end


function estimateB0Maps(imData,slices)

    b0Maps = Complex.(zeros(size(imData)[1:3]))

    for slice in slices
        b0Maps[:,:,slice] = rotl90(pcg_ml_est_fieldmap(imData[:,:,slice,1,1],imData[:,:,slice,2,1])) ./ ((7.38-4.92)/1000)
    end

    return real.(b0Maps)

end