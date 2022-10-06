using Flux

"""
ml_cost(x::Matrix{T},y::Matrix{Complex{T}},z::Matrix{Complex{T}}, β)
Calculates the ML estimator cost between an estimated phase map and the underlying multi-echo scan data
# Arguments
* `x::Matrix{T}` - fieldmap estimate (in radians)
* `y::Matrix{Complex{T}}` - Complex first-echo image data
* `z::Matrix{Complex{T}}` - Complex second-echo image data
* `β` - Regularization parameter controlling roughness penalty
"""
function ml_cost(x::Matrix{T},y::Matrix{Complex{T}},z::Matrix{Complex{T}}, β) where T

    Ψ = (sum(abs.(conj.(y) .* z) .* (1 .- cos.(angle.(z) .- angle.(y) .- x)),dims=[1,2]) + β * R(x))[1]

end

"""
R(x::Matrix{T})
Regularization function which penalizes roughness
# Arguments
* `x::Matrix{T}` - fieldmap estimate (in radians)
"""
function R(x::Matrix{T}) where T

    return T.(0.5) * sum(abs2,diff(x,dims=1),dims=[1,2]) + T.(0.5) * sum(abs2,diff(x,dims=2),dims=[1,2])

end

"""
pcg_ml_est_fieldmap(y::AbstractMatrix{Complex{T}},z::AbstractMatrix{Complex{T}},β) 
Estimates the fieldmap using the formalism presented by Funai and Fessler (Regularized Field map Estimation, IEEE 2008)
# Arguments
* `y::AbstractMatrix{Complex{T}}` - Complex first-echo image data
* `z::AbstractMatrix{Complex{T}}` - Complex second-echo image data
* `β` - Regularization parameter controlling roughness penalty
"""
function pcg_ml_est_fieldmap(y::AbstractMatrix{Complex{T}},z::AbstractMatrix{Complex{T}},β) where T

    d = 1
    κ = conj.(y) .* z
    x = angle.(κ)
    m = abs.(κ) ./ maximum(abs.(κ))

    reltol = 5e-3
    c = 0
    Δ = 1
    itcount = 0

    trust_step = 1.0 ./ (m .+ 4*β)

    while Δ > reltol

        gs = gradient(Flux.params(x)) do
            c = ml_cost(x, y,z,β) # need to interpolate the y z and β to have better performance 
        end

        x̄ = gs[x]
        x .-= trust_step .* x̄ 

        Δ = abs.(ml_cost(x,y,z,β) - c)/c

        itcount +=1

    end

    print("Required $itcount iterations to converge below tolerance $reltol ")

    return x

end


"""
estimateB0Maps(imData,slices, TE1,TE2,β,isrotated)
Processes 3D volume data as output from MRIReco.reconstruction to estimate fieldmaps using the method presented by Funai and Fessler
# Arguments
* `imData` - 5-D array with complex image data
* `slices` - vector of slices to process (must be within range of 3rd dimension of imData)
* `TE1` - Echo time 1 [ms]
* `TE2` - Echo time 2 [ms]
* `β` - Regularization parameter controlling roughness penalty
* `isrotated` - Boolean controlling whether to rotate the B0 maps to match the images or not (legacy feature)
"""
function estimateB0Maps(imData,slices, TE1,TE2,β,isrotated)

    b0Maps = Complex.(zeros(size(imData)[1:3]))
    for slice in slices
        b0Maps[:,:,slice] = pcg_ml_est_fieldmap(imData[:,:,slice,1,1],imData[:,:,slice,2,1],β) ./ ((TE2 - TE1)/1000)
        println("for slice $slice")
    end

    b0Maps = mapslices(isrotated ? x->x : x-> rotl90(x),b0Maps,dims=(1,2)) 

    return real.(b0Maps)

end