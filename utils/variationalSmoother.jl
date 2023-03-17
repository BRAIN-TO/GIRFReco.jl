using Flux

"""
cost_smooth(x::Matrix{T},y::Matrix{Complex{T}}, β)
Calculates the cost between an estimated sensitivity map and the unsmoothed map
# Arguments
* `x::Matrix{Complex{T}}` - smoothed map
* `y::Matrix{Complex{T}}` - unsmoothed map
* `β` - Regularization parameter controlling roughness penalty
"""
function cost_smooth(x::Matrix{T},y::Matrix{T}, β) where T

    (0.5 .* sum(abs2, x - y, dims=[1,2]) .+ β*R(x))[1]

end

"""
R(x::Matrix{T})
Regularization function which penalizes roughness
# Arguments
* `x::Matrix{Complex{T}}` - fieldmap estimate (in radians)
"""
function R(x::Matrix{T}) where T

    0.5 * sum(abs2,diff(x,dims=1),dims=[1,2]) + 0.5 * sum(abs2,diff(x,dims=2),dims=[1,2])

end

"""
sens_smooth(x::AbstractMatrix{Complex{T}},y::AbstractMatrix{Complex{T}},β) 
Estimates smooth sensitivity maps using gradient descent and a roughness penalty
# Required Arguments
* `x::AbstractMatrix{Complex{T}}` - Complex smoothed sensitivity map
* `y::AbstractMatrix{Complex{T}}` - Complex reference map

# Optional Arguments
* `β` - Regularization parameter controlling roughness penalty
* `reltol` - early stopping criteria (exit if subsequent cost function change < reltol)
"""
function sens_smooth(x_c::AbstractMatrix{Complex{T}},y_c::AbstractMatrix{Complex{T}},β = 1e-1, reltol = 5e-5) where T

    c_r = 0
    c_i = 0
    Δ = 1
    itcount = 0

    x_r = real.(x_c)
    y_r = real.(y_c)
    x_i = imag.(x_c)
    y_i = imag.(y_c)

    while Δ > reltol && itcount < 1000

        gs_r = gradient(Flux.params(x_r)) do
            c_r = cost_smooth(x_r, y_r,β) 
        end

        gs_i = gradient(Flux.params(x_i)) do
            c_i = cost_smooth(x_i, y_i,β)
        end

        x̄_r = gs_r[x_r]
        x_r .-=  0.01 .* x̄_r 

        x̄_i = gs_i[x_i]
        x_i .-=  0.01 .* x̄_i 

        Δ_r = abs.(cost_smooth(x_r,y_r,β) - c_r)/(c_r + eps())
        Δ_i = abs.(cost_smooth(x_i,y_i,β) - c_i)/(c_i + eps())

        Δ = sqrt.(abs2.(Δ_r) + abs2.(Δ_i))

        @info Δ

        itcount +=1

    end

    print("Required $itcount iterations to converge below tolerance $reltol ")

    return complex.(x_r,x_i)

end

"""
smooth_sensitivity_maps(smaps;β = 5e-4, reltol = 0.001)
Processes 3D volume data as output from MRIReco.reconstruction to estimate fieldmaps using the method presented by Funai and Fessler
# Required Arguments
* `smaps` - 4-D array with sensitivity maps (x,y,z,coils) 

# Keyword Arguments 
* `β` - Regularization parameter controlling roughness penalty (larger = smoother)
* `reltol` - early stopping criteria (exit if subsequent cost function change < reltol
"""
function smooth_sensitivity_maps(smaps; β = 5e-4, reltol = 0.001)

    nx,ny,slices,coils = size(smaps)
    smaps_ref = smaps
    smaps_out = Complex.(zeros(nx,ny,slices,coils))

    for coil in 1:coils
        @sync for slice in 1:slices
            Threads.@spawn begin
                smaps_out[:,:,slice,coil] = sens_smooth(smaps[:,:,slice,coil],smaps_ref[:,:,slice,coil],β,reltol)
                println("for slice $slice") 
            end
        end
        println("for coil $coil")
    end

    return smaps_out

end


function process_sens_maps(smap)

    d = size(smap)
    psize = 4
    thres = 0.005

    mask = collect(padarray(abs.(smap) .> thres, Fill(0, (psize,psize))))

    for i = 1:8

        mask = erode(mask)

    end

    mask = mask[(psize+1):(psize+d[1]),(psize+1):(psize+d[2])]

    mask = get_conv_hull_mask(mask)

    return 1im .* (mask) .* abs.(smap) .* exp.(1im .* angle.(smap)).^mask

end

function get_conv_hull_mask(mask)

    cord = convexhull(mask)
    nedge = length(cord)
    edges = hcat(1:nedge,2:(nedge+1))
    edges[end] = 1
    tol = 1e-1

    stat = inpoly2(vec(CartesianIndices((1:size(mask,1),1:size(mask,2)))),cord,edges,atol=tol)

    return reshape(stat[:,1],size(mask))

end

function standardize_sens_phase(sens)

    for i = 1:size(sens,4)

        sens[:,:,:,i] = sens[:,:,:,i] .* exp.(-1im .*mean(angle.(sensitivity[:,:,:,i])))

    end

    return sens

end