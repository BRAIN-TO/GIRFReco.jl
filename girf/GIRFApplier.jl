using DelimitedFiles, MRIReco, PyPlot, Dierckx, MAT, DSP

export applyGIRF

include("GIRFEssential.jl")

## Mutable struct definition for GirfApplier type: following convention of JV in MATLAB repo
# composite type of GirfEssential and a field for gamma, the gyromagnetic ratio
mutable struct GirfApplier

    # PROPERTIES (fields)

    # GirfEssential type containing majority of GIRF data
    essential::GirfEssential

    # Gamma (gyromagnetic ratio)
    gamma::Float64

end

## Function for predicting the gradient by applying the GIRF
# This is done in one dimension only. Repeat this function for multiple dimensions
# Inputs:
#   freq_GIRF: Vector (sampling frequencies for the GIRF data)
#   GIRF: Vector (GIRF data in the frequency domain)
#   time_in: Vector (Sampling times for the time domain gradient data)
#   gradient_in: Vector (time domain gradient data sampled at times in time_in)
#   time_out: Vector (desired sampling times for the output gradient)

# Outputs:
#   grad_OUT_resampled: Vector (filtered gradient data sampled at time_out)

function applyGIRF(g::GirfApplier, gradient_in, time_in, time_out, direction)

    # Input variables
    # time_in: input time vector
    # gradient_in: input gradient
    # directions: gradient correction directions
    # convolution_type: type of convolution

    # Output variables
    # gradient_out: output gradient vector (sampled on time_out)
    # time_out: output time vector

    # GET GIRF DATA FROM THE GIRF DATATYPE
    freq_GIRF = g.essential.freq
    GIRF = g.essential.girf[:,direction]

    # check dimensions
    ns_GIRF = length(GIRF)
    df_GIRF = abs.(freq_GIRF[2] - freq_GIRF[1])

    # make copies so that original variables and objects are not mutated
    gradient_in_original = deepcopy(gradient_in)
    time_in_original = deepcopy(time_in)

    t_or1, t_or2, t_or3 = time2freq(time_in_original)

    figure("RAW FFT")
    plot(t_or1, abs.(fftshift(fft(gradient_in_original))))
    semilogy()

    # Set length of prediction
    T_GIRF = 1 ./df_GIRF # Signal Length for nyquist fs = 2B and T = 1/2B
    T_O = max(time_in[end], time_out[end]) - min(time_in[1], time_out[1])
    T_ZF = min(T_GIRF, T_O) # Zero pad to the shortest time of the two

    # prepare input in time
    dt_in = abs.(time_in[2] - time_in[1])
    nZF = Int(ceil(T_ZF./dt_in)) # test this out

    # Do Zero filling
    gradient_in = vcat(zeros(nZF),  gradient_in, zeros(nZF))
    time_in = vcat(time_in[1] .+ dt_in*(-nZF:-1), time_in, time_in[end] .+ dt_in*(1:nZF))

    # get updated frequency vector corresponding to Zero padded time vector
    f_in, df, f_max = time2freq(time_in)

    # Get frequency domain gradient
    IN = fftshift(fft(gradient_in)).*dt_in
    ns_in = length(time_in)

    # interpolate GIRF onto the input grid
    GIRF_real_spline = Spline1D(freq_GIRF,real.(GIRF), w=ones(length(freq_GIRF)), k=3, bc = "zero")
    GIRF_imaginary_spline = Spline1D(freq_GIRF, imag.(GIRF), w=ones(length(freq_GIRF)), k=3, bc = "zero")

    # recombine interpolated reals and imaginary parts
    GIRF_ip = GIRF_real_spline(f_in) .+ 1im.*GIRF_imaginary_spline(f_in)

    ## filter the output to avoid aliasing
    figure("Magnitude FFT of Response")
    windowFunction = fftshift(tukey(length(GIRF_ip), 0.5;zerophase = true))
    plot(windowFunction)
    plot(GIRF_ip)
    GIRF_ip = windowFunction .* GIRF_ip
    plot(GIRF_ip)

    # Convolution in Time (mult in freq) to get output Gradient train in frequency domain
    OUT = GIRF_ip.*IN

    grad_OUT = real.(ifft(ifftshift(OUT))./dt_in)

    grad_OUT_spline = Spline1D(time_in,grad_OUT, w=ones(length(time_in)), k=3, bc = "zero")
    grad_OUT_resampled = grad_OUT_spline(time_out)

    figure("Difference between corrected and uncorrected gradients")
    plot(time_out, grad_OUT_resampled .- gradient_in_original)
    xlabel("Time [s]")
    ylabel("Pointwise Gradient Difference [a.u]")

    return grad_OUT_resampled

end
