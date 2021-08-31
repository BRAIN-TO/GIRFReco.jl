
using MAT, PyPlot, DSP, Dierckx, FFTW

export predictGrad_port, buildGIRF
## Port of gradient impulse response function correction code by Johanna Vannesjo
#  https://github.com/MRI-gradient/GIRF
function predictGrad_port(freq_GIRF, GIRF, time_in, gradient_in, time_out)

    # Input variables
    # time_in: input time vector
    # gradient_in: input gradient
    # directions: gradient correction directions
    # convolution_type: type of convolution

    # Output variables
    # gradient_out: output gradient vector (sampled on time_out)
    # time_out: output time vector

    # check dimensions

    ns_GIRF = length(GIRF)
    df_GIRF = abs.(freq_GIRF[2] - freq_GIRF[1])

    # Set length of prediction

    T_GIRF = 1 ./df_GIRF
    T_O = max(time_in[end], time_out[end]) - min(time_in[1], time_out[1])
    T_ZF = min(T_GIRF, T_O)

    # prepare input

    dt_in = abs.(time_in[2] - time_in[1])
    nZF = Int(ceil(T_ZF./dt_in))

    gradient_in = vcat(zeros(nZF), gradient_in, zeros(nZF))
    time_in = vcat(time_in[1] .+ dt_in*(-nZF:-1), time_in, time_in[end] .+ dt_in*(1:nZF))

    f_in, df, f_max = time2freq(time_in)

    IN = fftshift(fft(gradient_in)).*dt_in
    ns_in = length(IN)
    T_in = time_in[end] .- time_in[1]

    # interpolate GIRF onto the input grid

    GIRF_real_spline = Spline1D(freq_GIRF,real.(GIRF), w=ones(length(freq_GIRF)), k=3, bc = "zero")
    GIRF_imaginary_spline = Spline1D(freq_GIRF, imag.(GIRF), w=ones(length(freq_GIRF)), k=3, bc = "zero")

    GIRF_ip = GIRF_real_spline(f_in) .+ 1im.*GIRF_imaginary_spline(f_in)

    OUT = GIRF_ip.*IN

    grad_OUT = real.(ifft(ifftshift(OUT))./dt_in)

    grad_OUT_spline = Spline1D(time_in,grad_OUT, w=ones(length(time_in)), k=3, bc = "zero")
    grad_OUT_resampled = grad_OUT_spline(time_out)

    figure("Frequency Spectra")
    plot(f_in,abs.(IN))
    plot(f_in,abs.(OUT))
    semilogy()

    figure("Gradient Waveform Comparison")
    plot(time_in,gradient_in)
    plot(time_in, grad_OUT)
    plot(time_out,grad_OUT_resampled)

    return grad_OUT_resampled

end

## Port of Johanna Vannesjo's time2freq method
#  https://github.com/MRI-gradient/GIRF/blob/main/utils/time2freq.m

function time2freq(t)

    nrs = length(t)
    dt = t[2] - t[1]
    f_max = 1 ./dt
    df = f_max ./ nrs
    f = ((0:nrs-1) .- floor(nrs/2)).*df

    return f, df, f_max

end

##
#Testing
# NEEDS: freq_GIRF, GIRF, time_in, gradient_in
time_in = collect(0.009*(0:100000-1))
gradient_in = sin.(time_in.*60000)
time_out = time_in

gx_filename = "C:\\Users\\ajaff\\UHN\\Brain-TO - MRP-GIRF - MRP-GIRF\\results\\GIRF_5usDwellTime\\GIRFGx_CoilCombined_Dwell5us.mat"

dataFile = matread(gx_filename)

df = 20 #Hz
freq_GIRF = df.*(2:length(dataFile["GIRF_FT"])) .- 100000

GIRF = dataFile["GIRF_FT"][2:end]

predictGrad_port(freq_GIRF,GIRF,time_in,gradient_in,time_out)

##
function buildGIRF()

    girf_filename_x = "GIRFReco/data/GIRF/GIRF_X.mat"
    girf_filename_y = "GIRFReco/data/GIRF/GIRF_Y.mat"
    girf_filename_z = "GIRFReco/data/GIRF/GIRF_Z.mat"

    GIRF_file_x = matread(girf_filename_x)
    GIRF_file_y = matread(girf_filename_y)
    GIRF_file_z = matread(girf_filename_z)

    GIRF_length = length(GIRF_file_x["GIRF_FT"]) .- 1

    GIRF_data = Matrix{ComplexF64}(undef, GIRF_length,3)

    GIRF_data[:,1] = GIRF_file_x["GIRF_FT"][2:end]
    GIRF_data[:,2] = GIRF_file_y["GIRF_FT"][2:end]
    GIRF_data[:,3] = GIRF_file_z["GIRF_FT"][2:end]

    return GIRF_data

end
