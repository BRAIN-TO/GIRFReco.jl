using DelimitedFiles, MRIReco, PyPlot, Dierckx, MAT, DSP

export GirfEssential, convertDomain!, time2freq, loadGirf, formatGirf, setIdentifier!

include("Util.jl")

## GirfEssential
#  Struct definition as per Johanna's GirfEssential class in GIRF MATLAB repo
mutable struct GirfEssential

    ## INDEPENDENT PROPERTIES (i.e. PROPERTIES SET AT INSTANTIATION)

    # GIRF identifier or tag
    name::String

    # is GIRF calculated in freq or time domain?
    isFreqDomainGirf::Bool

    # Input gradient and shim channel names
    inChannels::AbstractVector

    # Output field basis
    outBasis::AbstractVecOrMat

    # [nSamples x nOutBasis x nInChannels] Calculated GIRF (self and cross-terms)
    girf::AbstractVecOrMat

    # [nSamples] frequency vector
    freq::AbstractVector

    # [nSamples x nOutBasis x nInChannels] if GIRF is in the time domain
    girfTime::AbstractVecOrMat

    # [nSamples] time vector
    time::AbstractVector

    ## DEPENDENT PROPERTIES (REQUIRES INSTANTIATION OF ABOVE PROPERTIES)

    # Sampling frequency
    df::Float64

    # Sampling time
    dt::Float64

    # Define self basis
    selfBasis::AbstractVector

end

## TODO
#  (Continue) Standard Method Overloading for GIRF

Base.size(g::GirfEssential) = size(g.girf)
Base.size(g::GirfEssential, i::Int) = size(g.girf,i)
Base.string(g::GirfEssential) = g.name
Base.print(g::GirfEssential) = displayGirf(g)
Base.length(g::GirfEssential) = size(g.girf, 1)

## TODO
#  Continue Get and Set methods
isFreqDomain(g::GirfEssential) = g.isFreqDomainGirf
numOutBasis(g::GirfEssential) = size(g.girf,2)
numInChannels(g::GirfEssential) = size(g.girf, 3)

## Function to create a girf structure with the same interface as the constructor in the MATLAB implementation
function GirfEssential(data::AbstractVecOrMat, vect::AbstractVector, isFreq::Bool, inChannels, outBasis)

    # Define unused and auxiliary parameters
    # inChannels = ["Empty"]
    # outBasis = [0]
    selfBasis = [0]
    girfTime = [0]
    time = [0]
    df = 0
    dt = 0
    identifier = "Identifier Not Set Yet"

    # Call to standard constructor
    g = GirfEssential(identifier, isFreq, inChannels, outBasis, data, vect, girfTime, time, df, dt, selfBasis)

    convertDomain!(g)

    @info "Created GirfEssential Structure with the parameters as follows:\n"
    print(g)

    return g

end

## Convert frequency-domain GIRF to time domain and vice versa
#
# IN
# domain    Specify conversion domains {'freq2time' 'time2freq'}
#           Default is determined by this.isFreqDomainGirf
#
# OUT
#
# EXAMPLE
#   convertDomain!(g);
#
#   See also GirfEssential
#
# Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
# Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
#               2016 FMRIB centre, University of Oxford
#
# This file is part of a code package for GIRF computation and application.
# The package is available under a BSD 3-clause license. Further info see:
# https://github.com/MRI-gradient/girf
#
# TODO edit comment to reflect in-place function on girf

function convertDomain!(g::GirfEssential)

    if !g.isFreqDomainGirf

        # Convert time domain GIRF to frequency domain and fill the girf struct properly
        g.girfTime = deepcopy(g.girf)
        g.girf = fftshift(fft(ifftshift(g.girfTime)))
        g.time = deepcopy(g.freq)
        g.dt = g.time[2] .- g.time[1]
        g.freq, g.df, f_max = time2freq(g.time)

        @info "Time-domain GIRF provided, converting to Frequency-domain"

    else

        # Convert frequency domain GIRF to time domain and fill girf struct properly
        g.girfTime = real(fftshift(ifft(fftshift(g.girf))))
        g.df = g.freq[2] .- g.freq[1]
        g.time, g.dt, t_max = time2freq(g.freq)

        @info "Frequency-domain GIRF provided, converting to Time-domain"

    end

    return g

end

## Function to name the GIRF with an identifier
function setIdentifier!(g::GirfEssential, identifier::String)

    g.name = identifier

    @info "Name has been set as: $identifier\n"

end


## Function to load girf data from file into GirfEssential object
#
# IN
# filename  Name of file to load data from
#
# OUT
#
# EXAMPLE
#   girfE.Load(mySavedGirfFilename);
#
#   See also GirfEssential
#
# Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
# Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
#               2016 FMRIB centre, University of Oxford
#
# This file is part of a code package for GIRF computation and application.
# The package is available under a BSD 3-clause license. Further info see:
# https://github.com/MRI-gradient/girf
#
# TODO

function loadGirf(g::GirfEssential, filename::String)

    # gFreq, gData = buildGIRF_PN()

    @info "Loaded GirfEssential data from $filename"

end

## Function to save data in GirfEssential object to file
#
# IN
# filename  Name of file to save data into
# overwrite Overwrite existing file?
#
# OUT
#
# EXAMPLE
#   girfE.Save(mySaveGirfFilename);
#
#   See also GirfEssential
#
# Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
# Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
#               2016 FMRIB centre, University of Oxford
#
# This file is part of a code package for GIRF computation and application.
# The package is available under a BSD 3-clause license. Further info see:
# https://github.com/MRI-gradient/girf
#
# TODO

function saveGirf(g::GirfEssential, filename::String)

    @info "Saved GirfEssential to $filename"

end

## Function for loading (AJAFFRAY Pre-Built)
function loadGirf(degree = 1)

    if degree == 1
        gFreq, gData = buildGIRF_PN()
    elseif degree == 0
        gFreq, gData = buildGIRF_K0()
    else
        @error "Tried to load unexpected GIRF degree"
    end

    @info "Loaded GIRF data\n"

    return GirfEssential(gData, gFreq, true, ["X", "Y", "Z"], [3])

end

## Function for building GIRF from the measurements taken by Tim in October 2020
# Takes in nothing, uses fixed filenames for now as we only have one set of GIRF measurements
# First sample of GIRF is an underflow, so use samples from 2:end
# Store the GIRF Fourier Representation in an Nx3 matrix, with second dimension indices 1,2,3 corresponding to directions x,y,z respectively
function buildGIRF()

    girf_filename_x = "C:\\Users\\ajaff\\UHN\\Brain-TO - MRP-GIRF - MRP-GIRF\\results\\GIRF_5usDwellTime\\GIRFGx_CoilCombined_Dwell5us.mat"
    girf_filename_y = "C:\\Users\\ajaff\\UHN\\Brain-TO - MRP-GIRF - MRP-GIRF\\results\\GIRF_5usDwellTime\\GIRFGy_CoilCombined_Dwell5us.mat"
    girf_filename_z = "C:\\Users\\ajaff\\UHN\\Brain-TO - MRP-GIRF - MRP-GIRF\\results\\GIRF_5usDwellTime\\GIRFGz_CoilCombined_Dwell5us.mat"

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

## Function for building GIRF from the measurements taken by Tim in June 2021
# Takes in nothing, uses fixed filenames for now as we only have one set of GIRF measurements
# First sample of GIRF is an underflow, so use samples from 2:end
# Store the GIRF Fourier Representation in an Nx3 matrix, with second dimension indices 1,2,3 corresponding to directions x,y,z respectively
# Store the frequencies corresponding to the GIRF as a separate output vector

function buildGIRF_PN(doPlot = true, doFiltering = true)

    girf_filename_x = "E:\\MRP_GIRF\\results\\GIRF_X.mat"
    girf_filename_y = "E:\\MRP_GIRF\\results\\GIRF_Y.mat"
    girf_filename_z = "E:\\MRP_GIRF\\results\\GIRF_Z.mat"

    GIRF_file_x = matread(girf_filename_x)
    GIRF_file_y = matread(girf_filename_y)
    GIRF_file_z = matread(girf_filename_z)

    GIRF_length = length(GIRF_file_x["GIRF_FT"]) .- 1

    GIRF_data = Matrix{ComplexF64}(undef, GIRF_length,3)

    GIRF_data[:,1] = GIRF_file_x["GIRF_FT"][2:end]
    GIRF_data[:,2] = GIRF_file_y["GIRF_FT"][2:end]
    GIRF_data[:,3] = GIRF_file_z["GIRF_FT"][2:end]
    GIRF_freq = GIRF_file_z["freq"][2:end]

    if doFiltering

        window = tukey(GIRF_length, 0.25, zerophase = false)

        for l = 1:3

            GIRF_data[:,l] = window .* GIRF_data[:,l]

        end

    end

    if doPlot

        figure("Gx GIRF Magnitude")
        plot(GIRF_freq, abs.(GIRF_data[:,1]))
        xlim([-30,30])
        ylim([0.0, 1.05])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gy GIRF Magnitude")
        plot(GIRF_freq, abs.(GIRF_data[:,2]))
        xlim([-30,30])
        ylim([0.0, 1.05])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gz GIRF Magnitude")
        plot(GIRF_freq, abs.(GIRF_data[:,3]))
        xlim([-30,30])
        ylim([0.0, 1.05])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gx GIRF Phase")
        plot(GIRF_freq, angle.(GIRF_data[:,1]))
        xlim([-30,30])
        ylim([-pi, pi])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Phase")

        figure("Gy GIRF Phase")
        plot(GIRF_freq, angle.(GIRF_data[:,2]))
        xlim([-30,30])
        ylim([-pi, pi])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Phase")

        figure("Gz GIRF Phase")
        plot(GIRF_freq, angle.(GIRF_data[:,3]))
        xlim([-30,30])
        ylim([-pi,pi])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Phase")

    end

    return GIRF_freq, GIRF_data

end

## Function for building GIRF from the measurements taken by Tim in June 2021
# Takes in nothing, uses fixed filenames for now as we only have one set of GIRF measurements
# First sample of GIRF is an underflow, so use samples from 2:end
# Store the GIRF Fourier Representation in an Nx3 matrix, with second dimension indices 1,2,3 corresponding to directions x,y,z respectively
# Store the frequencies corresponding to the GIRF as a separate output vector

function buildGIRF_K0(doPlot = true, doFiltering = true)

    girf_filename_x = "E:\\MRP_GIRF\\results\\GIRF_X.mat"
    girf_filename_y = "E:\\MRP_GIRF\\results\\GIRF_Y.mat"
    girf_filename_z = "E:\\MRP_GIRF\\results\\GIRF_Z.mat"

    GIRF_file_x = matread(girf_filename_x)
    GIRF_file_y = matread(girf_filename_y)
    GIRF_file_z = matread(girf_filename_z)

    GIRF_length = length(GIRF_file_x["b0ec_FT"]) .- 1

    GIRF_data = Matrix{ComplexF64}(undef, GIRF_length,3)

    GIRF_data[:,1] = GIRF_file_x["b0ec_FT"][2:end]
    GIRF_data[:,2] = GIRF_file_y["b0ec_FT"][2:end]
    GIRF_data[:,3] = GIRF_file_z["b0ec_FT"][2:end]
    GIRF_freq = GIRF_file_z["freq"][2:end]

    if doFiltering

        window = tukey(GIRF_length, 0.25, zerophase = false)

        for l = 1:3

            GIRF_data[:,l] = window .* GIRF_data[:,l]

        end

    end

    # ADD PREPROCESSING TO REMOVE HIGH FREQUENCY NOISE

    if doPlot

        figure("Gx GIRF Magnitude")
        plot(GIRF_freq, abs.(GIRF_data[:,1]))
        xlim([-3,3])
        ylim([0.0, 80])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gy GIRF Magnitude")
        plot(GIRF_freq, abs.(GIRF_data[:,2]))
        xlim([-3,3])
        ylim([0.0, 80])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gz GIRF Magnitude")
        plot(GIRF_freq, abs.(GIRF_data[:,3]))
        xlim([-3,3])
        ylim([0.0, 80])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gx GIRF Phase")
        plot(GIRF_freq, angle.(GIRF_data[:,1]))
        xlim([-3,3])
        ylim([-pi, pi])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Phase")


        figure("Gy GIRF Phase")
        plot(GIRF_freq, angle.(GIRF_data[:,2]))
        xlim([-3,3])
        ylim([-pi, pi])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Phase")

        figure("Gz GIRF Phase")
        plot(GIRF_freq, angle.(GIRF_data[:,3]))
        xlim([-3,3])
        ylim([-pi,pi])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Phase")

    end

    return GIRF_freq, GIRF_data

end


## TODO
#  function for setting the selfBasis field in the GirfEssential struct
function setSelfBasis(g::GirfEssential)

    #Implement the method selfBasis written in lines 84 to 99 of GirfEssential.m in Johanna's GIRF code

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

function applyGIRF(g::GirfEssential, gradient_in, time_in, time_out)

    # Input variables
    # time_in: input time vector
    # gradient_in: input gradient
    # directions: gradient correction directions
    # convolution_type: type of convolution

    # Output variables
    # gradient_out: output gradient vector (sampled on time_out)
    # time_out: output time vector

    # GET GIRF DATA FROM THE GIRF DATATYPE
    freq_GIRF = g.freq
    GIRF = g.girf


    # check dimensions
    ns_GIRF = length(GIRF)
    df_GIRF = abs.(freq_GIRF[2] - freq_GIRF[1])
    #@info df_GIRF

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
    #@info T_GIRF
    #@info T_ZF

    # prepare input in time
    dt_in = abs.(time_in[2] - time_in[1])
    nZF = Int(ceil(T_ZF./dt_in)) # test this out
    #@info nZF

    # Do Zero filling
    gradient_in = vcat(zeros(nZF),  gradient_in, zeros(nZF))
    time_in = vcat(time_in[1] .+ dt_in*(-nZF:-1), time_in, time_in[end] .+ dt_in*(1:nZF))

    # WRITE GRADIENTS TO MATLAB FILE FOR EXTERNAL TESTING
    # testDictionary = Dict{String,Any}()
    # testDictionary["grads"] = gradient_in
    #
    # matwrite("grad_in.mat",testDictionary)

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
    figure("Response")
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
    plot(grad_OUT_resampled .- gradient_in_original)

    return grad_OUT_resampled

end

##  Function for displaying girf in terminal reasonably
function displayGirf(g::GirfEssential)

    # @info "Pretty-print NOT implemented yet!"

    println(g.name)
    println(g.isFreqDomainGirf)
    println(g.inChannels)
    println(g.outBasis)

    if !isempty(g.girf)
        println("GIRF in Frequency Domain Set Up")
    else
        println("GIRF Empty")
    end

    println(g.freq[1]:g.df:g.freq[end])

    if !isempty(g.girfTime)
        println("GIRF in TIme Domain Set Up")
    else
        println("GIRF Empty")
    end

    println(g.time[1]:g.dt:g.time[end])
    println(g.dt)
    println(g.df)

end
