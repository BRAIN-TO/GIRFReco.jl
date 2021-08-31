using DelimitedFiles, MRIReco, PyPlot, Dierckx, MAT, DSP

export GirfEssential, convertDomain!, time2freq, loadGirf, setIdentifier!

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

    # [nSamples] frequency vector [units Hz]
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

    @info "Converted GIRF from kHz to Hz"
    gFreq = gFreq.*1000

    @info "Loaded GIRF data\n"

    return GirfEssential(gData, gFreq, true, ["X", "Y", "Z"], [3])

end

## Function for building GIRF from the measurements taken by Tim in October 2020
# Takes in nothing, uses fixed filenames for now as we only have one set of GIRF measurements
# First sample of GIRF is an underflow, so use samples from 2:end
# Store the GIRF Fourier Representation in an Nx3 matrix, with second dimension indices 1,2,3 corresponding to directions x,y,z respectively
function buildGIRF()

    girf_filename_x = "data/GIRF/GIRF_X.mat"
    girf_filename_y = "data/GIRF/GIRF_Y.mat"
    girf_filename_z = "data/GIRF/GIRF_Z.mat"

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

    girf_filename_x = "data/GIRF/GIRF_X.mat"
    girf_filename_y = "data/GIRF/GIRF_Y.mat"
    girf_filename_z = "data/GIRF/GIRF_Z.mat"

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

        figure("Gx GIRF Magnitude (k1)")
        plot(GIRF_freq, abs.(GIRF_data[:,1]))
        xlim([-30,30])
        ylim([0.0, 1.05])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gy GIRF Magnitude (k1)")
        plot(GIRF_freq, abs.(GIRF_data[:,2]))
        xlim([-30,30])
        ylim([0.0, 1.05])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gz GIRF Magnitude (k1)")
        plot(GIRF_freq, abs.(GIRF_data[:,3]))
        xlim([-30,30])
        ylim([0.0, 1.05])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gx GIRF Phase (k1)")
        plot(GIRF_freq, angle.(GIRF_data[:,1]))
        xlim([-30,30])
        ylim([-pi, pi])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Phase")

        figure("Gy GIRF Phase (k1)")
        plot(GIRF_freq, angle.(GIRF_data[:,2]))
        xlim([-30,30])
        ylim([-pi, pi])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Phase")

        figure("Gz GIRF Phase (k1)")
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

    girf_filename_x = "../data/GIRF/GIRF_X.mat"
    girf_filename_y = "../data/GIRF/GIRF_Y.mat"
    girf_filename_z = "../data/GIRF/GIRF_Z.mat"

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

        figure("Gx GIRF Magnitude (k0)")
        plot(GIRF_freq, abs.(GIRF_data[:,1]))
        xlim([-3,3])
        ylim([0.0, 80])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gy GIRF Magnitude (k0)")
        plot(GIRF_freq, abs.(GIRF_data[:,2]))
        xlim([-3,3])
        ylim([0.0, 80])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gz GIRF Magnitude (k0)")
        plot(GIRF_freq, abs.(GIRF_data[:,3]))
        xlim([-3,3])
        ylim([0.0, 80])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Magnitude")

        figure("Gx GIRF Phase (k0)")
        plot(GIRF_freq, angle.(GIRF_data[:,1]))
        xlim([-3,3])
        ylim([-pi, pi])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Phase")


        figure("Gy GIRF Phase (k0)")
        plot(GIRF_freq, angle.(GIRF_data[:,2]))
        xlim([-3,3])
        ylim([-pi, pi])
        xlabel("Frequency [kHz]")
        ylabel("GIRF Phase")

        figure("Gz GIRF Phase (k0)")
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
