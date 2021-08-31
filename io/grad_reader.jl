using DelimitedFiles, MRIReco, PyPlot, Dierckx, MAT, DSP

export read_gradient_txt_file, buildGIRF, predictGrad_port, time2freq, buildGIRF_PN

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

    # SUPPRESS PLOTTING
    doPlot = false 

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

    # SUPPRESS PLOTTING
    doPlot = false

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

## Port of Johanna Vannesjo's time2freq method
#  https://github.com/MRI-gradient/GIRF/blob/main/utils/time2freq.m

function time2freq(t)

    nrs = length(t)
    dt = t[2] - t[1]
    f_max = 1 ./dt
    df = f_max ./ nrs
    f = ((0:nrs-1) .- floor(nrs/2.0)).*df

    return f, df, f_max

end

function read_gradient_txt_file(fileName, reconSize, delay)

    gradientData = readdlm(fileName,'\n')

    ## Read in the header data of the gradient text file (lines 1 to 21)
    dataDict = Dict{Symbol,Any}()
    dataDict[:versionNr] = gradientData[1]
    dataDict[:numSamples] = gradientData[2]
    dataDict[:dwellTime] = gradientData[3] # [seconds]
    dataDict[:samplesPerInterleave] = gradientData[4]
    dataDict[:numInterleaves] = gradientData[5]
    dataDict[:numDims] = gradientData[6]
    dataDict[:timeToCenterKSpace] = gradientData[7] # [seconds]
    dataDict[:acqDuration] = gradientData[8]
    dataDict[:samplesPerAcq] = gradientData[9]
    dataDict[:numAcquisitions] = gradientData[10]
    dataDict[:acqTR] = gradientData[11]
    dataDict[:gradientAcqStartDelay] = gradientData[12]
    dataDict[:echoTimeShiftSamples] = gradientData[13]
    dataDict[:FOV] = gradientData[14:16] # [m]
    dataDict[:voxelDims] = gradientData[17:19] # [m]
    dataDict[:gradientStrengthFactor] = gradientData[20] # [mT/m]
    dataDict[:isBinary] = gradientData[21]
    dataDict[:gamma] = 42577.478 # [Hz/mT] CAN CHANGE
    dataDict[:fieldStrength] = 3.0 # [T] CAN CHANGE

    #print(dataDict)

    ## reading and data scaling of gradient data
    dataDict[:gradData] = gradientData[22:end]
    interleaveGradArray = dataDict[:gradientStrengthFactor]*reshape(dataDict[:gradData],dataDict[:samplesPerInterleave],dataDict[:numInterleaves],dataDict[:numDims]) #[mT/m]

    plannedTimes = dataDict[:dwellTime].*(0:(dataDict[:samplesPerInterleave]-1))
    delayedTimes = plannedTimes .- delay .- dataDict[:dwellTime]./2 # seconds (dwellTime/2 compensates for integration)

    interleaveGradArrayFlexible = Array{Float64,3}(undef,size(interleaveGradArray))

    #print(size(interleaveGradArrayFlexible))

    ## Loop over all of the unique excitation trajectories and create an interpolant of the gradient
    for dim in 1:dataDict[:numDims]

        for l in 1:dataDict[:numInterleaves]

            #print((dim,l),"\n")

            sp = Spline1D(plannedTimes,interleaveGradArray[:,l,dim],w=ones(length(plannedTimes)), k=3, bc="zero", s=0.0)

            # evaluate the interpolant at the sampling times of the kspace data
            interleaveGradArrayFlexible[:,l,dim] = sp(delayedTimes)

            #print(interleaveGradArrayFlexible[:,l,dim][end],"\n")

        end

    end


    ## cumulative summation and numerical integration of the gradient data, resulting in the kspace trajectory
    kSpaceTrajArrayFlexible = dataDict[:gamma]*dataDict[:dwellTime]*cumsum(interleaveGradArrayFlexible,dims=1) # [rad/m]

    ## Conversion to the trajectory scaling convention in MRIReco.jl
    #  Currently only 2d Trajectories
    convertedKSpaceArrayFlexible = kSpaceTrajArrayFlexible

    ## This is hardcoded until I can fix the header data mutability :)
    convertedKSpaceArrayFlexible[:,:,1] *= 0.22/reconSize[1]
    convertedKSpaceArrayFlexible[:,:,2] *= 0.22/reconSize[2]

    ## Construction of the trajectory object ##

    ## Reshaping of the array to the format expected by the Trajectory constructor in MRIReco.jl
    # - dim 1 = kspace dimension
    # - dim 2 = kspace position (with interleaves/profiles arranged consecutively)
    permutedTrajectory = permutedims(reshape(convertedKSpaceArrayFlexible,dataDict[:samplesPerInterleave]*dataDict[:numInterleaves],dataDict[:numDims]),[2,1])

    ## Construction of the trajectory
    # - Note: timing vectors are automatically generated - seems to be consistent with the dwell time
    trajectoryObject = Trajectory(permutedTrajectory,dataDict[:numInterleaves],dataDict[:samplesPerInterleave],TE=dataDict[:echoTimeShiftSamples],AQ=dataDict[:acqDuration], numSlices=9, cartesian=false,circular=false)

end

##
function read_gradient_txt_file(fileName, reconSize, delay, doGIRF)

    gradientData = readdlm(fileName,'\n')

    ## Read in the header data of the gradient text file (lines 1 to 21)
    dataDict = Dict{Symbol,Any}()
    dataDict[:versionNr] = gradientData[1]
    dataDict[:numSamples] = gradientData[2]
    dataDict[:dwellTime] = gradientData[3] # [seconds]
    dataDict[:samplesPerInterleave] = gradientData[4]
    dataDict[:numInterleaves] = gradientData[5]
    dataDict[:numDims] = gradientData[6]
    dataDict[:timeToCenterKSpace] = gradientData[7] # [seconds]
    dataDict[:acqDuration] = gradientData[8]
    dataDict[:samplesPerAcq] = gradientData[9]
    dataDict[:numAcquisitions] = gradientData[10]
    dataDict[:acqTR] = gradientData[11]
    dataDict[:gradientAcqStartDelay] = gradientData[12]
    dataDict[:echoTimeShiftSamples] = gradientData[13]
    dataDict[:FOV] = gradientData[14:16] # [m]
    dataDict[:voxelDims] = gradientData[17:19] # [m]
    dataDict[:gradientStrengthFactor] = gradientData[20] # [mT/m]
    dataDict[:isBinary] = gradientData[21]
    dataDict[:gamma] = 42577.478 # [Hz/mT]
    dataDict[:fieldStrength] = 3.0 # [T]

    #print(dataDict)

    ## reading and data scaling of gradient data
    dataDict[:gradData] = gradientData[22:end]
    interleaveGradArray = dataDict[:gradientStrengthFactor]*reshape(dataDict[:gradData],dataDict[:samplesPerInterleave],dataDict[:numInterleaves],dataDict[:numDims]) #[mT/m]

    #print(size(interleaveGradArray))

    ## Set timings to allow for delay compensation
    plannedTimes = dataDict[:dwellTime].*(0:(dataDict[:samplesPerInterleave]-1))
    delayedTimes = plannedTimes .- delay .- dataDict[:dwellTime]./2 # seconds

    interleaveGradArrayUpdated = Array{Float64,3}(undef,size(interleaveGradArray))

    k0_Array = Array{Float64,3}(undef,size(interleaveGradArray))
    k0_combined = Array{Float64,2}(undef,size(k0_Array,1), size(k0_Array,2))

    if doGIRF

        frequencyVec, GIRF = buildGIRF_PN()

        frequencyVec = frequencyVec .*1000

        figure("RAW GIRF")
        plot(frequencyVec, abs.(GIRF[:,2]))
        xlim([-8000,8000])
        ylim([0.85,1.02])


        for k = 1:dataDict[:numInterleaves]

            interleaveGradArrayUpdated[:,k,1] = predictGrad_port(frequencyVec,GIRF[:,1], plannedTimes, interleaveGradArray[:,k,1],  plannedTimes)
            interleaveGradArrayUpdated[:,k,2] = predictGrad_port(frequencyVec,GIRF[:,2], plannedTimes, interleaveGradArray[:,k,2],  plannedTimes)

        end

    else

        interleaveGradArrayUpdated = interleaveGradArray

    end

    interleaveGradArrayFlexible = Array{Float64,3}(undef,size(interleaveGradArray))

    #print(size(interleaveGradArrayFlexible))

    ## Loop over all of the unique excitation trajectories and create an interpolant of the gradient
    for dim in 1:dataDict[:numDims]

        for l in 1:dataDict[:numInterleaves]

            #print((dim,l),"\n")

            sp = Spline1D(plannedTimes,interleaveGradArrayUpdated[:,l,dim],w=ones(length(plannedTimes)), k=3, bc="zero", s=0.0)

            # evaluate the interpolant at the sampling times of the kspace data
            interleaveGradArrayFlexible[:,l,dim] = sp(delayedTimes)

            #print(interleaveGradArrayFlexible[:,l,dim][end],"\n")



        end

    end

    doK0 = true

    if doK0

        figure("K0")
        plot(delayedTimes, interleaveGradArrayFlexible[:,1,1] .+ interleaveGradArrayFlexible[:,1,2])
        xlabel("Time (s)")
        ylabel("k₀ Fluctuation (A.U)")

    end

    if doK0

        frequencyVecK0, GIRF_K0 = buildGIRF_K0()

        frequencyVecK0 = frequencyVecK0 .*1000

        figure("RAW GIRF")
        plot(frequencyVecK0, abs.(GIRF_K0[:,2]))
        xlim([-8000,8000])
        ylim([0.85,1.02])

        for k = 1:dataDict[:numInterleaves]

            k0_Array[:,k,1] = predictGrad_port(frequencyVecK0,GIRF_K0[:,1], plannedTimes, interleaveGradArrayUpdated[:,k,1],  plannedTimes)
            k0_Array[:,k,2] = predictGrad_port(frequencyVecK0,GIRF_K0[:,2], plannedTimes, interleaveGradArrayUpdated[:,k,2],  plannedTimes)

        end

        k0_combined = (k0_Array[:,:,1] .+ k0_Array[:,:,2])./1000

        figure("Rad/s Fluctuation")
        plot(k0_Array[:,:,1] .+ k0_Array[:,:,2])

        figure("Rad Fluctuation k₀")
        plot(k0_combined)

    end

    ## cumulative summation and numerical integration of the gradient data, resulting in the kspace trajectory
    kSpaceTrajArrayFlexible = dataDict[:gamma]*dataDict[:dwellTime]*cumsum(interleaveGradArrayFlexible,dims=1) # [rad/m]

    ## Conversion to the trajectory scaling convention in MRIReco.jl
    #  Currently only 2d Trajectories
    convertedKSpaceArrayFlexible = kSpaceTrajArrayFlexible

    ## This is hardcoded until I can fix the header data mutability :)
    convertedKSpaceArrayFlexible[:,:,1] *= 0.22/reconSize[1]
    convertedKSpaceArrayFlexible[:,:,2] *= 0.22/reconSize[2]

    ## Construction of the trajectory object ##

    ## Reshaping of the array to the format expected by the Trajectory constructor in MRIReco.jl
    # - dim 1 = kspace dimension
    # - dim 2 = kspace position (with interleaves/profiles arranged consecutively)
    permutedTrajectory = permutedims(reshape(convertedKSpaceArrayFlexible,dataDict[:samplesPerInterleave]*dataDict[:numInterleaves],dataDict[:numDims]),[2,1])

    ## Construction of the trajectory
    # - Note: timing vectors are automatically generated - seems to be consistent with the dwell time
    trajectoryObject = Trajectory(permutedTrajectory,dataDict[:numInterleaves],dataDict[:samplesPerInterleave],TE=dataDict[:echoTimeShiftSamples],AQ=dataDict[:acqDuration], numSlices=9, cartesian=false,circular=false)

    return trajectoryObject, k0_combined

end

## Testing
gradFile = "C:/Users//ajaff/Documents/BrainTO/Brain-TO - MRP-SPIDI - SPIDI/data/SPIDI_0003/gradients/gradients523.txt"

##
kSpaceTrajectory, k0_phase = read_gradient_txt_file(gradFile, (200,200), 0.00000, true)
kSpaceTrajectory_2 = read_gradient_txt_file(gradFile,(200,200),0.00000)

##
pulledTrajectory11 = kspaceNodes(kSpaceTrajectory)[1,1:3100]
pulledTrajectory12 = kspaceNodes(kSpaceTrajectory)[2,1:3100]

pulledTrajectory21 = kspaceNodes(kSpaceTrajectory_2)[1,1:3100]
pulledTrajectory22 = kspaceNodes(kSpaceTrajectory_2)[2,1:3100]

#
fig = figure(234, figsize=(10,10))
ax = fig.gca()
ax.scatter(pulledTrajectory11,pulledTrajectory12, label="GIRF Corrected")
ax.scatter(pulledTrajectory21,pulledTrajectory22, label="Nominal")
xlabel("kx")
ylabel("ky")
title("K-space Center")
xlim((-0.05,0.05))
ylim((-0.05,0.05))
legend()




#ax.plot(pulledTrajectory31[1:3100],pulledTrajectory32[1:3100])
#ax.scatter(pulledTrajectory21[1:3100],pulledTrajectory22[1:3100])
#ax.scatter(pulledTrajectory41[1:3100],pulledTrajectory42[1:3100])

# ax.set_aspect(aspect=1.0)
