using PyPlot, HDF5, MRIReco, LinearAlgebra, Dierckx, DSP, FourierTools, ImageBinarization, ImageEdgeDetection

## General Plotting function for the reconstruction

# Mosaic-plots reconstruction for selected slices and corresponding B0 map
function plotReconstruction(images, slices, b0)

    # Slice ordering check (show the correct order of slices as the images are read in not in geometrically sequential order)
    indexArray = slices

    # Plot magnitude images (normalize)
    figure("Magnitude Images")
    absData = abs.(images[:,:,slices,1,1])./maximum(abs.(images[:,:,slices,1,1]))
    absMosaic = mosaicview(absData, nrow=Int(floor(sqrt(length(slices)))),npad=5,rowmajor=true, fillvalue=0)

    PyPlot.imshow(absMosaic,cmap="gray")
    colorbar()

    gcf().suptitle("|Images|")

    # Plot phase images
    figure("Phase Images")

    phaseData = mapslices(x ->ROMEO.unwrap(x),angle.(images[:,:,slices,1,1]),dims=(1,2))
    phaseMosaic = mosaicview(phaseData,nrow=Int(floor(sqrt(length(slices)))),npad=5,rowmajor=true, fillvalue=0)

    PyPlot.imshow(phaseMosaic, cmap="inferno",vmax = 3*pi,vmin=-3*pi)
    colorbar()

    gcf().suptitle("∠Images")

    # Plot B0 maps
    figure("B₀ Map Images")

    b0Mosaic = mosaicview(b0[:,:,slices],nrow=Int(floor(sqrt(length(slices)))),npad = 5, rowmajor=true, fillvalue=0)

    PyPlot.imshow(b0Mosaic, cmap="inferno",vmax=500,vmin=-500)
    colorbar()

    gcf().suptitle("B₀ Maps")

end

# Function plots all profiles in the acquisition to check consistency with ISMRMRD file
function checkProfiles(rawData)

    numProfiles2 = 128 # Set to the number of profiles that you would like to see

    for l = 1:numProfiles2
        figure("Profile $l")
        subplot(2,1,1)
        plot(abs.(rawData.profiles[l].data[:,1]))
        subplot(2,1,2)
        plot(angle.(rawData.profiles[l].data[:,1]))
    end

end


# Create figure and plot the sensitivity maps for each coil composed as a collage
function plotSenseMaps(sense,n_channels)
    # Magnitude maps
    figure("Sensitivity Map Magnitude"); clf(); for ch in 1:n_channels; subplot(8,4,ch); PyPlot.imshow((abs.(sense[:,:,1,ch]))); end;
    subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
    gcf()

    # Phase maps
    figure("Sensitivity Map Phase"); clf(); for ch in 1:n_channels; subplot(8,4,ch); PyPlot.imshow(ROMEO.unwrap(angle.(sense[:,:,1,ch]))); end;
    subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
    gcf()

end

# TODO?
function plotTrajAndData(acq)

    for l in 1:length(acq.traj)

        freqEncode[l,:] = acq.traj[l].nodes[1,:]
        phaseEncode[l,:] = acq.traj[l].nodes[2,:]
        kSignal[l,:] = acq.kdata[l,:,1]

    end

end

# Mosaic-plots reconstruction resuls (abs and phase) for selected slices
function plotReconstruction(reco, slices)

    indexArray = slices

    figure(990)
    absData = 2 .*mapslices(x->abs.(x),reco.data[:,:,slices,1,1],dims=(1,2))
    absMosaic = mosaicview(absData, nrow=Int(floor(sqrt(length(slices)))),npad=5,rowmajor=true, fillvalue=0)

    PyPlot.imshow(absMosaic,cmap="gray", vmax = 0.5e-5, vmin = 0)
    colorbar()

    gcf().suptitle("|Images|")

    figure(1010)

    phaseData = mapslices(x -> ROMEO.unwrap(x),angle.(reco.data[:,:,slices,1,1]),dims=(1,2))
    phaseMosaic = mosaicview(phaseData,nrow=Int(floor(sqrt(length(slices)))),npad=5,rowmajor=true, fillvalue=0)

    PyPlot.imshow(phaseMosaic, cmap="inferno",vmax = 4*pi,vmin=-4*pi)
    colorbar()

    gcf().suptitle("∠Images")

end

## PREPROCESSING

# Function to synchronize trajectory and kspace data as they are not sampled 1:1
function syncTrajAndData!(rawData, traj, idx_crop, interleave)

    # get number of gradient samples
    numGradSamples = traj.numSamplingPerProfile

    # get vector of gradient samples which pertain to one interleave of the trajectory
    ilExtractionVector = numGradSamples*(interleave-1) .+ (1:numGradSamples)

    # Read the trajectory nodes into the rawAcquisitionData type field (.traj)
    for l = 1:length(rawData.profiles)
        rawData.profiles[l].traj = traj.nodes[:,ilExtractionVector]
    end

    # define dwell times for trajectory (dt_k) and signal sampling (dt_s)
    dt_s = 2*10^(-6) # [s]
    dt_k = 10*10^(-6) # [s]

    # Go through every profile (this means every slice in the MRIReco.jl convention for multislice, multiTE, and diffusion scans)
    for l = 1:length(rawData.profiles)

        # Get size of trajectory and signal vectors
        Ns = size(rawData.profiles[l].data,1)
        Nk = size(rawData.profiles[l].traj,2)

        # Define time vectors for signal and trajectory
        t_s = (0:Ns-1)*dt_s
        t_k = (0:Nk-1)*dt_k

        # Interpolate trajectory onto the same sample times as the sampled signal
        trajNodes_interpolated_X = Spline1D(t_k,rawData.profiles[l].traj[1,:], w=ones(length(rawData.profiles[l].traj[1,:])), k=3, bc = "zero")
        trajNodes_interpolated_Y = Spline1D(t_k,rawData.profiles[l].traj[2,:], w=ones(length(rawData.profiles[l].traj[2,:])), k=3, bc = "zero")

        # Concatenate the trajectory node kx and ky positions
        adjustedTraj = vcat(trajNodes_interpolated_X(t_s)',trajNodes_interpolated_Y(t_s)')

        # Crop the data and trajectory to avoid return-to-center of traj, and also set trajectory as upsampled trajectory adjustedTraj
        rawData.profiles[l].traj = adjustedTraj[:,1:idx_crop]
        rawData.profiles[l].data = rawData.profiles[l].data[1:idx_crop,:]

    end

    # Return the vector of sampling times
    return dt_s*(0:idx_crop-1)

end

## Function to correct k0 fluctuations during trajectory
function do_k0_correction!(rawData,k0_phase_modulation, interleave)



    # Get number of samples for the k0 phase modulation (should be same size as the trajectory BEFORE resampling)
    numk0_samples = size(k0_phase_modulation,1)

    # @info size(rawData.profiles) # DEBUG


    # define dwell times for phase modulation (dt_k) and signal sampling (dt_s)
    dt_s = 2*10^(-6) # [s]
    dt_k = 10*10^(-6) # [s]

    # Go through every profile (this means every slice in the MRIReco.jl convention for multislice, multiTE, and diffusion scans)
    for l = 1:length(rawData.profiles)

        # Get size of data (signal samples) and size of k0 modulation
        Ns = size(rawData.profiles[l].data,1)
        Nk = numk0_samples

        # Define time vectors for k0 and signal sampling times
        t_s = (0:Ns-1)*dt_s
        t_k = (0:Nk-1)*dt_k

        # interpolate k0 to the time basis of the signal
        k0_interpolant = Spline1D(t_k,k0_phase_modulation[:,interleave], w=ones(numk0_samples), k=3, bc = "zero")
        k0_interpolated = k0_interpolant(t_s)

        # modulate the data by the k0 modulation by multiplying with e^(i*k0) where k0 is in radians
        rawData.profiles[l].data = rawData.profiles[l].data.* exp.(1im .* k0_interpolated)

        # Visualization of Phase Modulation
        figure("Phase Modulation")
        plot(t_s, angle.(exp.(1im .* k0_interpolated)))
        xlabel("Time [s]")
        ylabel("k₀ [rad]")
        title("B₀ Eddy Current Fluctuation During Readout ")

    end

end

## Function to adjust the header data in the acquisition (requires mutable header)

function adjustHeader!(raw, reconSize, numSamples, interleaveNumber, singleSlice)

    # For every profile in the acquisition
    for l = 1:length(raw.profiles)

        # Set the discard post to 0 (don't discard any samples from the end of the acquisition)
        raw.profiles[l].head.discard_post = 0

        # Set the discard pre to 0 (don't discard any samples from the beginning of the acqusition)
        raw.profiles[l].head.discard_pre = 0

        # Set the contrast to 0 or raw.profiles[l].head.idx.repetition for diffusion directions
        raw.profiles[l].head.idx.contrast = 0 #raw.profiles[l].head.idx.repetition

        # Set the repetition to 0
        raw.profiles[l].head.idx.repetition = 0

        # Set the number of samples properly
        raw.profiles[l].head.number_of_samples = numSamples

        # Set the non-standard encode step (interleave dimension) into the encode step 1 field
        # raw.profiles[l].head.idx.kspace_encode_step_1 = 0
        raw.profiles[l].head.idx.kspace_encode_step_1 = interleaveNumber - 1 # IF MULTI-INTERLEAVE

        # Set slice to 0 for singleslice, if it is not 0 then there will be an error
        if singleSlice
            raw.profiles[l].head.idx.slice = 0
        end

        # Set center sample to 0 (only for spiral scans)

        raw.profiles[l].head.center_sample = 0

        # print("\n",raw.profiles[l].head.idx)

    end

    # Set encoding size to the reconSize
    raw.params["encodedSize"] = [reconSize[1],reconSize[2],1]

end

## PREPROCESSING
function mergeInterleaves(params)

    # constant flag for GIRF correction (maybe set up as a listener later on in the development)
    doGIRF = true

    # Get the other interleave indexes other than the one asked for
    otherInterleaveIndices = [x for x ∈ 1:params[:numInterleaves] if x ∉ params[:interleave]]

    @info "indices = $otherInterleaveIndices"

    # read in the data file from the ISMRMRD format
    dataFile = ISMRMRDFile(params[:interleaveDataFileNames][params[:interleave]])

    # Read in the gradient file and perform GIRF correction to calculate trajectory and calculate k0 phase modulation
    traj, k0_correction = read_gradient_txt_file(params[:trajFilename],params[:reconSize],params[:delay], doGIRF)

    # Read in raw data from the dataFile
    rawData = RawAcquisitionData(dataFile)

    # delete everything that is not a chosen excitation (for efficiency)
    indices = 1:length(rawData.profiles)
    ic = [x for x ∈ indices if x ∉ params[:excitations]]
    deleteat!(rawData.profiles,ic)

    # @info "indices = $ic" #DEBUG

    # set up time vector for tracking all of the interleaves
    timeTrack = []

    # synchronize trajectory data and the kspace data
    times = syncTrajAndData!(rawData, traj, params[:numSamples], params[:interleave])

    # Perform data correction to offset k0 modulation
    do_k0_correction!(rawData, k0_correction, params[:interleave])

    # adjust the header so that each diffusion direction is considered as a contrast instead of a repetition
    adjustHeader!(rawData, params[:reconSize], params[:numSamples], params[:interleave],params[:singleSlice])

    # add the times to the time tracking vector
    append!(timeTrack,times)

    # Repeat the above steps for each interleave, adjusting the times and headers appropriately
    if params[:doMultiInterleave]

        for l in otherInterleaveIndices

            # read in separate interleave data file
            dataFileTemp = ISMRMRDFile(params[:interleaveDataFileNames][l])
            rawDataTemp = RawAcquisitionData(dataFileTemp)
            deleteat!(rawDataTemp.profiles,ic) # delete profiles which aren't needed

            # read in the gradient file and perform the GIRF correction
            trajTemp, k0_correction = read_gradient_txt_file(params[:trajFilename],params[:reconSize],params[:delay], doGIRF)

            # synchronize the trajectory from the gradient file and the data from the raw data file for the interleave
            timesTemp = syncTrajAndData!(rawDataTemp, trajTemp, params[:numSamples], l)

            # perform eddy current correction to the data
            do_k0_correction!(rawDataTemp, k0_correction, l)

            # adjust the header to reflect the arrangement of data expected by MRIReco.jl's reconstruction function
            adjustHeader!(rawDataTemp,params[:reconSize], params[:numSamples], l, params[:singleSlice])

            # append the important data (the profile and the sampling times) to the raw Data file created out of this look
            append!(rawData.profiles,deepcopy(rawDataTemp.profiles))
            append!(timeTrack,deepcopy(timesTemp))

        end

    # if there is the choice to do odd or opposing interleaves, add the 2nd interleave
    elseif params[:doOddInterleave]

        dataFileTemp = ISMRMRDFile(params[:interleaveDataFileNames][3])
        rawDataTemp = RawAcquisitionData(dataFileTemp)
        deleteat!(rawDataTemp.profiles,ic)

        trajTemp, k0_correction = read_gradient_txt_file(params[:trajFilename],params[:reconSize],params[:delay], doGIRF)

        timesTemp = syncTrajAndData!(rawDataTemp, trajTemp, params[:numSamples], 3)

        do_k0_correction!(rawDataTemp, k0_correction, 3)

        adjustHeader!(rawDataTemp,params[:reconSize], params[:numSamples], 2, params[:singleSlice])

        append!(rawData.profiles,copy(rawDataTemp.profiles))
        append!(timeTrack,timesTemp)

    end

    # converting rawData to AcquisitionData
    @info "Converting RawAcquisitionData to AcquisitionData"
    acqData = AcquisitionData(rawData,estimateProfileCenter=false)

    ## Assume all of the slices share a trajectory
    for l = 1:length(acqData.traj)

        acqData.traj[l].times = timeTrack # set times to the total time vector
        acqData.traj[l].TE = 0.00 # set the TE to 0
        acqData.traj[l].AQ = times[end] # set the acquisition time to the last element of the time vector (should be the latest time)
        acqData.traj[l].circular = false # set whether to use a circular filter on the kspace data

    end

    # correctedRaw = RawAcquisitionData(acqData)
    # fout = ISMRMRDFile("data/nonCartesianProcessed.h5")
    # save(fout,correctedRaw)

    # return the acquisition data object with everything corrected
    return acqData

end

## Get gradients from the trajectory
function nodes_to_gradients(nodes::Matrix)

    ## Normalized Conversion (norm kspace to grads in T/m) is scalingFactor = reconSize/(gamma*dwellTime*FOV)

    gradients = diff(hcat([0; 0], nodes), dims = 2)
    gradients = gradients*10.6941176
    return gradients

end

## Convert gradients to trajectory nodes
function gradients_to_nodes(gradients::Matrix)

    ## Normalized Conversion (grads in T/m to normalized k-space) is scalingFactor = (gamma*dwellTime*FOV)/reconSize

    nodes = cumsum(gradients, dims = 2)
    nodes = nodes/10.6941176
    return nodes

end

function checkAcquisitionNodes!(a::AcquisitionData)

    a.traj[1].nodes[abs.(a.traj[1].nodes[:]) .> 0.5] .= 0.5

end

## Function for correcting conversion errors in MRD from Siemens data. Mainly related to the ordering of indices and fov units
function validateSiemensMRD!(r::RawAcquisitionData)

    @info "Validating Siemens converted data"

    ## FOV CHECK:

    if maximum(r.params["encodedFOV"]) > 0.8 # FOV should never be greater than the bore size in [m]

        @info "FOV was recorded in [mm]!Changing to [m]"
        r.params["encodedFOV"] = r.params["encodedFOV"]./1000

    end

end

function validateAcqData!(a::AcquisitionData)

    ## Dimensions CHECK:

    #/ TODO add dimension check that the k-space encoding counters are set properly:
    # kdata dimensions: dim1:=contrast/echo | dim2:=slices | dim3:=repetitions 
    # kdata element dimensions: dim1:=kspace nodes | dim2:=channels/coils

    permutedims(a.kdata,[3,2,1])

end

## Abstraction of the preprocessing of raw MRD data
# OUTPUT IS AN ISMRMRD FILE READY FOR READING WITHOUT FURTHER PROCESSING 
function preprocessCartesianData(r::RawAcquisitionData, fname)

    # Convert rawAcquisitionData object to an AcquisitionData object (these can be reconstructed)
    acqDataCartesian = AcquisitionData(r,estimateProfileCenter=true)

    ## Properly arrange data from the converted siemens file

    validateAcqData!(acqDataCartesian)

    # Fix the FOV (can be set incorrectly)
    #acqDataCartesian.fov[1] = 0.22

    # Need to permute the dimensions of kdata to match the convention of MRIReco.jl
    permutedims(acqDataCartesian.kdata,[3,2,1])

    raw = RawAcquisitionData(acqDataCartesian)

    fout = ISMRMRDFile(fname)

    save(fout, raw)

end

function preprocessNonCartesianData(r::RawAcquisitionData)



end

## PREPROCESSING
function mergeRawInterleaves(params)

    # Get the other interleave indexes other than the one asked for
    otherInterleaveIndices = [x for x ∈ 1:params[:numInterleaves] if x ∉ params[:interleave]]

    @info "indices = $otherInterleaveIndices"

    # read in the data file from the ISMRMRD format
    dataFile = ISMRMRDFile(params[:interleaveDataFileNames][params[:interleave]])

    # Read in the gradient file and perform GIRF correction to calculate trajectory and calculate k0 phase modulation
    traj = read_gradient_txt_file(params[:trajFilename],params[:reconSize],params[:delay])

    # Read in raw data from the dataFile
    rawData = RawAcquisitionData(dataFile)

    # delete everything that is not a chosen excitation (for efficiency)
    indices = 1:length(rawData.profiles)
    ic = [x for x ∈ indices if x ∉ params[:excitations]]
    deleteat!(rawData.profiles,ic)

    # @info "indices = $ic" #DEBUG

    # set up time vector for tracking all of the interleaves
    timeTrack = []

    # synchronize trajectory data and the kspace data
    times = syncTrajAndData!(rawData, traj, params[:numSamples], params[:interleave])

    # adjust the header so that each diffusion direction is considered as a contrast instead of a repetition
    adjustHeader!(rawData, params[:reconSize], params[:numSamples], params[:interleave],params[:singleSlice])

    # add the times to the time tracking vector
    append!(timeTrack,times)

    # Repeat the above steps for each interleave, adjusting the times and headers appropriately
    if params[:doMultiInterleave]

        for l in otherInterleaveIndices

            # read in separate interleave data file
            dataFileTemp = ISMRMRDFile(params[:interleaveDataFileNames][l])
            rawDataTemp = RawAcquisitionData(dataFileTemp)
            deleteat!(rawDataTemp.profiles,ic) # delete profiles which aren't needed

            # read in the gradient file and perform the GIRF correction
            trajTemp = read_gradient_txt_file(params[:trajFilename],params[:reconSize],params[:delay])

            # synchronize the trajectory from the gradient file and the data from the raw data file for the interleave
            timesTemp = syncTrajAndData!(rawDataTemp, trajTemp, params[:numSamples], l)

            # adjust the header to reflect the arrangement of data expected by MRIReco.jl's reconstruction function
            adjustHeader!(rawDataTemp,params[:reconSize], params[:numSamples], l, params[:singleSlice])

            # append the important data (the profile and the sampling times) to the raw Data file created out of this look
            append!(rawData.profiles,deepcopy(rawDataTemp.profiles))
            append!(timeTrack,deepcopy(timesTemp))

        end

    # if there is the choice to do odd or opposing interleaves, add the 2nd interleave
    elseif params[:doOddInterleave]

        dataFileTemp = ISMRMRDFile(params[:interleaveDataFileNames][3])
        rawDataTemp = RawAcquisitionData(dataFileTemp)
        deleteat!(rawDataTemp.profiles,ic)

        trajTemp = read_gradient_txt_file(params[:trajFilename],params[:reconSize],params[:delay])

        timesTemp = syncTrajAndData!(rawDataTemp, trajTemp, params[:numSamples], 3)

        adjustHeader!(rawDataTemp,params[:reconSize], params[:numSamples], 2, params[:singleSlice])

        append!(rawData.profiles,copy(rawDataTemp.profiles))
        append!(timeTrack,timesTemp)

    end

    # converting rawData to AcquisitionData
    @info "Converting RawAcquisitionData to AcquisitionData"
    acqData = AcquisitionData(rawData,estimateProfileCenter=true)

    ## Assume all of the slices share a trajectory
    for l = 1:length(acqData.traj)

        acqData.traj[l].times = timeTrack # set times to the total time vector
        acqData.traj[l].TE = 0.00 # set the TE to 0
        acqData.traj[l].AQ = times[end] # set the acquisition time to the last element of the time vector (should be the latest time)
        acqData.traj[l].circular = false # set whether to use a circular filter on the kspace data

    end

    for l = 1:length(acqData.subsampleIndices)

        acqData.subsampleIndices[l] = 1:size(acqData.traj[l].nodes,2)

    end


    # return the acquisition data object with everything corrected
    return acqData

end

## Currently only works for spiral trajectories but will have to be extended to cartesian!
function parseTrajectoryGradients(a::AcquisitionData)

    for l = 1:length(a.traj)
        
        nProfiles = a.traj[l].numProfiles
        nSamples = a.traj[l].numSamplingPerProfile
        nodes = a.traj[l].nodes

        for profile = 1:nProfiles
            
            ilExtractor = nSamples*(profile-1) .+ (1:nSamples)
            ilNodes = nodes[:,ilExtractor]

            figure()
            ilGrads = nodes_to_gradients(ilNodes)

            plot(ilGrads')

        end
    
    end
        
end

## ApplyGIRF to AcqData Object
function applyGIRF!(a::AcquisitionData, freq::AbstractVector, g_data::AbstractMatrix)

    for l = 1:length(a.traj)
        
        nProfiles = a.traj[l].numProfiles
        nSamples = a.traj[l].numSamplingPerProfile
        nodes = a.traj[l].nodes
        times = a.traj[l].times
        oldNodes = a.traj[l].nodes

        for profile = 1:nProfiles
            
            ilExtractor = nSamples*(profile-1) .+ (1:nSamples)
            ilNodes = nodes[:,ilExtractor]
            ilTimes = times[ilExtractor]

            ilGrads = nodes_to_gradients(ilNodes)

            for dim = 1:size(ilGrads,1)

                correctedGrads = predictGrad_port(freq,g_data[:,dim],ilTimes,ilGrads[dim,:], ilTimes) # THESE ARE ALL VECTORS SO INPUT orientation (column/row major ordering) doesn't matter
                ilGrads[dim,:] = correctedGrads'

            end

            ilNodes = gradients_to_nodes(ilGrads)
            nodes[:,ilExtractor] = ilNodes

        end

        a.traj[l].nodes = nodes
    
    end

end

## ApplyGIRF to AcqData Object
function applyK0!(a::AcquisitionData, freq::AbstractVector, k0_data::AbstractMatrix)

    for l = 1:length(a.traj)
        
        nProfiles = a.traj[l].numProfiles
        nSamples = a.traj[l].numSamplingPerProfile
        nodes = a.traj[l].nodes
        times = a.traj[l].times
        oldNodes = a.traj[l].nodes

        for profile = 1:nProfiles
            
            ilExtractor = nSamples*(profile-1) .+ (1:nSamples)
            ilNodes = nodes[:,ilExtractor]
            ilTimes = times[ilExtractor]

            ilGrads = nodes_to_gradients(ilNodes)

            k0_correction = ones(size(ilGrads))
            @info size(k0_correction)

            for dim = 1:size(ilGrads,1)

                k0_correction[dim,:] = predictGrad_port(freq,k0_data[:,dim],ilTimes,ilGrads[dim,:], ilTimes) # THESE ARE ALL VECTORS SO INPUT orientation (column/row major ordering) doesn't matter

            end

            finalCorrection = sum(k0_correction,dims=1) #back to radians!

            a.kdata[l][ilExtractor,:] = a.kdata[l][ilExtractor,:] .* exp.(1im .*finalCorrection')

            # Visualization of Phase Modulation
            figure("Phase Modulation 2")
            plot(vec(ilTimes), vec(angle.(exp.(1im .* finalCorrection))))
            xlabel("Time [s]")
            ylabel("k₀ [rad]")
            title("B₀ Eddy Current Fluctuation During Readout ")
            
        end
    
    end

end

## Calibrate the phase from individual interleaves
function calibrateAcquisitionPhase!(a::AcquisitionData)

    for l = 1:length(a.traj)
        
        nProfiles = a.traj[l].numProfiles
        nSamples = a.traj[l].numSamplingPerProfile
        nodes = a.traj[l].nodes
        times = a.traj[l].times
        oldNodes = a.traj[l].nodes

        initialInterleavePhase = angle.(a.kdata[l][1,:])'

        for profile = 2:nProfiles
            
            ilExtractor = nSamples*(profile-1) .+ (1:nSamples)

            initialProfilePhase = angle.(a.kdata[l][ilExtractor[1],:])'

            @info size(initialInterleavePhase)
            @info size(a.kdata[l][ilExtractor,:])
            
            a.kdata[l][ilExtractor,:] .*= exp.(-1im * (initialInterleavePhase - initialProfilePhase))
            
        end    
    end

end

function testConversion()

    gamma = 42577478
    dt = 2e-6

    # correction to k-space = kspace * fov (m) / reconSize

    0.5*200 / 0.22

end