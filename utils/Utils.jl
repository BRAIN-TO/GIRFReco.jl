export plotReconstruction, plotSenseMaps, calculateB0Maps, getSliceOrder, syncTrajAndData!, do_k0_correction!, adjustHeader!, checkAcquisitionNodes!, validateSiemensMRD!, validateAcqData!, preprocessCartesianData, removeOversampling!, mergeRawInterleaves, applyGIRF!, applyK0!, saveMap, loadMap, shiftksp!

## Choose plotting backend to be PlotlyJS!
# plotlyjs()

## General Plotting function for the reconstruction

# "Mosaic-plots reconstruction for selected slices and corresponding B0 map"
"""
    plotReconstruction(images, slicesIndex, b0; figHandles = [], isSliceInterleaved = false)
Plots the magnitude and phase of the reconstructed images for a given slice or slices, along with a B₀ map if applicable

# Arguments
* `images` - Complex-valued images reconstructed using MRIReco.jl
* `slicesIndex::Vector{Int}` - slices to plot
* `b0` - off-resonance map to plot along with images
* `figHandles` - String vectors in size of [3,1] for titles of three figures (Magnitude & Phase of reconstructed images, and B0 maps)
* `isSliceInterleaved::Bool` - for 2D scanning, indicate this value as `true` to make sure the slice order on the displayed results is correct
* `rotateAngle::Int` - Counterclock-wise rotation angle for each slice, should be a value from 0, 90, 180, 270 degrees
"""
function plotReconstruction(images, slicesIndex, b0; figHandles = [], isSliceInterleaved = false, rotateAngle = 0)
    # plot()
    ## If we need to re-order all slices
    sliceNum = length(slicesIndex)
    reorderSliceIndex = zeros(Int16, size(slicesIndex))

    if isSliceInterleaved && sliceNum > 1
        reorderSliceIndex[1:2:end] = slicesIndex[1:Int(ceil(sliceNum / 2))]
        reorderSliceIndex[2:2:end] = slicesIndex[Int(ceil(sliceNum / 2) + 1) : end]
    else
        reorderSliceIndex = slicesIndex
    end

    ## If we need to rotate each slice
    if mod(rotateAngle, 90) != 0 || rotateAngle < 0 || rotateAngle > 270
        error("rotateAngle must be 0, 90, 180 or 270 degrees.")
    end

    # Plot magnitude images (normalize)
    # if length(figHandles) < 1
    #     figure("Magnitude Images")
    # else
    #     figure(figHandles[1])
    # end

    # clf()

    absData = mapslices(x -> abs.(x) ./ maximum(abs.(x)), images[:, :, reorderSliceIndex], dims = [1,2])
    if rotateAngle == 90
        absData = mapslices(x -> rotr90(x), absData, dims = [1,2])
    elseif rotateAngle == 180
        absData = mapslices(x -> rot180(x), absData, dims = [1,2])
    else
        absData = mapslices(x -> rotl90(x), absData, dims = [1,2])
    end
    absMosaic = mosaicview(absData, nrow = Int(floor(sqrt(sliceNum))), npad = 5, rowmajor = true, fillvalue = 0)

    heatmap(absMosaic,show=true, plot_title="|Images|",plot_titlevspan=0.1,color=:grays,aspectratio=:equal)
    # display(plot!())

    phaseData = angle.(images[:, :, reorderSliceIndex, 1, 1])
    if rotateAngle == 90
        phaseData = mapslices(x -> rotr90(x), phaseData, dims = [1,2])
    elseif rotateAngle == 180
        phaseData = mapslices(x -> rot180(x), phaseData, dims = [1,2])
    else
        phaseData = mapslices(x -> rotl90(x), phaseData, dims = [1,2])
    end
    phaseMosaic = mosaicview(phaseData, nrow = Int(floor(sqrt(sliceNum))), npad = 5, rowmajor = true, fillvalue = 0)

    heatmap(phaseMosaic,show=true,plot_title="∠ Images",plot_titlevspan=0.1,color=:plasma,aspectratio=:equal)
    #colorbar()
    # display(plot!())

    #gcf().suptitle("∠Images")

    # # Plot B0 maps
    # if length(figHandles) < 3
    #     figure("B₀ Map Images")
    # else
    #     figure(figHandles[3])
    # end

    # clf()
    # plot()
    b0Data = mapslices(x -> x, b0, dims = [1,2])
    if rotateAngle == 90
        b0Data = mapslices(x -> rotr90(x), b0, dims = [1,2])
    elseif rotateAngle == 180
        b0Data = mapslices(x -> rot180(x), b0, dims = [1,2])
    else
        b0Data = mapslices(x -> rotl90(x), b0, dims = [1,2])
    end
    b0Mosaic = mosaicview(b0Data[:, :, reorderSliceIndex], nrow = Int(floor(sqrt(sliceNum))), npad = 5, rowmajor = true, fillvalue = 0)

    heatmap(b0Mosaic,show=true, plot_title="B₀ Map Images",plot_titlevspan=0.1,color=:plasma)
    # colorbar()
    # display(plot!())

    # gcf().suptitle("B₀ Maps [rad/s]")
    
    1

end

"Function plots all profiles in the acquisition to check consistency with ISMRMRD file"
function checkProfiles(rawData)

    numProfiles2 = 128 # Set to the number of profiles that you would like to see

    for l = 1:numProfiles2
        p1 = plot(abs.(rawData.profiles[l].data[:,1]))
        p2 = plot(angle.(rawData.profiles[l].data[:,1]))
    end

end

"""
    plotSenseMaps!(sense, n_channels)
Plots coil sensitivity maps from the channels, for a total of n_channels plots

# Arguments
* `sense` - sensitivity maps
* `n_channels::Int` - number of coils (usually last dimension of sense)
* `sliceIndex::Int` - The index of the slice to be displayed (if multislice)
"""
function plotSenseMaps(sense,n_channels; sliceIndex = 1)
    sliceNum = size(sense, 3)
    if sliceIndex > sliceNum
        errMsg = @sprintf("The index of slice to be displayed is %d, but total slice number is %d.", sliceIndex, sliceNum)
        error(errMsg)
    end

    # # Magnitude maps
    # figure(@sprintf("Sensitivity Map Magnitude of Slice %d / %d", sliceIndex, sliceNum)); clf(); for ch in 1:n_channels; subplot(8,4,ch); imshow((abs.(sense[:,:,sliceIndex,ch])), cmap = "gray"); end;
    # subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
    # gcf()

    magMosaic = mosaicview((abs.(sense[:,:,sliceIndex,:])), nrow = Int(floor(sqrt(n_channels))), npad = 5, rowmajor = true, fillvalue = 0)
    heatmap(magMosaic, show=true, plot_title="|Sensitivity|",plot_titlevspan=0.1,color=:gnuplot2)

    # # Phase maps
    # figure(@sprintf("Sensitivity Map Phase of Slice %d / %d", sliceIndex, sliceNum)); clf(); for ch in 1:n_channels; subplot(8,4,ch); imshow(angle.(sense[:,:,sliceIndex,ch]), cmap = "gray"); end;
    # subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
    # gcf()

    phaseMosaic = mosaicview((angle.(sense[:,:,sliceIndex,:])), nrow = Int(floor(sqrt(n_channels))), npad = 5, rowmajor = true, fillvalue = 0)
    heatmap(phaseMosaic,show=true, plot_title="∠ Sensitivity",plot_titlevspan=0.1,color=:plasma)

end

# "WIP: Plots trajectory and Data, doesn't work currently"
function plotTrajAndData(acq)

    for l in 1:length(acq.traj)

        freqEncode[l,:] = acq.traj[l].nodes[1,:]
        phaseEncode[l,:] = acq.traj[l].nodes[2,:]
        kSignal[l,:] = acq.kdata[l,:,1]

    end

end

## PREPROCESSING

"""
    calculateB0Maps(imData,slices,echoTime1,echoTime2)

Calculate  B0 map from the two images with different echo times via their phase difference (obtained from imTE2.*conj(imTE1))
TODO have the b0 map calculation be capable of handling variable echo times
TODO2: Do we need this basic B0 map calculation or is it superseded by estimateB0Maps?

# Arguments
* `imdata`                          - [nX nY nZ 2 nCoils] 5D image array, 4th dim echo time
* `slices::NTuple{nSlices,Int}`     - slice index vector (tuple?) for which map is computed
* `echoTime1::AbstractFloat`        - TE1 [ms]
* `echoTime2::AbstractFloat`        - TE2 [ms]
"""
function calculateB0Maps(imData,slices,echoTime1,echoTime2)

    # b0Maps = mapslices(x -> rotl90(x),ROMEO.unwrap(angle.(imData[:,:,slices,2,1].*conj(imData[:,:,slices,1,1]))),dims=(1,2))./((7.38-4.92)/1000)
    b0Maps = mapslices(x -> x, ROMEO.unwrap(angle.(imData[:,:,slices,2,1].*conj(imData[:,:,slices,1,1]))),dims=(1,2))./((echoTime2-echoTime1)/1000)

end

"""
    getSliceOrder(nSlices, isSliceInterleaved)

Returns array mapping from acquisition number to slice number (geometric position) (indexArray[slice = 1:9] = [acquisitionNumbers])
TODO: Add ascending/descending options

# Arguments
* `nSlices::Int`                    - number of slices in total acquired stack (FOV)
* `isSliceInterleaved::Bool=true`   - if true, interleaved slice order is created, otherwise ascending slice order is returned
"""
function getSliceOrder(nSlices; isSliceInterleaved::Bool=true)
    
    sliceIndexArray = 1:nSlices
    reorderedSliceIndexArray = zeros(Int16, size(sliceIndexArray))
    if isSliceInterleaved && nSlices > 1
        reorderedSliceIndexArray[1:2:end] = sliceIndexArray[1:Int(ceil(nSlices / 2))]
        reorderedSliceIndexArray[2:2:end] = sliceIndexArray[Int(ceil(nSlices / 2) + 1) : end]
    else
        reorderedSliceIndexArray = sliceIndexArray
    end

    return reorderedSliceIndexArray

end

"""
    syncTrajAndData!(a::AcquisitionData)
Synchronizes k-space trajectory and sampled data as they do not usually have a common sampling rate

# Arguments
* `rawData::RawAcquisitionData` - RawAcquisitionData object
* `traj::Trajectory` - Trajectory object to be synchronized with data contained in rawData
* `idx_crop::Int` - Trajectory and Data may contain samples we don't want in the recon, usually at the end of acquisition. Ignore samples after idx_crop
* `interleave::Int` - index of interleave
"""
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


"""
    do_k0_correction!(rawData, k0_phase_modulation, interleave)
Applies phase modulation due to 0th-order field fluctuations during the acquisition

# Arguments
* `rawData::RawAcquisitionData` - RawAcquisitionData object
* `k0_phase_modulation::Matrix{Complex{T}}` - Vector containing phase modulation measurements
* `interleave::Int` - index of interleave
"""
function do_k0_correction!(rawData,k0_phase_modulation, interleave)

    # Get number of samples for the k0 phase modulation (should be same size as the trajectory BEFORE resampling)
    numk0_samples = size(k0_phase_modulation,1)

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

        # # Visualization of Phase Modulation
        # figure("Phase Modulation")
        # plot(t_s, angle.(exp.(1im .* k0_interpolated)))
        # xlabel("Time [s]")
        # ylabel("k₀ [rad]")
        # title("B₀ Eddy Current Fluctuation During Readout ")

        plot(t_s, angle.(exp.(1im .* k0_interpolated)),show=true,title="B₀ Eddy Current Fluctuation During Readout ")

    end

end


"""
    adjustHeader!(raw::RawAcquisitionData, reconSize, numSamples, interleaveNumber, singleSlice)
Adjusts the header data for each interleave and slice of spiral diffusion RawAcquisitionData

# Arguments
* `raw::RawAcquisitionData` - RawAcquisitionData object
* `reconSize::Vector` - Reconstruction matrix size
* `numSamples::Int` - Number of samples per interleave
* `interleaveNumber::Int` - Index of interleave for multi-shot acquisitionNumbers
* `singleSlice::Bool` - flag for single-slice reconstruction/acquisition
"""
function adjustHeader!(raw, reconSize, numSamples, interleaveNumber, singleSlice)

    # For every profile in the acquisition
    for l = 1:length(raw.profiles)

        # Set the discard post to 0 (don't discard any samples from the end of the acquisition)
        raw.profiles[l].head.discard_post = 0

        # Set the discard pre to 0 (don't discard any samples from the beginning of the acqusition)
        raw.profiles[l].head.discard_pre = 0

        # Set the contrast to 0 or raw.profiles[l].head.idx.repetition for diffusion directions
        # raw.profiles[l].head.idx.contrast = raw.profiles[l].head.idx.repetition
        raw.profiles[l].head.idx.contrast = 0

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

    end

    # Set encoding size to the reconSize
    raw.params["encodedSize"] = [reconSize[1],reconSize[2],1]

end

"""
    checkAcquisitionNodes!(a::AcquisitionData)
Validates processed AcquisitionData object to make sure that |kᵢ| < 0.5 ∀ i ∈ [1, Nₛ]

# Arguments
* `a::AcquisitionData` - AcquisitionData object
"""
function checkAcquisitionNodes!(a::AcquisitionData)

    a.traj[1].nodes[abs.(a.traj[1].nodes[:]) .> 0.5] .= 0.5

end


"""
    validateSiemensMRD!(r::RawAcquisitionData)
Validates RawAcquisitionData object created from ISMRMRD format object

# Arguments
* `r::RawAcquisitionData` - RawAcquisitionData object
"""
function validateSiemensMRD!(r::RawAcquisitionData)

    @info "Validating Siemens converted data"

    ## FOV CHECK:

    if maximum(r.params["encodedFOV"]) > 0.8 # FOV should never be greater than the bore size in [m]

        @info "FOV was recorded in [mm]!Changing to [m]"
        r.params["encodedFOV"] = r.params["encodedFOV"]./1000

    end

end

"""
    validateAcqData!(a::AcquisitionData)
Validates processed AcquisitionData object after manipulation, etc...

# Arguments
* `a::AcquisitionData` - AcquisitionData object
"""
function validateAcqData!(a::AcquisitionData)

    ## Dimensions CHECK:

    # TODO add dimension check that the k-space encoding counters are set properly:
    # kdata dimensions: dim1:=contrast/echo | dim2:=slices | dim3:=repetitions 
    # kdata element dimensions: dim1:=kspace nodes | dim2:=channels/coils

    permutedims(a.kdata,[3,2,1])
    checkAcquisitionNodes!(a)

end

"""
    preprocessCartesianData!(raw::RawAcquisitionData; dims = 1)
Prepares Cartesian for reconstruction

# Arguments
* `r::RawAcquisitionData{T}`          - RawAcquisitionData object
* `fname`                             - filename to save the preprocessed data to
"""
function preprocessCartesianData(r::RawAcquisitionData, doSave; fname = "data/testFile.h5")

    removeOversampling!(r)

    # Convert rawAcquisitionData object to an AcquisitionData object (these can be reconstructed)
    acqDataCartesian = AcquisitionData(r,estimateProfileCenter=true)

    ## Properly arrange data from the converted siemens file
    validateAcqData!(acqDataCartesian)

    if doSave

        raw = RawAcquisitionData(acqDataCartesian)

        # Since the data should generally have 3D information when saved, we make sure 2D data is appropriately stored as 3D data with a singleton dimension
        if length(raw.params["encodedSize"]) == 2
            e_sz = raw.params["encodedSize"]
            raw.params["encodedSize"] = [e_sz[1], e_sz[2], 1]
        end
        
        # raw.params = headerCopy
        fout = ISMRMRDFile(fname)
        save(fout, raw)

    end

    return acqDataCartesian

end

"""
    removeOversampling!(raw::RawAcquisitionData; dims = 1)
Removes 2x readout oversampling in specified raw data dimensions by iFFT, cropping FOV and FFT

# Arguments
* `raw::RawAcquisitionData{T}`          - RawAcquisitionData object
* `dims`                                - dimension alongside which oversampling is removed (default: 1)
"""
function removeOversampling!(raw::RawAcquisitionData; dims = [1])

    idxDim = dims[1]
    Ns = raw.params["encodedSize"][idxDim]
    idxCropFov = convert(Vector{Int32}, [1:floor(Ns/4); ceil(3/4*Ns+1):Ns])
    
    # For every profile in the acquisition
    for iProfile = 1:length(raw.profiles)
    
        # IFFT to image space, crop, FFT back to k-space
        ifft!(raw.profiles[iProfile].data, idxDim)
        raw.profiles[iProfile].data = fft!(raw.profiles[iProfile].data[idxCropFov,:], idxDim)

    end

    # halve encoding size of first dimension
    raw.params["encodedSize"][idxDim] /= 2
    raw.params["encodedFOV"][idxDim] /= 2

end


"""
    mergeRawInterleaves(params)
Merges multiple interleave data together from individually acquired interleave scans

# Arguments
* `params`          - Dictionary
"""
function mergeRawInterleaves(params)

    # Get the other interleave indexes other than the one asked for
    otherInterleaveIndices = [x for x ∈ 1:params[:numInterleaves] if x ∉ params[:interleave]]

    # @info "indices = $otherInterleaveIndices" #DEBUG

    # read in the data file from the ISMRMRD format
    dataFile = ISMRMRDFile(params[:interleaveDataFileNames][params[:interleave]])

    # Read in the gradient file
    trajAll = read_gradient_txt_file(params[:trajFilename],params[:reconSize],params[:delay])

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
    times = syncTrajAndData!(rawData, trajAll, params[:numSamples], params[:interleave])

    # adjust the header so that each diffusion direction is considered as a contrast instead of a repetition
    # adjustHeader!(rawData, params[:reconSize], params[:numSamples], params[:interleave],params[:singleSlice])
    adjustHeader!(rawData, params[:reconSize], params[:numSamples], 1, params[:singleSlice])

    # add the times to the time tracking vector
    append!(timeTrack,times)

    # Repeat the above steps for each interleave, adjusting the times and headers appropriately
    if params[:doMultiInterleave]

        for l in otherInterleaveIndices

            # read in separate interleave data file
            dataFileTemp = ISMRMRDFile(params[:interleaveDataFileNames][l])
            rawDataTemp = RawAcquisitionData(dataFileTemp)
            deleteat!(rawDataTemp.profiles,ic) # delete profiles which aren't needed

            # synchronize the trajectory from the gradient file and the data from the raw data file for the interleave
            timesTemp = syncTrajAndData!(rawDataTemp, trajAll, params[:numSamples], l)

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

        timesTemp = syncTrajAndData!(rawDataTemp, trajAll, params[:numSamples], 3)

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
        acqData.traj[l].circular = true # set whether to use a circular filter on the kspace data

    end

    for l = 1:length(acqData.subsampleIndices)

        acqData.subsampleIndices[l] = 1:size(acqData.traj[l].nodes,2)

    end

    # return the acquisition data object with everything corrected
    return acqData

end

"""
    applyGIRF!(raw::RawAcquisitionData, freq::AbstractVector, g_data::AbstractMatrix)
Applies the GIRF to the trajectories inside of a::AcquisitionData

# Arguments
* `a::AcquisitionData{T}`          - AcquisitionData object
* `freq::AbstractVector`           - Vector containing frequencies of GIRF data
* `g_data::AbstractMatrix`         - Matrix of size N x length(freq) containing complex GIRF data
"""
function applyGIRF!(a::AcquisitionData{T}, g::GirfApplier) where T

    # Read parameters for gradient and node conversion
    S = a.encodingSize
    F = a.fov

    # Check dimensions of the acquisition data and ensure encoding size and FOV are consistent
    if length(S) == 2
        S = (S[1], S[2], 1)
    end

    if length(F) == 2
        F = Float32.(F[1], F[2], 1.0)
    end

    # loop over all contained trajectories
    for l = 1:length(a.traj)
        
        nProfiles = a.traj[l].numProfiles
        nSamples = a.traj[l].numSamplingPerProfile
        nodes = a.traj[l].nodes
        times = a.traj[l].times
        oldNodes = a.traj[l].nodes

        # loop over all profiles in a trajectory
        for profile = 1:nProfiles
            
            ilExtractor = nSamples*(profile-1) .+ (1:nSamples)
            ilNodes = nodes[:,ilExtractor]
            ilTimes = times[ilExtractor]

            DT = ilTimes[1] - ilTimes[2]

            ilGrads = nodes_to_gradients(ilNodes; dwellTime=DT, reconSize=S, FOV = F)

            # loop over trajectory dimensions
            for dim = 1:size(ilGrads,1)

                correctedGrads = apply_girf(g,ilGrads[dim,:], ilTimes ,ilTimes, dim) # THESE ARE ALL VECTORS SO INPUT orientation (column/row major ordering) doesn't matter
                ilGrads[dim,:] = correctedGrads'

            end

            ilNodes = gradients_to_nodes(ilGrads; dwellTime=DT, reconSize=S, FOV = F)
            nodes[:,ilExtractor] = ilNodes

        end

        a.traj[l].nodes = nodes
    
    end

end

"""
    applyK0!(raw::RawAcquisitionData, freq::AbstractVector, g_data::AbstractMatrix)
Applies the K0 modulation due to imaging gradients to the data inside of a::AcquisitionData

# Arguments
* `a::AcquisitionData{T}`          - AcquisitionData object
* `freq::AbstractVector`           - Vector containing frequencies of GIRF data
* `k0_data::AbstractMatrix`         - Matrix of size N x length(freq) containing complex k0 function data
"""
function applyK0!(a::AcquisitionData{T},g::GirfApplier) where T

    # Read parameters for gradient and node conversion
    S = a.encodingSize
    F = a.fov

    if length(S) == 2
        S = (S[1], S[2], 1)
    end
    if length(F) == 2
        F = (F[1], F[2], 1.0)
    end

    # loop over all contained trajectories
    for l = 1:length(a.traj)
        
        nProfiles = a.traj[l].numProfiles
        nSamples = a.traj[l].numSamplingPerProfile
        nodes = a.traj[l].nodes
        times = a.traj[l].times
        oldNodes = a.traj[l].nodes

        # loop over all profiles in a trajectory
        for profile = 1:nProfiles
            
            ilExtractor = nSamples*(profile-1) .+ (1:nSamples)
            ilNodes = nodes[:,ilExtractor]
            ilTimes = times[ilExtractor]

            DT = ilTimes[1] - ilTimes[2]

            ilGrads = nodes_to_gradients(ilNodes; dwellTime=DT, reconSize=S, FOV = F)

            k0_correction = ones(size(ilGrads))
            
            # loop over all trajectory dims
            for dim = 1:size(ilGrads,1)

                k0_correction[dim,:] = apply_girf(g,ilGrads[dim,:], ilTimes, ilTimes, dim) # THESE ARE ALL VECTORS SO INPUT orientation (column/row major ordering) doesn't matter

            end

            finalCorrection = sum(k0_correction,dims=1) #back to radians!

            a.kdata[l][ilExtractor,:] = a.kdata[l][ilExtractor,:] .* exp.(-1im .*finalCorrection')

            # # Visualization of Phase Modulation
            # figure("Phase Modulation 2")
            # plot(vec(ilTimes), vec(angle.(exp.(1im .* finalCorrection))))
            # xlabel("Time [s]")
            # ylabel("k₀ [rad]")
            # title("B₀ Eddy Current Fluctuation During Readout ")

            # plot(ilTimes, angle.(exp.(1im .* finalCorrection')),show=true,title="B₀ Eddy Current Fluctuation During Readout ") #DEBUG
            
        end
    
    end

end

# ## Calibrate the phase from individual interleaves
# function calibrateAcquisitionPhase!(a::AcquisitionData)

#     for l = 1:length(a.traj)
        
#         nProfiles = a.traj[l].numProfiles
#         nSamples = a.traj[l].numSamplingPerProfile
#         nodes = a.traj[l].nodes
#         times = a.traj[l].times
#         oldNodes = a.traj[l].nodes

#         initialInterleavePhase = angle.(a.kdata[l][1,:])'

#         for profile = 2:nProfiles
            
#             ilExtractor = nSamples*(profile-1) .+ (1:nSamples)

#             initialProfilePhase = angle.(a.kdata[l][ilExtractor[1],:])'

#             @info size(initialInterleavePhase)
#             @info size(a.kdata[l][ilExtractor,:])
            
#             a.kdata[l][ilExtractor,:] .*= exp.(-1im * (initialInterleavePhase - initialProfilePhase))
            
#         end    
#     end

# end

## Input/Output, File handling

"""
    saveMap(filename, calib_map, resolution_mm; offset_mm = [0.0, 0.0, 0.0])
Saves calibration maps (sensitivity or B0) as 4D NIfTI file(s)

For complex-valued data, magnitude and phase can be split into separate files
# Arguments
* `filename::String`            - string filename with extension .nii, example "sensemap.nii"
* `calib_map`                   - [nX nY nZ {nChannels}] 4-D sensitivity or 3D B0 map array 
* `resolution_mm`               - resolution in mm, 3 element vector, e.g., [1.0, 1.0, 2.0]
* `offset_mm`                   - isocenter offset in mm, default: [0.0, 0.0, 0.0]
* `doSplitPhase::Bool=false`    - if true, data is saved in two nifti files with suffix "_magn" and "_phase", respectively
                                  to enable display in typical NIfTI viewers
"""
function saveMap(filename, calib_map, resolution_mm; offset_mm = [0.0, 0.0, 0.0], doSplitPhase::Bool=false, doNormalize::Bool=true)
    
    # multiplication with 1000 should no longer be necessary after MRIReco 0.7.1
    spacing = 1000.0 .*resolution_mm .*Unitful.mm
    offset = 1000.0 .*offset_mm .*Unitful.mm

    if ndims(calib_map) >= 4 # multi-coil calib_map, e.g., sensitivity, or recon, but we can only store the first 4 dims in a Nifti
        I = reshape(calib_map, size(calib_map,1), size(calib_map,2), size(calib_map,3), size(calib_map,4), 1, 1);
    else
        I = reshape(calib_map, size(calib_map,1), size(calib_map,2), size(calib_map,3), 1, 1, 1);
    end

    # scale to max 1
    if doNormalize
        I /= maximum(abs.(I))
    end

    # AxisArray Constructor
    im = AxisArray(I,
    Axis{:x}(range(offset[1], step=spacing[1], length=size(I, 1))),
    Axis{:y}(range(offset[2], step=spacing[2], length=size(I, 2))),
    Axis{:z}(range(offset[3], step=spacing[3], length=size(I, 3))),
    Axis{:coils}(1:size(I, 4)),
    Axis{:echos}(1:size(I, 5)),
    Axis{:repetitions}(1:size(I, 6)))

    # if separate mag and phase are desired, save them separately
    if doSplitPhase

        filename_magn = splitext(filename)[1] * "_magn.nii"
        saveImage(filename_magn, map(abs,im)) # map is needed, because abs.(im) would convert AxisArray back into basic array

        filename_phase = splitext(filename)[1] * "_phase.nii"
        saveImage(filename_phase, map(angle,im))

    else

        saveImage(filename, im)

    end

end

"""
    loadMap(filename, calib_map, resolution_mm; offset_mm = [0.0, 0.0, 0.0])
Saves calibration maps (sensitivity or B0) as 4D NIfTI file(s)

For complex-valued data, magnitude and phase can be split into separate files
# Arguments
* `filename::String`            - string filename with extension .nii, example "sensemap.nii"
* `doSplitPhase::Bool=false`    - if true, data is saved in two nifti files with suffix "_magn" and "_phase", respectively
                                  to enable display in typical NIfTI viewers
# Output
* `calib_map`                    - [nX nY nZ {nChannels}] 4-D sensitivity or 3D B0 map array 
"""
function loadMap(filename; doSplitPhase::Bool=false)
    
    # if separate mag and phase are saved, load and combine them
    if doSplitPhase

        filename_magn = splitext(filename)[1] * "_magn.nii"
        I_magn = loadImage(filename_magn) # map is needed, because abs.(im) would convert AxisArray back into basic array

        filename_phase = splitext(filename)[1] * "_phase.nii"
        I_phase = loadImage(filename_phase)

        calib_map = (I_magn.data).*exp.(1im.*(I_phase.data))

    else

        I = loadImage(filename)
        calib_map = I.data

    end
    
    # squeeze singleton dimensions of 6-dim array
    calib_map = dropdims(calib_map, dims = tuple(findall(size(calib_map) .== 1)...))

    return calib_map

end

function shiftksp!(acqData,shift)

    numSl = numSlices(acqData)
    numRep, numContr = numRepetitions(acqData), numContrasts(acqData)

    smat = prod(exp.(1im .* acqData.traj[1].nodes[:,acqData.subsampleIndices[1]] .* shift .* 2 .* pi),dims=1)

    for slice = 1:numSl
        for contr = 1:numContr
            for rep = 1:numRep                    
                acqData.kdata[contr,slice,rep] = acqData.kdata[contr,slice,rep] .* smat'                    
            end
        end
    end

end