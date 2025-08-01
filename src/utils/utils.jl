"""
    plot_reconstruction(images, slices_index, b0; is_slice_interleaved = false, rotation = 0)
Plots the magnitude and phase of the reconstructed images for a given slice or slices, along with a B₀ map if applicable

# Arguments
* `images` - Complex-valued images reconstructed using MRIReco.jl
* `slices_index::Vector{Int}` - slices to plot
* `b0` - off-resonance map to plot along with images
* `is_slice_interleaved::Bool` - for 2D scanning, indicate this value as `true` to make sure the slice order on the displayed results is correct
* `rotation::Int` - Counterclock-wise rotation angle for each slice, should be a value from 0, 90, 180, 270 degrees
"""
function plot_reconstruction(images, slices_index, b0; is_slice_interleaved = false, rotation = 0)

    num_slices = length(slices_index)
    reordered_slice_indices = zeros(Int16, size(slices_index))

    if is_slice_interleaved && num_slices > 1
        reordered_slice_indices[1:2:end] = slices_index[1:Int(ceil(num_slices / 2))]
        reordered_slice_indices[2:2:end] = slices_index[Int(ceil(num_slices / 2) + 1):end]
    else
        reordered_slice_indices = slices_index
    end

    ## If we need to rotate each slice
    if mod(rotation, 90) != 0 || rotation < 0 || rotation > 270
        error("rotation must be 0, 90, 180 or 270 degrees.")
    end

    magnitude_images = mapslices(x -> abs.(x) ./ maximum(abs.(x)), images[:, :, reordered_slice_indices], dims = [1, 2])
    if rotation == 90
        magnitude_images = mapslices(x -> rotr90(x), magnitude_images, dims = [1, 2])
    elseif rotation == 180
        magnitude_images = mapslices(x -> rot180(x), magnitude_images, dims = [1, 2])
    else
        magnitude_images = mapslices(x -> rotl90(x), magnitude_images, dims = [1, 2])
    end
    magnitude_mosaic = mosaicview(magnitude_images, nrow = Int(floor(sqrt(num_slices))), npad = 5, rowmajor = true, fillvalue = 0)

    p1 = Plots.heatmap(magnitude_mosaic, show = false, title = "|Images|", color = :grays, aspect_ratio = 1, yflip=true, colorbar=:none, grid=:none, showaxis=:false)
    display(p1)

    phase_images = angle.(images[:, :, reordered_slice_indices, 1, 1])
    if rotation == 90
        phase_images = mapslices(x -> rotr90(x), phase_images, dims = [1, 2])
    elseif rotation == 180
        phase_images = mapslices(x -> rot180(x), phase_images, dims = [1, 2])
    else
        phase_images = mapslices(x -> rotl90(x), phase_images, dims = [1, 2])
    end
    phase_mosaic = mosaicview(phase_images, nrow = Int(floor(sqrt(num_slices))), npad = 5, rowmajor = true, fillvalue = 0)

    p2 = Plots.heatmap(phase_mosaic, show = false, title = "∠ Images", color = :plasma, aspect_ratio = 1, yflip=true, colorbar=:none, grid=:none, showaxis=:false)
    display(p2)

    b0_map_images = mapslices(x -> x, b0, dims = [1, 2])
    if rotation == 90
        b0_map_images = mapslices(x -> rotr90(x), b0, dims = [1, 2])
    elseif rotation == 180
        b0_map_images = mapslices(x -> rot180(x), b0, dims = [1, 2])
    else
        b0_map_images = mapslices(x -> rotl90(x), b0, dims = [1, 2])
    end
    b0_map_mosaic = mosaicview(b0_map_images[:, :, reordered_slice_indices], nrow = Int(floor(sqrt(num_slices))), npad = 5, rowmajor = true, fillvalue = 0)

    p3 = Plots.heatmap(b0_map_mosaic, show = false, title = "B₀ Map Images", color = :plasma, aspect_ratio = 1, yflip=true, colorbar=:right, grid=:none, showaxis=:false)
    display(p3)
end

"""
    check_profiles(raw_data::RawAcquisitionData)
Sanity check of RawAcqData object by ploting all profiles to confirm its consistency with ISMRMRD file

# Arguments
* `raw_data::RawAcquisitionData` - RawAcquisitionData object
* `num_profiles_display` - The number of profiles to be displayed
"""
function check_profiles(raw_data; num_profiles_display = 128)

    for l = 1:num_profiles
        p1 = plot(abs.(raw_data.profiles[l].data[:, 1]))
        p2 = plot(angle.(raw_data.profiles[l].data[:, 1]))
    end

end

"""
    plot_sense_maps(sensitivity, num_channels; slice_index = 1)
Plots coil sensitivity maps from the channels, for a given number of num_channels plots on a given slice index.

# Arguments
* `sensitivity` - sensitivity maps, a 4D array: [nX, nY, nZ, nCoil]
* `num_channels` - number of coils to be displayed.
* `slice_index` - The index of the slice to be displayed (if multislice)
"""
function plot_sense_maps(sensitivity, num_channels; slice_index = 1)
    if ndims(sensitivity) == 4
        num_slices = size(sensitivity, 3)
        num_channels_total = size(sensitivity, 4)
    else
        err_msg = @sprintf("sensitivity must be a 4D array with a size of [nX, nY, nZ, nCoil]. Current input has %d dimensions.", ndims(sensitivity))
        error(err_msg)
    end

    if num_channels > num_channels_total
        err_msg = @sprintf("The number of coils to be displayed is %d, but total available coil number is %d.", num_channels, num_channels_total)
        error(err_msg)
    end

    if slice_index > num_slices
        err_msg = @sprintf("The index of slice to be displayed is %d, but total slice number is %d.", slice_index, num_slices)
        error(err_msg)
    end

    magnitude_mosaic = mosaicview((abs.(sensitivity[:, :, slice_index, :])), nrow = Int(floor(sqrt(num_channels))), npad = 5, rowmajor = true, fillvalue = 0)
    p4 = Plots.heatmap(magnitude_mosaic, show = false, title = "|Sensitivity|", color = :gnuplot2, aspect_ratio = 1, yflip=true, colorbar=:none, grid=:none, showaxis=:false)

    phase_mosaic = mosaicview((angle.(sensitivity[:, :, slice_index, :])), nrow = Int(floor(sqrt(num_channels))), npad = 5, rowmajor = true, fillvalue = 0)
    p5 = Plots.heatmap(phase_mosaic, show = false, title = "∠ Sensitivity", color = :plasma, aspect_ratio = 1, yflip=true, colorbar=:none, grid=:none, showaxis=:false)

    display(p4)
    display(p5)

end

# "WIP: Plots trajectory and Data, doesn't work currently"
function plot_traj_and_data(acq)

    for l = 1:length(acq.traj)

        freq_encode[l, :] = acq.traj[l].nodes[1, :]
        phase_encode[l, :] = acq.traj[l].nodes[2, :]
        k_space_signal[l, :] = acq.kdata[l, :, 1]

    end

end

## PREPROCESSING

"""
    calculate_b0_maps(me_data, slices, echotime_1, echotime_2)

Calculate B₀ map from the two images with different echo times via their phase difference (phase of img_TE2.*conj(img_TE1))

# Arguments
* `me_data`                          - [nX nY nZ 2 num_coils] 5D image array, 4th dim is echo time
* `slices::NTuple{num_slices,Int}`   - slice index vector (tuple?) for which map is computed
* `echotime_1::AbstractFloat`        - TE1 [ms]
* `echotime_2::AbstractFloat`        - TE2 [ms]
"""
function calculate_b0_maps(me_data, slices, echotime_1, echotime_2)

    # b0_maps = mapslices(x -> rotl90(x),ROMEO.unwrap(angle.(me_data[:,:,slices,2,1].*conj(me_data[:,:,slices,1,1]))),dims=(1,2))./((7.38-4.92)/1000)
    b0_maps =
        mapslices(x -> x, ROMEO.unwrap(angle.(me_data[:, :, slices, 2, 1] .* conj(me_data[:, :, slices, 1, 1]))), dims = (1, 2)) ./
        ((echotime_2 - echotime_1) / 1000)

end

"""
    get_slice_order(r::RawAcquisitionData, sliceNum::Int, startProfile::Int, incProfile:Int)
Return a array of slice order index with ascending order of Z position.

e.g. For an interleaved pattern of slice position in RawAcquisitionData given below (in unit of mm):
[-7, -3, 1, 5, 9, -9, -5, -1, 3, 7], the output will be [6, 1, 7, 2, 8, 3, 9, 4, 10, 5]

# Arguments
* `r::RawAcquisitionData`       - A RawAcquisitionData that directly reads from original MRD file
* `sliceNum::Int`               - Total slice number that included in the RawAcquisitionData
* `startProfile::Int`           - Starting index of the profile in the RawAcqData for the first valid slice to be processed
* `incProfile::Int`             - Increment of profile index for the next valid slices

# Output
* `orderedIndex`                    - array with slice index of RawAcquisitionData with ascending order of position in Z.
"""
function get_slice_order(r::RawAcquisitionData, sliceNum::Int, startProfile::Int, incProfile::Int)
    origZPos = zeros(Float32, sliceNum)
    for m = 1 : sliceNum
        profileIndex = startProfile + (m - 1) * incProfile
        origZPos[m] = r.profiles[profileIndex].head.position[3]
    end
    return sortperm(origZPos)
end

"""
    sync_traj_and_data!(a::AcquisitionData)
Synchronizes k-space trajectory and sampled data as they do not usually have a common sampling rate

# Arguments
* `raw_data::RawAcquisitionData` - RawAcquisitionData object
* `traj::Trajectory` - Trajectory object to be synchronized with data contained in raw_data
* `idx_crop::Int` - Trajectory and Data may contain samples we don't want in the recon, usually at the end of acquisition. Ignore samples after idx_crop
* `interleave::Int` - index of interleave
"""
function sync_traj_and_data!(raw_data, traj, idx_crop, interleave)

    # get number of gradient samples
    num_gradient_samples = traj.numSamplingPerProfile # cannot avoid camelcase as it is in MRIBase

    # get vector of gradient samples which pertain to one interleave of the trajectory
    interleave_extraction_vector = num_gradient_samples * (interleave - 1) .+ (1:num_gradient_samples)

    # Read the trajectory nodes into the rawAcquisitionData type field (.traj)
    for l = 1:length(raw_data.profiles)
        raw_data.profiles[l].traj = traj.nodes[:, interleave_extraction_vector]
    end

    # define dwell times for trajectory (dt_k) and signal sampling (dt_s)
    dt_s = 2 * 10^(-6) # [s]
    dt_k = 10 * 10^(-6) # [s]

    # Go through every profile (this means every slice in the MRIReco.jl convention for multislice, multiTE, and diffusion scans)
    for l = 1:length(raw_data.profiles)

        # Get size of trajectory and signal vectors
        num_data_samples = size(raw_data.profiles[l].data, 1)
        num_kspace_samples = size(raw_data.profiles[l].traj, 2)

        # Define time vectors for signal and trajectory
        t_s = (0:num_data_samples-1) * dt_s
        t_k = (0:num_kspace_samples-1) * dt_k

        # Interpolate trajectory onto the same sample times as the sampled signal
        trajectory_interpolated_x = Spline1D(t_k, raw_data.profiles[l].traj[1, :], w = ones(length(raw_data.profiles[l].traj[1, :])), k = 3, bc = "zero")
        trajectory_interpolated_y = Spline1D(t_k, raw_data.profiles[l].traj[2, :], w = ones(length(raw_data.profiles[l].traj[2, :])), k = 3, bc = "zero")

        # Concatenate the trajectory node kx and ky positions
        adjusted_trajectory = vcat(trajectory_interpolated_x(t_s)', trajectory_interpolated_y(t_s)')

        # Crop the data and trajectory to avoid return-to-center of traj, and also set trajectory as upsampled trajectory adjusted_trajectory
        raw_data.profiles[l].traj = adjusted_trajectory[:, 1:idx_crop]
        raw_data.profiles[l].data = raw_data.profiles[l].data[1:idx_crop, :]

    end

    # Return the vector of sampling times
    return dt_s * (0:idx_crop-1)

end


"""
    do_k0_correction!(raw_data, k0_phase_modulation, interleave)
Applies phase modulation due to 0th-order field fluctuations during the acquisition

# Arguments
* `raw_data::RawAcquisitionData` - RawAcquisitionData object
* `k0_phase_modulation::Matrix{Complex{T}}` - Vector containing phase modulation measurements
* `interleave::Int` - index of interleave
"""
function do_k0_correction!(raw_data, k0_phase_modulation, interleave)

    # Get number of samples for the k0 phase modulation (should be same size as the trajectory BEFORE resampling)
    num_k0_samples = size(k0_phase_modulation, 1)

    # define dwell times for phase modulation (dt_k) and signal sampling (dt_s)
    dt_s = 2 * 10^(-6) # [s]
    dt_k = 10 * 10^(-6) # [s]

    # Go through every profile (this means every slice in the MRIReco.jl convention for multislice, multiTE, and diffusion scans)
    for l = 1:length(raw_data.profiles)

        # Get size of data (signal samples) and size of k0 modulation
        num_data_samples = size(raw_data.profiles[l].data, 1)
        num_kspace_samples = num_k0_samples

        # Define time vectors for k0 and signal sampling times
        t_s = (0:num_data_samples-1) * dt_s
        t_k = (0:num_kspace_samples-1) * dt_k

        # interpolate k0 to the time basis of the signal
        k0_interpolant = Spline1D(t_k, k0_phase_modulation[:, interleave], w = ones(num_k0_samples), k = 3, bc = "zero")
        k0_interpolated = k0_interpolant(t_s)

        # modulate the data by the k0 modulation by multiplying with e^(i*k0) where k0 is in radians
        raw_data.profiles[l].data = raw_data.profiles[l].data .* exp.(1im .* k0_interpolated)

        # # Visualization of Phase Modulation
        # figure("Phase Modulation")
        # plot(t_s, angle.(exp.(1im .* k0_interpolated)))
        # xlabel("Time [s]")
        # ylabel("k₀ [rad]")
        # title("B₀ Eddy Current Fluctuation During Readout ")

        # plot(t_s, angle.(exp.(1im .* k0_interpolated)), show = true, title = "B₀ Eddy Current Fluctuation During Readout ")

    end

end


"""
    adjust_header!(raw::RawAcquisitionData, recon_size, num_samples, interleave_number, single_slice)
Adjusts the header data for each interleave and slice of spiral diffusion RawAcquisitionData

# Arguments
* `raw::RawAcquisitionData` - RawAcquisitionData object
* `recon_size::Vector` - Reconstruction matrix size
* `num_samples::Int` - Number of samples per interleave
* `interleave_number::Int` - Index of interleave for multi-shot acquisitionNumbers
* `single_slice::Bool` - flag for single-slice reconstruction/acquisition
"""
function adjust_header!(raw, recon_size, num_samples, interleave_number, single_slice)

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
        raw.profiles[l].head.number_of_samples = num_samples

        # Set the non-standard encode step (interleave dimension) into the encode step 1 field
        # raw.profiles[l].head.idx.kspace_encode_step_1 = 0
        raw.profiles[l].head.idx.kspace_encode_step_1 = interleave_number - 1 # IF MULTI-INTERLEAVE

        # Set the slice index header for the first slice as 0, otherwise the trajectory will be all zeros after converting from RawData to AcqData
        if l == 1
            raw.profiles[l].head.idx.slice = 0
        end

        # Set center sample to 0 (only for spiral scans)
        raw.profiles[l].head.center_sample = 0

    end

    # Set encoding size to the recon_size
    raw.params["encodedSize"] = [recon_size[1], recon_size[2], 1]

end

"""
    check_acquisition_nodes!(a::AcquisitionData)
Validates processed AcquisitionData object to make sure that |kᵢ| < 0.5 ∀ i ∈ [1, Nₛ]

# Arguments
* `a::AcquisitionData` - AcquisitionData object
"""
function check_acquisition_nodes!(a::AcquisitionData)

    a.traj[1].nodes[abs.(a.traj[1].nodes[:]).>0.5] .= 0.5

end


"""
    validate_siemens_mrd!(r::RawAcquisitionData)
Validates RawAcquisitionData object created from ISMRMRD format object

# Arguments
* `r::RawAcquisitionData` - RawAcquisitionData object
"""
function validate_siemens_mrd!(r::RawAcquisitionData)

    @info "Validating Siemens converted data"

    ## FOV CHECK:

    if maximum(r.params["encodedFOV"]) > 0.8 # FOV should never be greater than the bore size in [m]

        @info "FOV was recorded in [mm]! Changing to [m]!"
        r.params["encodedFOV"] = r.params["encodedFOV"] ./ 1000

    end

end

"""
    validate_acq_data!(a::AcquisitionData)
Validates processed AcquisitionData object after manipulation, etc...

# Arguments
* `a::AcquisitionData` - AcquisitionData object
"""
function validate_acq_data!(a::AcquisitionData)

    ## Dimensions CHECK:

    # TODO add dimension check that the k-space encoding counters are set properly:
    # kdata dimensions: dim1:=contrast/echo | dim2:=slices | dim3:=repetitions 
    # kdata element dimensions: dim1:=kspace nodes | dim2:=channels/coils

    # permutedims(a.kdata, [3, 2, 1])
    check_acquisition_nodes!(a)

end

"""
    preprocess_cartesian_data(r::RawAcquisitionData, do_save; filename = "data/processed_cartesian_file.h5")
Prepares Cartesian for reconstruction

# Arguments
* `r::RawAcquisitionData{T}`   - RawAcquisitionData object
* `do_save::Boolean`           - Save the processed Cartesian data as a HDF5 file
* `filename`                   - filename to save the preprocessed data
"""
function preprocess_cartesian_data(r::RawAcquisitionData, do_save; filename = "data/processed_cartesian_file.h5")

    remove_oversampling!(r)

    # Convert rawAcquisitionData object to an AcquisitionData object (these can be reconstructed)
    cartesian_acq_data = AcquisitionData(r, estimateProfileCenter = true) # cannot avoid camel case as defined by MRIBase

    ## Properly arrange data from the converted siemens file
    validate_acq_data!(cartesian_acq_data)

    if do_save

        raw = RawAcquisitionData(cartesian_acq_data)

        # Since the data should generally have 3D information when saved, we make sure 2D data is appropriately stored as 3D data with a singleton dimension
        if length(raw.params["encodedSize"]) == 2
            e_sz = raw.params["encodedSize"]
            raw.params["encodedSize"] = [e_sz[1], e_sz[2], 1]
        end

        raw.params["TE"] = r.params["TE"]

        # raw.params = headerCopy
        fout = ISMRMRDFile(filename)
        save(fout, raw)

    end

    return cartesian_acq_data

end

"""
    remove_oversampling!(raw::RawAcquisitionData)
Removes 2x readout oversampling in raw data along read-out dimension.

# Arguments
* `raw::RawAcquisitionData{T}`          - RawAcquisitionData object
"""
function remove_oversampling!(raw::RawAcquisitionData)

    dimension_index = 1
    num_data_samples = raw.params["encodedSize"][dimension_index]
    index_crop_fov = convert(Vector{Int32}, [1:floor(num_data_samples / 4); ceil(3 / 4 * num_data_samples + 1):num_data_samples])

    # For every profile in the acquisition
    for profile_index = 1:length(raw.profiles)

        # IFFT to image space, crop, FFT back to k-space
        ifft!(raw.profiles[profile_index].data, dimension_index)
        raw.profiles[profile_index].data = fft!(raw.profiles[profile_index].data[index_crop_fov, :], dimension_index)

    end

    # halve encoding size of first dimension
    raw.params["encodedSize"][dimension_index] /= 2
    raw.params["encodedFOV"][dimension_index] /= 2

end


"""
    merge_raw_interleaves(params, output_raw)
Merges multiple interleave data together from individually acquired interleave scans

# Arguments
* `params`          - Dictionary
* `output_raw`      - Bool
"""
function merge_raw_interleaves(params, output_raw)

    # Get the other interleave indexes other than the one asked for
    interleave_complement = [x for x ∈ 1:params[:num_interleaves] if x ∉ params[:interleave]]

    # @info "indices = $interleave_complement" #DEBUG

    # read in the data file from the ISMRMRD format
    data_file = ISMRMRDFile(params[:interleave_data_filenames][params[:interleave]])

    # Read in the gradient file
    input_trajectory = read_gradient_text_file(params[:traj_filename], params[:recon_size], params[:delay])

    # Read in raw data from the data_file
    raw_data = RawAcquisitionData(data_file)

    # delete everything that is not a chosen excitation (for efficiency)
    indices = 1:length(raw_data.profiles)
    ic = [x for x ∈ indices if x ∉ params[:excitations]]
    deleteat!(raw_data.profiles, ic)

    # @info "indices = $ic" #DEBUG

    # set up time vector for tracking all of the interleaves
    time_track_vector = []

    # synchronize trajectory data and the kspace data
    times = sync_traj_and_data!(raw_data, input_trajectory, params[:num_samples], params[:interleave])

    # adjust the header so that each diffusion direction is considered as a contrast instead of a repetition
    # adjust_header!(raw_data, params[:recon_size], params[:num_samples], params[:interleave],params[:single_slice])
    adjust_header!(raw_data, params[:recon_size], params[:num_samples], 1, params[:single_slice])

    # add the times to the time tracking vector
    append!(time_track_vector, times)

    # Repeat the above steps for each interleave, adjusting the times and headers appropriately
    if params[:do_multi_interleave]

        for l in interleave_complement

            # read in separate interleave data file
            data_file_temp = ISMRMRDFile(params[:interleave_data_filenames][l])
            raw_data_temp = RawAcquisitionData(data_file_temp)
            deleteat!(raw_data_temp.profiles, ic) # delete profiles which aren't needed

            # synchronize the trajectory from the gradient file and the data from the raw data file for the interleave
            times_temp = sync_traj_and_data!(raw_data_temp, input_trajectory, params[:num_samples], l)

            # adjust the header to reflect the arrangement of data expected by MRIReco.jl's reconstruction function
            adjust_header!(raw_data_temp, params[:recon_size], params[:num_samples], l, params[:single_slice])

            # append the important data (the profile and the sampling times) to the raw Data file created out of this look
            append!(raw_data.profiles, deepcopy(raw_data_temp.profiles))
            append!(time_track_vector, deepcopy(times_temp))

        end

        # if there is the choice to do odd or opposing interleaves, add the 2nd interleave
    elseif params[:do_odd_interleave]

        data_file_temp = ISMRMRDFile(params[:interleave_data_filenames][3])
        raw_data_temp = RawAcquisitionData(data_file_temp)
        deleteat!(raw_data_temp.profiles, ic)

        times_temp = sync_traj_and_data!(raw_data_temp, input_trajectory, params[:num_samples], 3)

        adjust_header!(raw_data_temp, params[:recon_size], params[:num_samples], 2, params[:single_slice])

        append!(raw_data.profiles, copy(raw_data_temp.profiles))
        append!(time_track_vector, times_temp)

    end

    # Correct raw data reading if necessary
    for ll = 1:length(raw_data.profiles)

        raw_data.profiles[ll].head.trajectory_dimensions = 2

    end

    if output_raw

        # return the RawAcquisition data object (missing some corrections but maybe better for some users)
        return raw_data

    else

        # converting raw_data to AcquisitionData
        @info "Converting RawAcquisitionData to AcquisitionData"
        acq_data = AcquisitionData(raw_data, estimateProfileCenter = true)

        ## Assume all of the slices share a trajectory
        for l = 1:length(acq_data.traj)

            acq_data.traj[l].times = time_track_vector # set times to the total time vector
            acq_data.traj[l].TE = 0.00 # set the TE to 0
            acq_data.traj[l].AQ = times[end] # set the acquisition time to the last element of the time vector (should be the latest time)
            acq_data.traj[l].circular = true # set whether to use a circular filter on the kspace data

        end

        for l = 1:length(acq_data.subsampleIndices) # Cannot avoid camelcase

            acq_data.subsampleIndices[l] = 1:size(acq_data.traj[l].nodes, 2) 

        end

        # return the acquisition data object with everything corrected
        return acq_data

    end

end

"""
    apply_girf!(a::AcquisitionData{T}, g::GirfApplier)
Applies the GIRF to the trajectories inside of a::AcquisitionData

# Arguments
* `a::AcquisitionData{T}`          - AcquisitionData object
* `g::GirfApplier`                 - GirfApplier object containing GIRF definition
"""
function apply_girf!(a::AcquisitionData{T}, g::GirfApplier) where {T}

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

        num_profiles = a.traj[l].numProfiles
        num_samples = a.traj[l].numSamplingPerProfile
        nodes = a.traj[l].nodes
        times = a.traj[l].times
        old_nodes = a.traj[l].nodes

        # loop over all profiles in a trajectory
        for profile = 1:num_profiles

            interleave_extractor = num_samples * (profile - 1) .+ (1:num_samples)
            interleave_nodes = nodes[:, interleave_extractor]
            interleave_times = times[interleave_extractor]

            dt = interleave_times[1] - interleave_times[2]

            interleave_gradients = nodes_to_gradients(interleave_nodes; dwellTime = dt, reconSize = S, FOV = F)

            # loop over trajectory dimensions
            for dim = 1:size(interleave_gradients, 1)

                corrected_gradients = apply_girf(g, interleave_gradients[dim, :], interleave_times, interleave_times, dim) # THESE ARE ALL VECTORS SO INPUT orientation (column/row major ordering) doesn't matter
                interleave_gradients[dim, :] = corrected_gradients'

            end

            interleave_nodes = gradients_to_nodes(interleave_gradients; dwellTime = dt, reconSize = S, FOV = F)
            nodes[:, interleave_extractor] = interleave_nodes

        end

        a.traj[l].nodes = nodes

    end

end

"""
    apply_k0!(a::AcquisitionData{T}, g::GirfApplier)
Applies the K0 modulation due to imaging gradients to the data inside of a::AcquisitionData

# Arguments
* `a::AcquisitionData{T}`          - AcquisitionData object
* `g::GirfApplier`                 - GirfApplier containing GIRF definition
"""
function apply_k0!(a::AcquisitionData{T}, g::GirfApplier) where {T}

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

        num_profiles = a.traj[l].numProfiles
        num_samples = a.traj[l].numSamplingPerProfile
        nodes = a.traj[l].nodes
        times = a.traj[l].times
        old_nodes = a.traj[l].nodes

        # loop over all profiles in a trajectory
        for profile = 1:num_profiles

            interleave_extractor = num_samples * (profile - 1) .+ (1:num_samples)
            interleave_nodes = nodes[:, interleave_extractor]
            interleave_times = times[interleave_extractor]

            dt = interleave_times[1] - interleave_times[2]

            interleave_gradients = nodes_to_gradients(interleave_nodes; dwellTime = dt, reconSize = S, FOV = F)

            k0_correction = ones(size(interleave_gradients))

            # loop over all trajectory dims
            for dim = 1:size(interleave_gradients, 1)

                k0_correction[dim, :] = apply_girf(g, interleave_gradients[dim, :], interleave_times, interleave_times, dim) # THESE ARE ALL VECTORS SO INPUT orientation (column/row major ordering) doesn't matter

            end

            final_correction = sum(k0_correction, dims = 1) #back to radians!

            a.kdata[l][interleave_extractor, :] = a.kdata[l][interleave_extractor, :] .* exp.(-1im .* final_correction')

            # # Visualization of Phase Modulation
            # figure("Phase Modulation 2")
            # plot(vec(interleave_times), vec(angle.(exp.(1im .* final_correction))))
            # xlabel("Time [s]")
            # ylabel("k₀ [rad]")
            # title("B₀ Eddy Current Fluctuation During Readout ")

            # plot(interleave_times, angle.(exp.(1im .* final_correction')),show=true,title="B₀ Eddy Current Fluctuation During Readout ") #DEBUG

        end

    end

end

# ## Calibrate the phase from individual interleaves
# function calibrateAcquisitionPhase!(a::AcquisitionData)

#     for l = 1:length(a.traj)

#         num_profiles = a.traj[l].numProfiles
#         num_samples = a.traj[l].numSamplingPerProfile
#         nodes = a.traj[l].nodes
#         times = a.traj[l].times
#         old_nodes = a.traj[l].nodes

#         initialInterleavePhase = angle.(a.kdata[l][1,:])'

#         for profile = 2:num_profiles

#             interleave_extractor = num_samples*(profile-1) .+ (1:num_samples)

#             initialProfilePhase = angle.(a.kdata[l][interleave_extractor[1],:])'

#             @info size(initialInterleavePhase)
#             @info size(a.kdata[l][interleave_extractor,:])

#             a.kdata[l][interleave_extractor,:] .*= exp.(-1im * (initialInterleavePhase - initialProfilePhase))

#         end    
#     end

# end

## Input/Output, File handling

"""
    save_map(filename, calib_map, resolution_mm; offset_mm = [0.0, 0.0, 0.0], do_split_phase::Bool = false, do_normalize::Bool = true)
Saves calibration maps (sensitivity or B₀) as 4D NIfTI file(s)

For complex-valued data, magnitude and phase can be split into separate files
# Arguments
* `filename::String`              - string filename with extension .nii, example "sensemap.nii"
* `calib_map`                     - [nX nY nZ {nChannels}] 4-D sensitivity or 3D B₀ map array 
* `resolution_mm`                 - resolution in mm, 3 element vector, e.g., [1.0, 1.0, 2.0]
* `offset_mm`                     - isocenter offset in mm, default: [0.0, 0.0, 0.0]
* `do_split_phase::Bool = false`  - if true, data is saved in two nifti files with suffix "_magn" and "_phase", respectively
                                  to enable display in typical NIfTI viewers
* `do_normalize::Bool = true`     - if true, normalize the image by its magnitude maxima
"""
function save_map(filename, calib_map, resolution_mm; offset_mm = [0.0, 0.0, 0.0], do_split_phase::Bool = false, do_normalize::Bool = true)

    # multiplication with 1000 should no longer be necessary after MRIReco 0.7.1
    spacing = 1000.0 .* resolution_mm .* Unitful.mm
    offset = 1000.0 .* offset_mm .* Unitful.mm

    if ndims(calib_map) >= 4 # multi-coil calib_map, e.g., sensitivity, or recon, but we can only store the first 4 dims in a Nifti
        I = reshape(calib_map, size(calib_map, 1), size(calib_map, 2), size(calib_map, 3), size(calib_map, 4), 1, 1)
    else
        I = reshape(calib_map, size(calib_map, 1), size(calib_map, 2), size(calib_map, 3), 1, 1, 1)
    end

    # scale to max 1
    if do_normalize
        I /= maximum(abs.(I))
    end

    # AxisArray Constructor
    im = AxisArray(
        I,
        Axis{:x}(range(offset[1], step = spacing[1], length = size(I, 1))),
        Axis{:y}(range(offset[2], step = spacing[2], length = size(I, 2))),
        Axis{:z}(range(offset[3], step = spacing[3], length = size(I, 3))),
        Axis{:coils}(1:size(I, 4)),
        Axis{:echos}(1:size(I, 5)),
        Axis{:repetitions}(1:size(I, 6)),
    )

    # if separate mag and phase are desired, save them separately
    if do_split_phase

        magnitude_filename = splitext(filename)[1] * "_magn.nii"
        saveImage(magnitude_filename, map(abs, im)) # map is needed, because abs.(im) would convert AxisArray back into basic array

        phase_filename = splitext(filename)[1] * "_phase.nii"
        saveImage(phase_filename, map(angle, im))

    else

        saveImage(filename, im)

    end

end

"""
    load_map(filename; do_split_phase::Bool = false)
Load calibration maps (sensitivity or B₀) from 4D NIfTI file(s)

For complex-valued data, magnitude and phase parts are stored in two files with suffix "_magn" and "_phase".

# Arguments
* `filename::String`            - string filename with extension .nii, example "sensemap.nii"
* `do_split_phase::Bool=false`    - if true, data is saved in two nifti files with suffix "_magn" and "_phase", respectively
                                  to enable display in typical NIfTI viewers
# Output
* `calib_map`                    - [nX nY nZ {nChannels}] 4-D sensitivity or 3D B₀ map array 
"""
function load_map(filename; do_split_phase::Bool = false)

    # if separate mag and phase are saved, load and combine them
    if do_split_phase

        magnitude_filename = splitext(filename)[1] * "_magn.nii"
        magnitude_image = loadImage(magnitude_filename) # map is needed, because abs.(im) would convert AxisArray back into basic array

        phase_filename = splitext(filename)[1] * "_phase.nii"
        phase_image = loadImage(phase_filename)

        calib_map = (magnitude_image.data) .* exp.(1im .* (phase_image.data))

    else

        I = loadImage(filename)
        calib_map = I.data

    end

    # squeeze singleton dimensions of 6-dim array
    calib_map = dropdims(calib_map, dims = tuple(findall(size(calib_map) .== 1)...))

    return calib_map

end

"""
    shift_kspace!(acqdata, shift)
This function applys additional phase ramps to k-space data to achive a given shift of center of image FOV in X and Y directions.

Perhaps this should be called shift_fov; however, since this function is modifying kspace data, it is named shift_kspace for now.

# Arguments
* `acqdata::AcquisitionData{T}`     - AcquisitionData object to be modified
* `shift::AbstractVector`           - Vector containing shift with size [shift_X, shift_Y]
"""
function shift_kspace!(acqdata, shift)

    num_slices = numSlices(acqdata)
    num_repetitions, num_contrasts = numRepetitions(acqdata), numContrasts(acqdata)

    smat = prod(exp.(1im .* acqdata.traj[1].nodes[:, acqdata.subsampleIndices[1]] .* shift .* 2 .* pi), dims = 1)

    for slice = 1:num_slices
        for contr = 1:num_contrasts
            for rep = 1:num_repetitions
                acqdata.kdata[contr, slice, rep] = acqdata.kdata[contr, slice, rep] .* smat'
            end
        end
    end

end

function run_cartesian_recon(params_general)

    @info "Loading Data Files"

    b0_filename = params_general[:map_scan_fullpath];

    # filename for preprocessed data (remove oversampling, permute dimensions wrt MRIReco)
    processed_filename = params_general[:processed_map_scan_fullpath]

    # Set the data file name (Change this for your own system)
    cartesian_data_file = ISMRMRDFile(b0_filename)

    # read in the raw data from the ISMRMRD file into a RawAcquisitionData object
    raw = RawAcquisitionData(cartesian_data_file)

    if params_general[:do_process_map_scan] 
        # Preprocess Data and save!
        preprocess_cartesian_data(raw::RawAcquisitionData, true; filename = processed_filename)

    end

    # Load preprocessed data!
    processed_cartesian_data_file = ISMRMRDFile(processed_filename)

    raw_new = RawAcquisitionData(processed_cartesian_data_file)
    cartesian_acq_data = AcquisitionData(raw_new, estimateProfileCenter = true)

    # Define coils and slices
    num_coils = size(cartesian_acq_data.kdata[1], 2)
    num_slices = numSlices(cartesian_acq_data)

    # Get the order of slices from the RawAcqData header
    slice_idx_array_cartesian = get_slice_order(raw, num_slices, 1, 1)

    # Define TEs 
    # Echo times for field map raw data, in ms
    TE1 = raw_new.params["TE"][1]
    TE2 = raw_new.params["TE"][2]

    # shift FOV to middle :) 
    #TODO: in MRIReco v0.7, try: correctOffset(cartesian_acq_data, [0 -20 0])
    shift_kspace!(cartesian_acq_data, params_general[:fov_shift]) # amount of FOV shift; in unit of number of voxels in [x,y] direction

    #changeFOV!(cartesian_acq_data,[1.5,1.5])

    ## Don't have to recalculate sense maps for both scans but possibly it could make a
    #  difference in Diffusion scans

    @info "Calculating Sense Maps"

    # from reading docstring, min thresholding
    cartesian_sensitivity = espirit(cartesian_acq_data, (6, 6), 30, eigThresh_1 = 0.01, eigThresh_2 = 0.9)

    # from Lars (7T spirals)
    #cartesian_sensitivity = espirit(cartesian_acq_data,(6,6),30,eigThresh_1=0.02, eigThresh_2=0.98)
    # from Alexander (Phantom?)
    #cartesian_sensitivity = espirit(cartesian_acq_data,(4,4),12,eigThresh_1=0.01, eigThresh_2=0.98)

    # normalize for consistency with saving/loading and better ranges of reconstruction values
    cartesian_sensitivity /= maximum(abs.(cartesian_sensitivity))

    res_x = fieldOfView(cartesian_acq_data)[1] ./ size(cartesian_sensitivity)[1]
    res_y = fieldOfView(cartesian_acq_data)[2] ./ size(cartesian_sensitivity)[2]
    res_z = fieldOfView(cartesian_acq_data)[3] .* (1 + params_general[:slice_distance_factor_percent] ./ 100.0) # for 2D only, since FOV[3] is slice thickness then, but gap has to be observed
    resolution_mm = (res_x, res_y, res_z)

    # save SENSE maps
    if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
        # TODO: use correct slice order everywhere, e.g., when saving/loading maps for spiral recon
        save_map(params_general[:sensitivity_save_fullpath], cartesian_sensitivity[:, :, slice_idx_array_cartesian, :], resolution_mm; do_split_phase = true)
    end

    ## Parameter dictionary definition for reconstruction

    @info "Setting Parameters"
    params_cartesian = Dict{Symbol,Any}() # instantiate dictionary
    params_cartesian[:reco] = "multiCoil" # choose multicoil reconstruction
    params_cartesian[:reconSize] = (cartesian_acq_data.encodingSize[1], cartesian_acq_data.encodingSize[2]) # set recon size to be the same as encoded size
    params_cartesian[:regularization] = "L2" # choose regularization for the recon algorithm
    params_cartesian[:λ] = 1.e-2 # recon parameter (there may be more of these, need to dig into code to identify them for solvers other than cgnr)
    params_cartesian[:iterations] = params_general[:num_recon_iterations] # number of CG iterations
    params_cartesian[:solver] = "cgnr" # inverse problem solver method
    params_cartesian[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false) # turn on store solutions if you want to see the reconstruction convergence (uses more memory)
    params_cartesian[:senseMaps] = ComplexF32.(cartesian_sensitivity) # set sensitivity map array


    ## Call the reconstruction function

    @info "Performing Reconstruction"
    @time cartesian_reco = reconstruction(cartesian_acq_data, params_cartesian)

    # save Map recon (multi-echo etc.)
    if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
        save_map(params_general[:map_save_fullpath], cartesian_reco.data[:, :, slice_idx_array_cartesian, :, :, :], resolution_mm; do_split_phase = true)
    end

    ## Calculate B0 maps from the acquired images (if two TEs)
    @info "Calculating B0 Maps"
    slices = 1:length(slice_idx_array_cartesian)
    b0_maps = zeros(200, 200, num_slices)
    b0_method = "2D_2008" # Can be "Simple","2D_2008" or "3D_2020" (How do we incorporate this into the recon demo?)

    if b0_method == "2D_2008"

        b0_maps = estimate_b0_maps(cartesian_reco.data, slices, TE1, TE2, true; β = params_general[:b0_map_beta], reltol = 1e-4)

    elseif b0_method == "3D_2020"

        b0_maps, times, out = b0map(
            reshape(cartesian_reco.data[:, :, :, :, 1, 1], (200, 200, num_slices,1, 2)),
            (TE1 / 1000.0, TE2 / 1000.0); 
            smap = ones(ComplexF32,size(b0_maps,1),size(b0_maps,2),size(b0_maps,3),1),
            track =true,
            chat = true,
            chat_iter = 2
        )

        b0_maps = 2 * pi * b0_maps

    end

    # save B0 map
    if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
        save_map(params_general[:b0_map_save_fullpath], b0_maps[:, :, slice_idx_array_cartesian], resolution_mm; do_normalize = false) # no normalization, we want absolute values for offres maps
    end


    if params_general[:do_plot_recon]
        @info "Plotting Cartesian Results (Sensitivity Maps and B0 Maps)"
        # plot_sense_maps(sensitivity)
        plot_reconstruction(cartesian_reco[:, :, slice_idx_array_cartesian, 1], 1:size(cartesian_reco, 3), b0_maps[:, :, slice_idx_array_cartesian], is_slice_interleaved = false, rotation = 180)
    end

    # cleanup unused file
    if !params_general[:do_save_processed_map_scan] && params_general[:do_process_map_scan]
        rm(processed_filename)
    end

    @info "Successfully Completed cartesian_reconstruction"
    
end