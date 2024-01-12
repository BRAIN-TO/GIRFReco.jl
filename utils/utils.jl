export plot_reconstruction,
    plot_sense_maps,
    calculate_b0_maps,
    get_slice_order,
    sync_traj_and_data!,
    do_k0_correction!,
    adjust_header!,
    check_acquisition_nodes!,
    validate_siemens_mrd!,
    validate_acq_data!,
    preprocess_cartesian_data,
    remove_oversampling!,
    merge_raw_interleaves,
    apply_girf!,
    apply_k0!,
    save_map,
    load_map,
    shift_kspace!

## Choose plotting backend to be PlotlyJS!
# plotlyjs()

## General Plotting function for the reconstruction

# "Mosaic-plots reconstruction for selected slices and corresponding B0 map"
"""
    plot_reconstruction(images, slices_index, b0; fig_handles = [], is_slice_interleaved = false)
Plots the magnitude and phase of the reconstructed images for a given slice or slices, along with a B₀ map if applicable

# Arguments
* `images` - Complex-valued images reconstructed using MRIReco.jl
* `slices_index::Vector{Int}` - slices to plot
* `b0` - off-resonance map to plot along with images
* `fig_handles` - String vectors in size of [3,1] for titles of three figures (Magnitude & Phase of reconstructed images, and B0 maps)
* `is_slice_interleaved::Bool` - for 2D scanning, indicate this value as `true` to make sure the slice order on the displayed results is correct
* `rotation::Int` - Counterclock-wise rotation angle for each slice, should be a value from 0, 90, 180, 270 degrees
"""
function plot_reconstruction(images, slices_index, b0; fig_handles = [], is_slice_interleaved = false, rotation = 0)
    # plot()
    ## If we need to re-order all slices
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

    # Plot magnitude images (normalize)
    # if length(fig_handles) < 1
    #     figure("Magnitude Images")
    # else
    #     figure(fig_handles[1])
    # end

    # clf()

    magnitude_images = mapslices(x -> abs.(x) ./ maximum(abs.(x)), images[:, :, reordered_slice_indices], dims = [1, 2])
    if rotation == 90
        magnitude_images = mapslices(x -> rotr90(x), magnitude_images, dims = [1, 2])
    elseif rotation == 180
        magnitude_images = mapslices(x -> rot180(x), magnitude_images, dims = [1, 2])
    else
        magnitude_images = mapslices(x -> rotl90(x), magnitude_images, dims = [1, 2])
    end
    magnitude_mosaic = mosaicview(magnitude_images, nrow = Int(floor(sqrt(num_slices))), npad = 5, rowmajor = true, fillvalue = 0)

    heatmap(magnitude_mosaic, show = true, plot_title = "|Images|", plot_titlevspan = 0.1, color = :grays, aspectratio = :equal)
    # display(plot!())

    phase_images = angle.(images[:, :, reordered_slice_indices, 1, 1])
    if rotation == 90
        phase_images = mapslices(x -> rotr90(x), phase_images, dims = [1, 2])
    elseif rotation == 180
        phase_images = mapslices(x -> rot180(x), phase_images, dims = [1, 2])
    else
        phase_images = mapslices(x -> rotl90(x), phase_images, dims = [1, 2])
    end
    phase_mosaic = mosaicview(phase_images, nrow = Int(floor(sqrt(num_slices))), npad = 5, rowmajor = true, fillvalue = 0)

    heatmap(phase_mosaic, show = true, plot_title = "∠ Images", plot_titlevspan = 0.1, color = :plasma, aspectratio = :equal)
    #colorbar()
    # display(plot!())

    #gcf().suptitle("∠Images")

    # # Plot B0 maps
    # if length(fig_handles) < 3
    #     figure("B₀ Map Images")
    # else
    #     figure(fig_handles[3])
    # end

    # clf()
    # plot()
    b0_map_images = mapslices(x -> x, b0, dims = [1, 2])
    if rotation == 90
        b0_map_images = mapslices(x -> rotr90(x), b0, dims = [1, 2])
    elseif rotation == 180
        b0_map_images = mapslices(x -> rot180(x), b0, dims = [1, 2])
    else
        b0_map_images = mapslices(x -> rotl90(x), b0, dims = [1, 2])
    end
    b0_map_mosaic = mosaicview(b0_map_images[:, :, reordered_slice_indices], nrow = Int(floor(sqrt(num_slices))), npad = 5, rowmajor = true, fillvalue = 0)

    heatmap(b0_map_mosaic, show = true, plot_title = "B₀ Map Images", plot_titlevspan = 0.1, color = :plasma)
    # colorbar()
    # display(plot!())

    # gcf().suptitle("B₀ Maps [rad/s]")

    1

end

"Function plots all profiles in the acquisition to check consistency with ISMRMRD file"
function check_profiles(raw_data)

    num_profiles = 128 # Set to the number of profiles that you would like to see

    for l = 1:num_profiles
        p1 = plot(abs.(raw_data.profiles[l].data[:, 1]))
        p2 = plot(angle.(raw_data.profiles[l].data[:, 1]))
    end

end

"""
    plot_sense_maps!(sensitivity, num_channels)
Plots coil sensitivity maps from the channels, for a total of num_channels plots

# Arguments
* `sensitivity` - sensitivity maps
* `num_channels::Int` - number of coils (usually last dimension of sensitivity)
* `slice_index::Int` - The index of the slice to be displayed (if multislice)
"""
function plot_sense_maps(sensitivity, num_channels; slice_index = 1)
    num_slices = size(sensitivity, 3)
    if slice_index > num_slices
        err_msg = @sprintf("The index of slice to be displayed is %d, but total slice number is %d.", slice_index, num_slices)
        error(err_msg)
    end

    # # Magnitude maps
    # figure(@sprintf("Sensitivity Map Magnitude of Slice %d / %d", slice_index, num_slices)); clf(); for ch in 1:num_channels; subplot(8,4,ch); imshow((abs.(sensitivity[:,:,slice_index,ch])), cmap = "gray"); end;
    # subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
    # gcf()

    magnitude_mosaic = mosaicview((abs.(sensitivity[:, :, slice_index, :])), nrow = Int(floor(sqrt(num_channels))), npad = 5, rowmajor = true, fillvalue = 0)
    heatmap(magnitude_mosaic, show = true, plot_title = "|Sensitivity|", plot_titlevspan = 0.1, color = :gnuplot2)

    # # Phase maps
    # figure(@sprintf("Sensitivity Map Phase of Slice %d / %d", slice_index, num_slices)); clf(); for ch in 1:num_channels; subplot(8,4,ch); imshow(angle.(sensitivity[:,:,slice_index,ch]), cmap = "gray"); end;
    # subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
    # gcf()

    phase_mosaic = mosaicview((angle.(sensitivity[:, :, slice_index, :])), nrow = Int(floor(sqrt(num_channels))), npad = 5, rowmajor = true, fillvalue = 0)
    heatmap(phase_mosaic, show = true, plot_title = "∠ Sensitivity", plot_titlevspan = 0.1, color = :plasma)

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
    calculate_b0_maps(me_data,slices,echotime_1,echotime_2)

Calculate  B0 map from the two images with different echo times via their phase difference (obtained from imTE2.*conj(imTE1))
TODO have the b0 map calculation be capable of handling variable echo times
TODO2: Do we need this basic B0 map calculation or is it superseded by estimate_b0_maps?

# Arguments
* `me_data`                          - [nX nY nZ 2 num_coils] 5D image array, 4th dim echo time
* `slices::NTuple{num_slices,Int}`     - slice index vector (tuple?) for which map is computed
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
    get_slice_order(num_slices, is_slice_interleaved)

Returns array mapping from acquisition number to slice number (geometric position) (index_array[slice = 1:9] = [acquisitionNumbers])
TODO: Add ascending/descending options

# Arguments
* `num_slices::Int`                    - number of slices in total acquired stack (FOV)
* `is_slice_interleaved::Bool=true`   - if true, interleaved slice order is created, otherwise ascending slice order is returned
"""
function get_slice_order(num_slices; is_slice_interleaved::Bool = true)

    slice_index_array = 1:num_slices
    reordered_slice_indices = zeros(Int16, size(slice_index_array))
    if is_slice_interleaved && num_slices > 1
        reordered_slice_indices[1:2:end] = slice_index_array[1:Int(ceil(num_slices / 2))]
        reordered_slice_indices[2:2:end] = slice_index_array[Int(ceil(num_slices / 2) + 1):end]
    else
        reordered_slice_indices = slice_index_array
    end

    return reordered_slice_indices

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

        plot(t_s, angle.(exp.(1im .* k0_interpolated)), show = true, title = "B₀ Eddy Current Fluctuation During Readout ")

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

        # Set slice to 0 for singleslice, if it is not 0 then there will be an error
        if single_slice
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

    permutedims(a.kdata, [3, 2, 1])
    check_acquisition_nodes!(a)

end

"""
    preprocess_cartesian_data!(raw::RawAcquisitionData; dims = 1)
Prepares Cartesian for reconstruction

# Arguments
* `r::RawAcquisitionData{T}`          - RawAcquisitionData object
* `filename`                             - filename to save the preprocessed data to
"""
function preprocess_cartesian_data(r::RawAcquisitionData, do_save; filename = "data/testFile.h5")

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
    remove_oversampling!(raw::RawAcquisitionData; dims = 1)
Removes 2x readout oversampling in specified raw data dimensions by iFFT, cropping FOV and FFT

# Arguments
* `raw::RawAcquisitionData{T}`          - RawAcquisitionData object
* `dims`                                - dimension alongside which oversampling is removed (default: 1)
"""
function remove_oversampling!(raw::RawAcquisitionData; dims = [1])

    dimension_index = dims[1]
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
    apply_girf!(raw::RawAcquisitionData, freq::AbstractVector, g_data::AbstractMatrix)
Applies the GIRF to the trajectories inside of a::AcquisitionData

# Arguments
* `a::AcquisitionData{T}`          - AcquisitionData object
* `freq::AbstractVector`           - Vector containing frequencies of GIRF data
* `g_data::AbstractMatrix`         - Matrix of size N x length(freq) containing complex GIRF data
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
    apply_k0!(raw::RawAcquisitionData, freq::AbstractVector, g_data::AbstractMatrix)
Applies the K0 modulation due to imaging gradients to the data inside of a::AcquisitionData

# Arguments
* `a::AcquisitionData{T}`          - AcquisitionData object
* `freq::AbstractVector`           - Vector containing frequencies of GIRF data
* `k0_data::AbstractMatrix`         - Matrix of size N x length(freq) containing complex k0 function data
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
    save_map(filename, calib_map, resolution_mm; offset_mm = [0.0, 0.0, 0.0])
Saves calibration maps (sensitivity or B0) as 4D NIfTI file(s)

For complex-valued data, magnitude and phase can be split into separate files
# Arguments
* `filename::String`            - string filename with extension .nii, example "sensemap.nii"
* `calib_map`                   - [nX nY nZ {nChannels}] 4-D sensitivity or 3D B0 map array 
* `resolution_mm`               - resolution in mm, 3 element vector, e.g., [1.0, 1.0, 2.0]
* `offset_mm`                   - isocenter offset in mm, default: [0.0, 0.0, 0.0]
* `do_split_phase::Bool=false`    - if true, data is saved in two nifti files with suffix "_magn" and "_phase", respectively
                                  to enable display in typical NIfTI viewers
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
    load_map(filename, calib_map, resolution_mm; offset_mm = [0.0, 0.0, 0.0])
Saves calibration maps (sensitivity or B0) as 4D NIfTI file(s)

For complex-valued data, magnitude and phase can be split into separate files
# Arguments
* `filename::String`            - string filename with extension .nii, example "sensemap.nii"
* `do_split_phase::Bool=false`    - if true, data is saved in two nifti files with suffix "_magn" and "_phase", respectively
                                  to enable display in typical NIfTI viewers
# Output
* `calib_map`                    - [nX nY nZ {nChannels}] 4-D sensitivity or 3D B0 map array 
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
Shifts image to different location by applying phase ramp to image

Perhaps this should be called shift_fov
# Arguments
* `acqdata::AcquisitionData{T}`          - AcquisitionData object
* `shift::AbstractVector`           - Vector containing shift
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