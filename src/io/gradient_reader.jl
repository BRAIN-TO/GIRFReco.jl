

"""
    read_gradient_text_file(filename, reconsize, delay)
Reads in text file containing gradient waveform information

# Arguments
* `filename` - filename (with full path) of text file with gradient waveform information
* `reconsize::Tuple{Int64,Int64,Int64}` - size of reconstructed image (trailing dimension 1 for 2D acquisitions)
* `delay` - delay in seconds from the nominal first sampling point to the actual first sampling point
"""
function read_gradient_text_file(filename, reconsize, delay)

    gradient_data = readdlm(filename, '\n')

    ## Read in the header data of the gradient text file (lines 1 to 21)
    gradient_dict = Dict{Symbol,Any}()
    gradient_dict[:version_number] = gradient_data[1]
    gradient_dict[:num_samples] = gradient_data[2]
    gradient_dict[:dwell_time] = gradient_data[3] # [seconds]
    gradient_dict[:samples_per_interleave] = gradient_data[4]
    gradient_dict[:num_interleaves] = gradient_data[5]
    gradient_dict[:num_dims] = gradient_data[6]
    gradient_dict[:time_to_center_kspace] = gradient_data[7] # [seconds]
    gradient_dict[:acq_duration] = gradient_data[8]
    gradient_dict[:samples_per_acq] = gradient_data[9]
    gradient_dict[:num_acq] = gradient_data[10]
    gradient_dict[:acq_TR] = gradient_data[11]
    gradient_dict[:gradient_acq_start_delay] = gradient_data[12]
    gradient_dict[:echo_time_shift_samples] = gradient_data[13]
    gradient_dict[:fov] = gradient_data[14:16] # [m]
    gradient_dict[:voxel_dims] = gradient_data[17:19] # [m]
    gradient_dict[:gradient_strength_factor] = gradient_data[20] # [mT/m]
    gradient_dict[:is_binary] = gradient_data[21]
    gradient_dict[:gamma] = 42577.478 # [Hz/mT] CAN CHANGE
    gradient_dict[:field_strength] = 3 # [T] CAN CHANGE

    #print(gradient_dict)

    ## reading and data scaling of gradient data
    gradient_dict[:gradient_vector] = gradient_data[22:end]
    interleave_gradient_array =
        gradient_dict[:gradient_strength_factor] * reshape(gradient_dict[:gradient_vector], gradient_dict[:samples_per_interleave], gradient_dict[:num_interleaves], gradient_dict[:num_dims]) #[mT/m]

    planned_times = gradient_dict[:dwell_time] .* (0:(gradient_dict[:samples_per_interleave]-1))
    delayed_times = planned_times .- delay .- gradient_dict[:dwell_time] ./ 2 # seconds (dwell_time/2 compensates for integration)

    interleave_gradient_array_new = Array{Float64,3}(undef, size(interleave_gradient_array))

    #print(size(interleave_gradient_array_new))

    ## Loop over all of the unique excitation trajectories and create an interpolant of the gradient
    for dim = 1:gradient_dict[:num_dims]

        for l = 1:gradient_dict[:num_interleaves]

            #print((dim,l),"\n")

            sp = Spline1D(planned_times, interleave_gradient_array[:, l, dim], w = ones(length(planned_times)), k = 1, bc = "zero", s = 0.0)

            # evaluate the interpolant at the sampling times of the kspace data
            interleave_gradient_array_new[:, l, dim] = sp(delayed_times)

            #print(interleave_gradient_array_new[:,l,dim][end],"\n")

        end

    end

    ## cumulative summation and numerical integration of the gradient data, resulting in the kspace trajectory
    kspace_trajectory_array_new = gradient_dict[:gamma] * gradient_dict[:dwell_time] * cumsum(interleave_gradient_array_new, dims = 1) # [rad/m]

    ## Conversion to the trajectory scaling convention in MRIReco.jl
    #  Currently only 2d Trajectories
    converted_kspace_trajectory_array_new = kspace_trajectory_array_new
    converted_kspace_trajectory_array_new[:, :, 1] *= gradient_dict[:fov][1] ./ reconsize[1]
    converted_kspace_trajectory_array_new[:, :, 2] *= gradient_dict[:fov][2] ./ reconsize[2]
    
    if gradient_dict[:num_dims] == 3
        converted_kspace_trajectory_array_new[:,:,3] *= gradient_dict[:FOV][3] ./ reconsize[3] # Normalized to -0.5 to 0.5, no unit.
    end

    ## Construction of the trajectory object ##

    ## Reshaping of the array to the format expected by the Trajectory constructor in MRIReco.jl
    # - dim 1 = kspace dimension
    # - dim 2 = kspace position (with interleaves/profiles arranged consecutively)
    permuted_trajectory =
        permutedims(reshape(converted_kspace_trajectory_array_new, gradient_dict[:samples_per_interleave] * gradient_dict[:num_interleaves], gradient_dict[:num_dims]), [2, 1])

    ## Construction of the trajectory
    # - Note: timing vectors are automatically generated - seems to be consistent with the dwell time
    trajectory_object = Trajectory(
        permuted_trajectory,
        gradient_dict[:num_interleaves],
        gradient_dict[:samples_per_interleave],
        TE = gradient_dict[:echo_time_shift_samples],
        AQ = gradient_dict[:acq_duration],
        numSlices = 9,
        cartesian = false,
        circular = false,
    )

    return trajectory_object

end