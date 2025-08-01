using HDF5, MRIReco, LinearAlgebra, DSP, FourierTools, ROMEO, MRIGradients

include("../utils/utils.jl")
include("../utils/fieldmap_estimator.jl")

## Dictionary of frequently changed parameters
include("recon_config.jl")

## Load data files

# Echo times for field map raw data, in ms
TE1 = params_general[:mapTEs_ms][1]
TE2 = params_general[:mapTEs_ms][2]

@info "Loading Data Files"

b0_filename = params_general[:map_scan_fullpath];

# filename for preprocessed data (remove oversampling, permute dimensions wrt MRIReco)
processed_filename = params_general[:processed_map_scan_fullpath]

if params_general[:do_process_map_scan]

    # Set the data file name (Change this for your own system)
    cartesian_data_file = ISMRMRDFile(b0_filename)

    # read in the raw data from the ISMRMRD file into a RawAcquisitionData object
    r = RawAcquisitionData(cartesian_data_file)

    # does not change anything...
    # r.params["reconFOV"] = [230, 230, 2]

    # Preprocess Data and save!
    preprocess_cartesian_data(r::RawAcquisitionData, true; fname = processed_filename)

end

# Load preprocessed data!
new_data_file = ISMRMRDFile(processed_filename)


new_raw_data = RawAcquisitionData(new_data_file)
cartesian_acq_data = AcquisitionData(new_raw_data, estimateProfileCenter = true)

# Define coils and slices
num_coils = size(cartesian_acq_data.kdata[1], 2)
num_slices = numSlices(cartesian_acq_data)


slice_index_array = get_slice_order(num_slices, isSliceInterleaved = true)
# shift FOV to middle :) 
#TODO: in MRIReco v0.7, try: correctOffset(cartesian_acq_data, [0 -20 0])
shift_kspace!(cartesian_acq_data, params_general[:fov_shift])
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
sensitivity = cartesian_sensitivity

resolution_mm = fieldOfView(cartesian_acq_data) ./ size(sensitivity)[1:3]
resolution_mm[3] = fieldOfView(cartesian_acq_data)[3] * (1 + params_general[:slice_distance_factor_percent] / 100.0); # for 2D only, since FOV[3] is slice thickness then, but gap has to be observed


# save SENSE maps
if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    # TODO: use correct slice order everywhere, e.g., when saving/loading maps for spiral recon
    save_map(params_general[:sensitivity_save_fullpath], sensitivity[:, :, slice_index_array, :], resolution_mm; doSplitPhase = true)
end

plot_sense_maps(sensitivity, num_coils)

#cartesian_acq_data.traj[1].cartesian = false
#cartesian_acq_data.traj[2].cartesian = false
## Parameter dictionary definition for reconstruction

@info "Setting Parameters"
params_cartesian = Dict{Symbol,Any}() # instantiate dictionary
params_cartesian[:reco] = "multiCoil" # choose multicoil reconstruction

# TODO: make recon size and FOV variable!
params_cartesian[:reconSize] = (cartesian_acq_data.encodingSize[1], cartesian_acq_data.encodingSize[2]) # set recon size to be the same as encoded size
params_cartesian[:regularization] = "L2" # choose regularization for the recon algorithm
params_cartesian[:λ] = 1.e-2 # recon parameter (there may be more of these, need to dig into code to identify them for solvers other than cgnr)
params_cartesian[:iterations] = params_general[:num_recon_iterations] # number of CG iterations
params_cartesian[:solver] = "cgnr" # inverse problem solver method
params_cartesian[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false) # turn on store solutions if you want to see the reconstruction convergence (uses more memory)
params_cartesian[:senseMaps] = ComplexF32.(sensitivity) # set sensitivity map array


## Call the reconstruction function

@info "Performing Reconstruction"
@time cartesian_reco = reconstruction(cartesian_acq_data, params_cartesian)

# save Map recon (multi-echo etc.)
if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    save_map(params_general[:map_save_fullpath], cartesian_reco.data[:, :, slice_index_array, :, :, :], resolution_mm; doSplitPhase = true)
end

## Calculate B0 maps from the acquired images (if two TEs)

slices = 1:length(slice_index_array)

@info "Calculating B0 Maps"
# b0_maps = calculate_b0_maps(cartesian_reco.data,slices, TE1, TE2)
b0_maps = estimate_b0_maps(cartesian_reco.data, slices, TE1, TE2, true; β = params_general[:b0_map_beta], reltol = 1e-4)

# save B0 map
if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    save_map(params_general[:b0_map_save_fullpath], b0_maps[:, :, slice_index_array], resolution_mm; doNormalize = false) # no normalization, we want absolute values for offres maps
end


if params_general[:do_plot_recon]
    @info "Plotting Cartesian Results (Sensitivity Maps and B0 Maps)"
    pygui(true) # Leave this code till we need plotting.
    # plot_sense_maps(sensitivity,num_coils)
    plot_reconstruction(cartesian_reco[:, :, :, 1], 1:size(cartesian_reco, 3), b0_maps, isSliceInterleaved = true, rotateAngle = 270)
end

# cleanup unused file
if !params_general[:do_save_processed_map_scan] && params_general[:do_process_map_scan]
    rm(processed_filename)
end

@info "Successfully Completed cartesian_reconstruction"

