using HDF5, MRIReco, LinearAlgebra, Dierckx, DSP, FourierTools, ImageBinarization, ImageEdgeDetection, MRIGradients

# %%
# Include tools and reader functions for running the spiral reconstruction recipe
# Note: the files are found relative of the location of the folder, not the
# environment current folder
include("../io/gradient_reader.jl")
include("../utils/utils.jl")

## Executing Cartesian recon from which B0/sensitivity maps have been computed
@info "Running cartesian_recon to retrieve maps (cartesian_sensitivity and b0_maps)"
include("../recon/cartesian_recon.jl")

## Set figures to be unlocked from the window (i.e use matplotlib backend with controls)

## Choose Slice (can be [single number] OR [1,2,3,4,5,6,7,8,9]
slice_choice = [1, 2, 3, 4, 5, 6, 7, 8, 9] # UNCOMMENT FOR MULTISLICE
# slice_choice = [7] # UNCOMMENT FOR SINGLESLICE (SLICES 3, 7 and 8 are good examples)
diffusion_direction = 0 # CAN BE FROM 0 (b=0) to 6 (1-6 are 6 directions of b=1000)

## Spiral Reconstruction Recipe Starts Here
@info "Starting Spiral Reconstruction Pipeline"

## Default to single slice selection. Choose multi-slice only if computer is capable.
multi_slice = false

if length(slice_choice) > 1
    multi_slice = true
end

if !multi_slice
    selected_slice = slice_choice
else
    selected_slice = sort(vec(slice_choice))
end

## The ISMRMRD File contains more than one excitation, so we choose the set corresponding to the b-value 0 images
excitation_list = vec(20:2:36) .+ diffusion_direction * 18 # DATASET SPECIFIC INDEXING
slice_selection = excitation_list[selected_slice]

@info "Slice Chosen = $selected_slice: \n \nExcitations Chosen = $excitation_list "

# params_general is the dictionary that sets the information for correct data loading and trajectory and data synchronization
params_general = Dict{Symbol,Any}()
params_general[:recon_size] = (200, 200)
params_general[:interleave] = 1
params_general[:slices] = 1
params_general[:coils] = 40
params_general[:num_samples] = 15475
params_general[:delay] = 0.00000 # naive delay correction

params_general[:interleave_data_filenames] = ["data/Spirals/523_96_2.h5", "data/Spirals/523_98_2.h5", "data/Spirals/523_100_2.h5", "data/Spirals/523_102_2.h5"]
params_general[:traj_filename] = "data/Gradients/gradients523.txt"
params_general[:excitations] = slice_selection

params_general[:do_multi_interleave] = true
params_general[:do_odd_interleave] = true
params_general[:num_interleaves] = 4

params_general[:single_slice] = !multi_slice

@info "Using Parameters:\n\nrecon size = $(params_general[:recon_size]) \n interleave = $(params_general[:interleave]) \n slices = $(params_general[:slices]) \n coils = $(params_general[:coils]) \n numSamples = $(params_general[:num_samples])\n\n"
# define recon size and parameters for data loading

## Convert raw to AcquisitionData

@info "Merging interleaves and reading data \n"
acq_data_imaging = merge_raw_interleaves(params_general)

@info "Loading Gradient Impulse Response Functions \n"
## Load GIRFs!
girf_k1 = loadGirf(1, 1)
girf_applier_k1 = GirfApplier(girf_k1, 42577478)

@info "Correcting For GIRF \n"
apply_girf!(acq_data_imaging, girf_applier_k1)

# Load K₀ GIRF
girf_k0 = loadGirf(0, 1)
girf_applier_k0 = GirfApplier(girf_k0, 42577478)

@info "Correcting For k₀ \n"
apply_k0!(acq_data_imaging, girf_applier_k0)

## Check the k-space nodes so they don't exceed frequency limits [-0.5, 0.5] (inclusive)
check_acquisition_nodes!(acq_data_imaging)

do_coil_compression = false

num_virtual_coils = size(sensitivity, 4)

## Do coil compression to make recon faster
if do_coil_compression
    num_virtual_coils = 8
    acq_data_imaging, cartesian_sensitivity = geometricCC_2d(acq_data_imaging, cartesian_sensitivity, num_virtual_coils)
end

## Sense Map loading
@info "Recalculating Sense Maps \n"
cartesian_sensitivity = espirit(
    cartesian_acq_data,
    (4, 4),
    12,
    (acq_data_imaging.encodingSize[1], acq_data_imaging.encodingSize[2]),
    eigThresh_1 = 0.01,
    eigThresh_2 = 0.98,
    match_acq_size = false,
)
sensitivity = mapslices(rotl90, cartesian_sensitivity, dims = [1, 2])

# ## Plot the sensitivity maps of each coil
@info "Plotting SENSE Maps \n"
plot_sense_maps(sensitivity, num_virtual_coils, sliceIndex = 6)

## B0 Maps (Assumes we have a B0 map from gradient echo scan named b0)
@info "Validating B0 Maps \n"
resized_b0_maps = mapslices(x -> imresize(x, (acq_data_imaging.encodingSize[1], acq_data_imaging.encodingSize[2])), b0_maps, dims = [1, 2])

## Define Parameter Dictionary for use with reconstruction
# CAST TO ComplexF32 if you're using current MRIReco.jl

@info "Setting Parameters \n"
params = Dict{Symbol,Any}()
params[:reco] = "multiCoil"
params[:reconSize] = params_general[:recon_size]
params[:regularization] = "L2"
params[:λ] = 1e-2 # CHANGE THIS TO GET BETTER OR WORSE RECONSTRUCTION RESULTS
params[:iterations] = 20
params[:solver] = "cgnr"
params[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false)
params[:senseMaps] = ComplexF32.(sensitivity[:, :, selected_slice, :])
params[:correctionMap] = ComplexF32.(-1im .* resized_b0_maps[:, :, selected_slice])

## Call to reconstruction
@info "Performing Reconstruction \n"
@time reco = reconstruction(acq_data_imaging, params)

#totalRecon = sum(abs2,reco.data,dims=5)
@info "Plotting Reconstruction \n"
plot_reconstruction(reco, 1:length(selected_slice), resized_b0_maps[:, :, selected_slice])

## Plot the image edges (feature comparison)

# img_edges₁ = detect_edges(slice1,Canny(spatial_scale = 2.6))
# img_edges₂ = detect_edges(slice2,Canny(spatial_scale = 2.7))

# imEdges = cat(img_edges₁,img_edges₂,zeros(size(img_edges₁)),dims=3)

# figure("Edge Differences")
# PyPlot.imshow(imEdges)

@info "Successfully Completed Spiral Reconstruction \n"
