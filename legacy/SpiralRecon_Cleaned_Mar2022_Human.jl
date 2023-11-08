using HDF5, MRIReco, LinearAlgebra, Dierckx, DSP, FourierTools, ImageBinarization, ImageEdgeDetection, MRIGradients

# %%
# Include tools and reader functions for running the spiral reconstruction recipe
# Note: the files are found relative of the location of the folder, not the
# environment current folder
include("../io/GradientReader.jl")
include("../utils/Utils.jl")

## Set true if we need to reload Cartesian and/or spiral data compulsively.
reload_cartesian_data = true
reload_spiral_data = true
do_coil_compression = true

## Gyromagnetic ratio, in unit of Hz
gamma = 42577478

## Only calculate sensitivity and B0 maps when they have not been done yet, or it's specifically required.
if reload_cartesian_data || !((@isdefined cartesian_sensitivity) && (@isdefined b0_maps))
    ## Executing Cartesian recon from which B0/sensitivity maps have been computed
    @info "Running CartesianRecon to retrieve maps (cartesian_sensitivity and b0_maps)"
    include("CartesianRecon_Mar2022_Human.jl")
end

## Set figures to be unlocked from the win9ow (i.e use matplotlib backend with controls)
#pygui(true)

## Choose Slice (can be [single number] OR [1,2,3,4,5,6,7,8,9]
slice_choice = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # UNCOMMENT FOR MULTISLICE
# slice_choice = [6] # UNCOMMENT FOR SINGLESLICE (SLICES 3, 7 and 8 are good examples)
diffusion_direction = 0 # CAN BE FROM 0 (b=0) to 6 (e.g. for 6 direction MDDW, 1-6 are 6 directions)

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
# excitation_list = vec(20:2:36) .+ diffusion_direction * 9 * 2 # DATASET SPECIFIC INDEXING
excitation_list = vec(32:2:62) .+ diffusion_direction * 15 * 2 # DATASET SPECIFIC INDEXING: 15 slices, starting from profile 32
slice_selection = excitation_list[selected_slice]

@info "Slice Chosen = $selected_slice: \n \nExcitations Chosen = $excitation_list "

fname_spiralIntlv1 = "data/Spirals/508_124_2.h5"
fname_spiralIntlv2 = "data/Spirals/508_126_2.h5"
fname_spiralIntlv3 = "data/Spirals/508_128_2.h5"
fname_spiralIntlv4 = "data/Spirals/508_130_2.h5"
fname_gradient = "data/Gradients/gradients508.txt"
fname_girfGx = "data/GIRF/GIRF_ISMRM2022/2021Nov_PosNeg_Gx.mat"
fname_girfGy = "data/GIRF/GIRF_ISMRM2022/2021Nov_PosNeg_Gy.mat"
fname_girfGz = "data/GIRF/GIRF_ISMRM2022/2021Nov_PosNeg_Gz.mat"

# params_general is the dictionary that sets the information for correct data loading and trajectory and data synchronization
params_general = Dict{Symbol,Any}()
params_general[:recon_size] = (200, 200)
params_general[:interleave] = 1
params_general[:slices] = 1
params_general[:coils] = 20
# params_general[:num_samples] = 16084 # Total Number of ADC event, including the period of gradient rewinder
params_general[:num_samples] = 15504 # Total Number of readout before gradient rewinder
params_general[:delay] = 0.00000 # naive delay correction

params_general[:interleave_data_filenames] = [fname_spiralIntlv1, fname_spiralIntlv2, fname_spiralIntlv3, fname_spiralIntlv4]
params_general[:traj_filename] = fname_gradient
params_general[:excitations] = slice_selection

params_general[:do_multi_interleave] = true
params_general[:do_odd_interleave] = true
params_general[:num_interleaves] = 4

params_general[:single_slice] = !multi_slice

# Defined recon size and parameters for data loading
@info "Using Parameters:\n\nreconSize = $(params_general[:recon_size]) \n interleave = $(params_general[:interleave]) \n slices = $(params_general[:slices]) \n coils = $(params_general[:coils]) \n numSamples = $(params_general[:num_samples])\n\n"

## Only load data when it has not been done yet, or it's specifically required.
if reload_spiral_data || !(@isdefined acq_data_imaging)
    ## Convert raw to AcquisitionData

    @info "Merging interleaves and reading data \n"
    acq_data_imaging = merge_raw_interleaves(params_general)

    @info "Loading Gradient Impulse Response Functions \n"
    ## Load GIRFs!
    # Tim Wu, use new read GIRF function
    #girf_k1 = loadGirf(1,1)
    girf_k1 = readGIRFFile(fname_girfGx, fname_girfGy, fname_girfGz, "GIRF_FT", false)
    girf_applier_k1 = GirfApplier(girf_k1, gamma)

    @info "Correcting For GIRF \n"
    apply_girf!(acq_data_imaging, girf_applier_k1)

    # Load K₀ GIRF
    # Tim Wu, use new read GIRF function
    #girf_k0 = loadGirf(0,1)
    girf_k0 = readGIRFFile(fname_girfGx, fname_girfGy, fname_girfGz, "b0ec_FT", true)
    girf_applier_k0 = GirfApplier(girf_k0, gamma)

    @info "Correcting For k₀ \n"
    apply_k0!(acq_data_imaging, girf_applier_k0)

    ## Check the k-space nodes so they don't exceed frequency limits [-0.5, 0.5] (inclusive)
    check_acquisition_nodes!(acq_data_imaging)

end

## Sense Map loading
@info "Recalculating Sense Maps \n"
sensitivity = espirit(cartesian_acq_data, (4, 4), 12, params_general[:recon_size], eigThresh_1 = 0.01, eigThresh_2 = 0.98)

# shift FOV to middle :) 
shift_kspace!(acq_data_imaging, [0, -20])
#changeFOV!(acq_data_imaging,[1.5,1.5])

num_virtual_coils = size(sensitivity, 4)

do_coil_compression = false

## Do coil compression to make recon faster
if do_coil_compression
    num_virtual_coils = 4
    acq_data_imaging, sensitivity = geometricCC_2d(acq_data_imaging, sensitivity, num_virtual_coils)
end

# ## Plot the sensitivity maps of each coil
if params_general[:do_plot_recon]
    @info "Plotting SENSE Maps"
    plot_sense_maps(sensitivity, num_virtual_coils, sliceIndex = 10)
end

## B0 Maps (Assumes we have a B0 map from gradient echo scan named b0)
@info "Resizing B0 Maps"
resized_b0_maps = mapslices(x -> imresize(x, params_general[:recon_size]), b0_maps, dims = [1, 2])

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
plot_reconstruction(cartesian_reco, 1:length(selected_slice), resized_b0_maps[:, :, selected_slice], isSliceInterleaved = true, rotateAngle = 270)

## Plot the image edges (feature comparison)

# img_edges₁ = detect_edges(slice1,Canny(spatial_scale = 2.6))
# img_edges₂ = detect_edges(slice2,Canny(spatial_scale = 2.7))

# imEdges = cat(img_edges₁,img_edges₂,zeros(size(img_edges₁)),dims=3)

# figure("Edge Differences")
# imshow(imEdges)

@info "Successfully Completed SpiralRecon \n"
