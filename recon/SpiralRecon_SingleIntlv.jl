using HDF5,
    MRIReco,
    LinearAlgebra,
    Dierckx,
    DSP,
    FourierTools,
    ImageBinarization,
    ImageEdgeDetection,
    MRIGradients,
    FileIO,
    MRIFiles,
    MRICoilSensitivities,
    RegularizedLeastSquares,
    GIRFReco,
    MosaicViews,
    Plots,
    Images

# All data-specific recon parameters
include("ReconConfig_SPIDI_0007.jl")

## ----------------------------- User-defined Variables -------------------------- ##

## Set true if we need to reload raw data compulsively.
reload_spiral_data = true
reloadGIRFData = true

# Choose Slice (can be [single number] OR [1,2,3,...])
# Leave empty ([]) or remove this line to later select all slices
slice_choice = []; # TODO: read from ISMRMRD itself

## Gyromagnetic ratio, in unit of Hz
gamma = 42577478

## Choose diffusion direction; starting from 0 (b=0) to the total number in MDDW protocol, e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs)
diffusion_direction = selector[:dif]
idxAverage = selector[:avg]
num_total_diffusion_directions = params_general[:num_total_diffusion_directions] # TODO: Read from ISMRMRD itself

# Which interleave to be reconstructed. For single-interleave data, it will always be set as 1; for multi-interleave data, the value set here will be used.
# For multi-interleaved data, this value is ranging from [1:TotNumIntlv] (total number of interleaves), indicating which interleave to be reconstructed
startIndexIntlv = selector[:seg]

## Determine to reconstruct single-interleave data, or one interleave out of multi-interleave data.
isDataSingleIntlv = isa(params_general[:scan_fullpath], String)


## ------------------------------------------------ Calculation Starts Here ---------------------------------------------------------- ##

if params_general[:do_load_maps] && isfile(params_general[:b0_map_save_fullpath]) # # TODO ask for sense map (but split in magn/phase)
    @info "Loading SENSE and B0 maps from $(params_general[:sensitivity_save_fullpath]) and $(params_general[:b0_map_save_fullpath])"
    # load maps, permute slice, sice files have geometric slice order
    b0_maps = load_map(params_general[:b0_map_save_fullpath])

    num_slices = size(b0_maps, 3)
    slice_index_array = get_slice_order(num_slices, isSliceInterleaved = true)

    b0_maps = b0_maps[:, :, invperm(slice_index_array)]
    cartesian_sensitivity = load_map(params_general[:sensitivity_save_fullpath]; doSplitPhase = true)[:, :, invperm(slice_index_array), :]
else
    ## Only calculate sensitivity and B0 maps when they have not been done yet, or it's specifically required.
    ## Executing Cartesian recon from which B0/sensitivity maps have been computed
    @info "Running CartesianRecon to retrieve maps (cartesian_sensitivity and b0_maps)"
    include("CartesianRecon.jl")
    num_slices = size(b0_maps, 3)
end

## Spiral Reconstruction Recipe Starts Here
@info "Starting Spiral Reconstruction Pipeline"

if isempty(slice_choice) || !(@isdefined slice_choice)
    slice_choice = collect(1:num_slices)
end

isMultiSlice = length(slice_choice) > 1

if !isMultiSlice
    selected_slice = slice_choice
else
    selected_slice = sort(vec(slice_choice))
end

## The ISMRMRD File contains more than one excitation, so we choose the set corresponding to the b-value 0 images
excitation_list = collect(num_slices*2+2:2:num_slices*4) .+ diffusion_direction * num_slices * 2 .+ (idxAverage - 1) * num_slices * (num_total_diffusion_directions + 1) * 2 # DATASET SPECIFIC INDEXING: 15 slices, starting from profile 32
slice_selection = excitation_list[selected_slice]

@info "Slice Chosen = $selected_slice: \n \nExcitations Chosen = $excitation_list "

# params_general is the dictionary that sets the information for correct data loading and trajectory and data synchronization
params_general = Dict{Symbol,Any}()
params_general[:recon_size] = Tuple(params_general[:recon_size])
params_general[:interleave] = startIndexIntlv
params_general[:num_samples] = params_general[:num_adc_samples]
params_general[:delay] = 0.00000 # naive delay correction

params_general[:interleave_data_filenames] = params_general[:scan_fullpath]

params_general[:traj_filename] = params_general[:gradient_fullpath]
params_general[:excitations] = slice_selection

params_general[:do_multi_interleave] = !isDataSingleIntlv
params_general[:do_odd_interleave] = false
params_general[:num_interleaves] = isDataSingleIntlv ? 1 : length(params_general[:interleave_data_filenames]) # one interleaf per file, count files, if filenames are array of strings (not only one string)

params_general[:single_slice] = !isMultiSlice

# Defined recon size and parameters for data loading
@info "Using Parameters:\n\nreconSize = $(params_general[:recon_size]) \n interleave = $(params_general[:interleave]) \n slices = $(slice_choice) \n coils = $(size(cartesian_sensitivity, 4)) \n numSamples = $(params_general[:num_samples])\n\n"

if reloadGIRFData || !(@isdefined girf_k1) || !(@isdefined gAK1) || !(@isdefined girf_k0) || !(@isdefined gAK0)
    @info "Loading Gradient Impulse Response Functions"

    ## Load GIRFs (K1)
    girf_k1 = readGIRFFile(params_general[:girf_fullpath][1], params_general[:girf_fullpath][2], params_general[:girf_fullpath][3], "GIRF_FT", false)
    girf_applier_k1 = GirfApplier(girf_k1, gamma)

    ## Load K₀ GIRF
    girf_k0 = readGIRFFile(params_general[:girf_fullpath][1], params_general[:girf_fullpath][2], params_general[:girf_fullpath][3], "b0ec_FT", true)
    girf_applier_k0 = GirfApplier(girf_k0, gamma)
end

## Only load data when it has not been done yet, or it's specifically required.
if reload_spiral_data || !(@isdefined acq_data_imaging)
    ## Convert raw to AcquisitionData

    @info "Reading spiral data and merging interleaves"
    acq_data_imaging = merge_raw_interleaves(params_general)

    if params_general[:do_correct_with_girf_k1]
        @info "Correcting For GIRF"
        apply_girf!(acq_data_imaging, girf_applier_k1)
    end

    if params_general[:do_correct_with_girf_k0]
        @info "Correcting For k₀"
        apply_k0!(acq_data_imaging, girf_applier_k0)
    end

    ## Check the k-space nodes so they don't exceed frequency limits [-0.5, 0.5] (inclusive)
    check_acquisition_nodes!(acq_data_imaging)

end

## Sense Map loading
@info "Resizing Sense Maps"

# Resize sense maps to match encoding size of data matrix
sensitivity = mapslices(x -> imresize(x, params_general[:recon_size][1], params_general[:recon_size][2]), cartesian_sensitivity, dims = [1, 2])

# Plot the sensitivity maps of each coil
@info "Plotting SENSE Maps"

if params_general[:do_plot_recon]
    plot_sense_maps(sensitivity, size(sensitivity, 4), sliceIndex = 10)
end


# shift FOV to middle :) 
shift_kspace!(acq_data_imaging, params_general[:fov_shift])
# changeFOV!(acq_data_imaging,[0.99, 0.99])


## Do coil compression to make recon faster
if params_general[:do_coil_compression]
    acq_data_imaging, sensitivity = geometricCC_2d(acq_data_imaging, sensitivity, params_general[:num_virtual_coils])
end


## B0 Maps (Assumes we have a B0 map from gradient echo scan named b0)
@info "Resizing B0 Maps"
resized_b0_maps = mapslices(x -> imresize(x, params_general[:recon_size][1], params_general[:recon_size][2]), b0_maps, dims = [1, 2])

## Define Parameter Dictionary for use with reconstruction
# CAST TO ComplexF32 if you're using current MRIReco.jl

@info "Setting Parameters"
params = Dict{Symbol,Any}()
params[:reco] = "multiCoil"
params[:reconSize] = params_general[:recon_size][1:2]
params[:regularization] = "L2"
params[:λ] = 1e-2 # CHANGE THIS TO GET BETTER OR WORSE RECONSTRUCTION RESULTS
params[:iterations] = params_general[:num_recon_iterations]
params[:solver] = "cgnr"
params[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false)
params[:senseMaps] = ComplexF32.(sensitivity[:, :, selected_slice, :])

if params_general[:do_correct_with_b0_map]
    params[:correctionMap] = ComplexF32.(-1im .* resized_b0_maps[:, :, selected_slice])
end

## Call to reconstruction
@info "Performing Reconstruction"
@time reco = reconstruction(acq_data_imaging, params)
#totalRecon = sum(abs2,reco.data,dims=5)

# save Map recon (multi-echo etc.)
if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    resolution_tmp = fieldOfView(acq_data_imaging)[1:2] ./ encodingSize(acq_data_imaging)
    resolution_mm = (resolution_tmp[1], resolution_tmp[2], fieldOfView(acq_data_imaging)[3] * (1 + params_general[:slice_distance_factor_percent] / 100.0)) # for 2D only, since FOV[3] is slice thickness then, but gap has to be observed

    # TODO: use slice ordering from cartesian scan directly!
    num_slices = numSlices(acq_data_imaging)
    slice_index_array = get_slice_order(num_slices, isSliceInterleaved = true)
    save_map(
        params_general[:recon_save_fullpath],
        params_general[:saving_scalefactor] * reco.data[:, :, slice_index_array],
        resolution_mm;
        doSplitPhase = true,
        doNormalize = params_general[:do_normalize_recon],
    )
end

if params_general[:do_plot_recon]
    @info "Plotting Reconstruction"
    #pygui(true)
    plot_reconstruction(
        reco,
        1:length(selected_slice),
        resized_b0_maps[:, :, selected_slice],
        figHandles = ["Original Magnitude", "Original Phase", "B0"],
        isSliceInterleaved = true,
        rotateAngle = 270,
    )
end

@info "Successfully Completed SpiralRecon"
