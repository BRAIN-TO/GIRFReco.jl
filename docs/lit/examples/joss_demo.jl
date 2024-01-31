#-----------------------------------------------------------------------------------
# # [GIRFReco.jl Example Script](@id example_script)
#-----------------------------------------------------------------------------------

#=
This page demonstrates an example script for using GIRFReco.jl

This page was generated from the following Julia file: [`joss_demo.jl`](@__REPO_ROOT_URL__/docs/lit/examples/joss_demo.jl)

The configuration file is [`recon_config_joss_demo.jl`](@__REPO_ROOT_URL__/docs/lit/examples/recon_config_joss_demo.jl)
=#

#=
## 1. Setup

The necessary Julia packages needed for spiral reconstruction.
=#

#Our developed packages
using GIRFReco, MRIGradients

#MRIReco and its sub-packages
using MRIReco, FileIO, MRIFiles, MRIBase, MRICoilSensitivities

using RegularizedLeastSquares, Flux

using ImageTransformations

using PlotlyJS, Plots

#=
## 2. Configurations for reconstruction

The following file, [`recon_config_joss_demo.jl`](@__REPO_ROOT_URL__/docs/lit/examples/recon_config_joss_demo.jl),
includes general configuration for spiral reconstruction.
It is necessary to execute this file to make sure all parameters are loaded.
Sample Data that works with this script can be found [here](https://doi.org/10.5281/zenodo.7779044).
Please download, extract and set the `root_project_path` as the top level folder (should be something like `/your/path/joss_data_zenodo/`)
=#

root_project_path = "/your/path/joss_data_zenodo/" # Root path of the data extracted from Zenodo
include("recon_config_joss_demo.jl")

plotlyjs()

# Two parameters determining whether we want to reload the raw data of spiral and GIRF.
reload_spiral_data = true; # Set true if we need to reload raw data compulsively.
reload_girf_data = true; # Set true if we need to reload GIRF data compulsively.

#=
Choose Slice ([single number] OR [1, 2, 31, ...]）
Leave an empty array ([]) or remove this line to select all slices
=#
slice_choice = [];

#=
Choose which diffusion directions and averages to be processed. 
Diffusion direction index starts from 0 (b=0) to the total number in MDDW protocol (e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs). 
Index for average starts from 1.
=#
diffusion_direction = 0
idx_average = 1
num_total_diffusion_directions = params_general[:num_total_diffusion_directions]

## Determine to reconstruct single-interleave data, or one interleave out of multi-interleave data.
is_single_interleave = isa(params_general[:scan_fullpath], String)

#=
Choose which interleave to be reconstructed. 

For multi-interleave data, the range of this value is [1:num\_total\_interleaves].

For single-interleave data, it should always be set as 1; for multi-interleave data, the value set here will be used, indicating which interleaves to be merged and reconstructed.
=#
start_idx_interleave = 1;

#===================================================
## 3. Image Reconstruction

The steps of image reconstruction starts here.

### 3.1 Calculation of B₀ and Coil Sensitivity Maps

The first step in reconstruction pipeline is to calculate the off-resonance (B₀) maps `b0_map` 
and coil sensitivity maps `cartesian_sensitivity` through the Cartesian reconstruction script 
[cartesian_recon.jl](@__REPO_ROOT_URL__/recon/cartesian_recon.jl).

Ideally this script is execute once and the calculated maps are 
saved into files, which are loaded for future usage to save calculation time. 
This is controlled by `do_load_maps` in general parameters. 
===================================================#

if params_general[:do_load_maps] && isfile(params_general[:b0_map_save_fullpath])
    @info "Loading SENSE and B₀ maps from $(params_general[:sensitivity_save_fullpath]) and $(params_general[:b0_map_save_fullpath])"
    cartesian_sensitivity = load_map(params_general[:sensitivity_save_fullpath]; do_split_phase = true)
    b0_maps = load_map(params_general[:b0_map_save_fullpath])
    num_slices = size(b0_maps, 3)
else
    @info "Running cartesian_recon to retrieve maps (cartesian_sensitivity and b0_maps)"
    include("../../../recon/cartesian_recon.jl")
    num_slices = size(b0_maps, 3)
end

#=
### 3.2 Preparation of Spiral Reconstruction

With off-resonance (B₀) maps and coil sensitivity maps calculated, 
before the reconstruction of spiral images, there are necessary steps to prepare for 
the related data. 

#### 3.2.1 Data Selection

The first step is to select the part of spiral k-space data that we 
would like to reconstruct. This include selecting slices, diffusion directions, 
and averages that we want.

First we sort the slice index that we selected to reconstruct.
=#

if isempty(slice_choice) || !(@isdefined slice_choice)
    slice_choice = collect(1:num_slices)
end

is_multislice = length(slice_choice) > 1

if !is_multislice
    selected_slice = slice_choice
else
    selected_slice = sort(vec(slice_choice))
end

#=
Next we select the data we would like to reconstruct from the ISMRMRD file. 

The ISMRMRD data are stored in the following loops:

Slice 1, Slice 2 ... Slice N   Slice 1, Slice 2 ... Slice N     Slice 1, Slice 2 ... Slice N ... 

|______ Diff Dir 1 ______|   |______ Diff Dir 2 ______| ... |______ Diff Dir N ______| ... 

|_________________________________ Average 1 ___________________________________| ... |___ Average N___| 

Here we chose the set corresponding to the b-value = 0 images under the first average as the example.

Note that (1) The raw data file begins with a series of pre-scan profiles with a length of `num_slices*2` 
and we want to skip them; (2) There is a B1 measurement data profile between each k-space readout profile 
which also need to be skipped. Thus the reading of data profiles starts from `num_slices*2 + 2` with 
an increment of `2`.
=#
excitation_list = collect(num_slices*2+2:2:num_slices*4) .+ diffusion_direction * num_slices * 2 .+ (idx_average - 1) * num_slices * (num_total_diffusion_directions + 1) * 2
slice_selection = excitation_list[selected_slice]

#=
#### 3.2.2 Synchronizing and Merging of k-space Data and Trajectory

Since the k-space data and spiral k-space trajectories are sampled under different sampling rates 
and stored in separate files, they need to be first synchronized into the frequency of k-space data 
and then merged into a single object before final spiral image reconstruction.

Here we use a dictionary `params_spiral` to hold the parameters for this k-space data/trajectory synchronization and merging.
=#

params_spiral = Dict{Symbol,Any}()
params_spiral[:recon_size] = Tuple(params_general[:recon_size])
params_spiral[:interleave] = start_idx_interleave
params_spiral[:num_samples] = params_general[:num_adc_samples]
params_spiral[:delay] = 0.00000
params_spiral[:interleave_data_filenames] = params_general[:scan_fullpath]
params_spiral[:traj_filename] = params_general[:gradient_fullpath]
params_spiral[:excitations] = slice_selection
params_spiral[:do_multi_interleave] = !is_single_interleave
params_spiral[:do_odd_interleave] = false
params_spiral[:num_interleaves] = is_single_interleave ? 1 : length(params_spiral[:interleave_data_filenames]) # one interleaf per file, count files, if filenames are array of strings (not only one string)
params_spiral[:single_slice] = !is_multislice

#=
Here we synchronize the spiral k-space data with trajectory by upsampling the trajectory. 
Subsequently, data of all the selected spiral interleaves and the corresponding trajectories 
are merged into `imaging_acq_data`. 
This step is done through the function `merge_raw_interleaves`, which can be viewed in 
[utils.jl](@__REPO_ROOT_URL__/src/utils/utils.jl).

Since the loaded/calculated sens maps and B₀ maps are in ascending slice order,
they need to be reordered according to the slice order in the spiral RawAcqData.

We only do these steps when they have not been done yet or it's specifically required.
=#
if reload_spiral_data || !(@isdefined imaging_acq_data) || !(@isdefined slice_idx_array_spiral)
    @info "Reading spiral data and merging interleaves"
    imaging_acq_data = merge_raw_interleaves(params_spiral, false)
    raw_temp  = RawAcquisitionData(ISMRMRDFile(params_general[:scan_fullpath][1]))
    slice_idx_array_spiral = get_slice_order(raw_temp, num_slices, (num_slices+1)*2, 2)
    b0_maps = b0_maps[:, :, invperm(slice_idx_array_spiral)]
    cartesian_sensitivity = cartesian_sensitivity[:, :, invperm(slice_idx_array_spiral), :]
end

#=
#### 3.2.3 Correction of k-space Trajectory Using Gradient Impulse Response Function

The previously calculated GIRFs are loaded. 

The spiral trajectory is corrected by the 1st and 0th order of GIRF.

=#

#Correct trajectory with the first order GIRFs (K1)
girf_k1 = readGIRFFile(params_general[:girf_fullpath][1], params_general[:girf_fullpath][2], params_general[:girf_fullpath][3], "GIRF_FT", false)
girf_applier_k1 = GirfApplier(girf_k1, params_general[:gamma])

#Correct trajectory with the zeroth order GIRFs (K0)
girf_k0 = readGIRFFile(params_general[:girf_fullpath][1], params_general[:girf_fullpath][2], params_general[:girf_fullpath][3], "b0ec_FT", true)
girf_applier_k0 = GirfApplier(girf_k0, params_general[:gamma])

if params_general[:do_correct_with_girf_k1]
    @info "Correcting For GIRF"
    apply_girf!(imaging_acq_data, girf_applier_k1)
end

if params_general[:do_correct_with_girf_k0]
    @info "Correcting For k₀"
    apply_k0!(imaging_acq_data, girf_applier_k0)
end

# Check if the k-space trajectory is normalized to the range `[-0.5, 0.5]`.
check_acquisition_nodes!(imaging_acq_data)

#=
#### 3.2.4 Center the Object to the Field-of-View (FOV)

If the scanned object is not in the center of the FOV, we need to shift FOV 
to place the object in the center. This is achieved through adding linear phases 
on all dimensions.
=#
shift_kspace!(imaging_acq_data, params_general[:fov_shift])

#=
#### 3.2.5 Processing Coil Sensitivity Maps

We need to preprocess the coil sensitivity maps before reconstruction. 
This includes resizing the coil maps to the size of output encoding matrix size; 
compress the channels according to user's setting to achieve a faster reconstruction.
=#
sensitivity = mapslices(x -> imresize(x, params_spiral[:recon_size][1], params_spiral[:recon_size][2]), cartesian_sensitivity, dims = [1, 2])

# Optional: Plot the sensitivity maps of each coil on a given slice.
if params_general[:do_plot_recon]
    plotlyjs()
    plot_sense_maps(sensitivity,20)
end

# Optional: Coil compression to further reduce the time of recon
if params_general[:do_coil_compression]
    imaging_acq_data, ccMat_vec = softwareCoilCompression(imaging_acq_data, params_general[:num_virtual_coils])
    sensitivity = applyCoilCompressionSensitivityMaps(sensitivity, ccMat_vec)
end

#=
#### 3.2.6 Processing Off-Resonance (B₀) Maps

We need to resize the B₀ maps to the size of output encoding matrix size.

=#
resized_b0_maps = mapslices(x -> imresize(x, params_spiral[:recon_size][1], params_spiral[:recon_size][2]), b0_maps, dims = [1, 2])

#=
#### 3.2.7 Alignment of Off-Resonance, Sensitivity, and Spiral Data

We need to make sure that the axes line up so we rotate the sensitivities and the off-resonance maps.

Depending on your geometry, this might not be necessary but in case you need them:

```
resized_b0_maps = mapslices(x->rotl90(x),resized_b0_maps,dims=[1,2])
sensitivity = mapslices(x->rotl90(x),sensitivity,dims=[1,2])
```
=#

#=
### 3.3 Spiral Image Reconstruction

Here we start the spiral image reconstruction.

First we need to set necessary parameters for reconstruction, 
including iterative solver's setting, coil maps, B₀ maps, etc. 
These parameters are held under the dictionary `params_recon`.

Note that it is safer to cast B₀ maps to ComplexF32 if the current version of MRIReco.jl is used.
=#

@info "Setting Reconstruction Parameters"
params_recon = Dict{Symbol,Any}()
params_recon[:reco] = "multiCoil"
params_recon[:reconSize] = params_spiral[:recon_size][1:2] # cannot avoid camel-case here since it is defined by MRIReco.jl and RegularizedLeastSquares.jl
params_recon[:regularization] = "L2"
params_recon[:λ] = 1e-3
params_recon[:iterations] = params_general[:num_recon_iterations]
params_recon[:solver] = "cgnr"
params_recon[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false)
params_recon[:senseMaps] = ComplexF32.(sensitivity[:, :, selected_slice, :]) # cannot avoid camel-case here since it is defined by MRIReco.jl and RegularizedLeastSquares.jl

if params_general[:do_correct_with_b0_map]
    params_recon[:correctionMap] = ComplexF32.(-1im .* resized_b0_maps[:, :, selected_slice]) # cannot avoid camel-case here since it is defined by MRIReco.jl and RegularizedLeastSquares.jl
end

#= 
Finally we can call reconstruction function of the package `MRIReco.jl` 
to perform final spiral image reconstruction.
=#
@info "Performing Spiral Reconstruction"
@time reco = reconstruction(imaging_acq_data, params_recon)

GC.gc() # Recommended to force triger garbage collection especially when encountering memory issues.

# Reorder slices of the reconstructed images to an ascending order
reco = reco[:,:,slice_idx_array_spiral[selected_slice]]
resized_b0_maps = resized_b0_maps[:, :, slice_idx_array_spiral[selected_slice]]

#=
## 4. Save and Plot the Results (Optional)

All results could be saved into NIfTI files using the `save_map` function 
and be plotted using the `plot_reconstruction` function, both located in 
the file [utils.jl](@__REPO_ROOT_URL__/src/utils/utils.jl).

=#

if params_general[:do_save_recon]
    resolution_tmp = fieldOfView(imaging_acq_data)[1:2] ./ encodingSize(imaging_acq_data)
    resolution_mm = (resolution_tmp[1], resolution_tmp[2], fieldOfView(imaging_acq_data)[3] * (1 + params_general[:slice_distance_factor_percent] / 100.0)) #for 2D only, since FOV[3] is slice thickness then, but gap has to be observed
    num_slices = numSlices(imaging_acq_data)
    save_map(
        params_general[:recon_save_fullpath],
        params_general[:saving_scalefactor] * reco.data[:, :, slice_idx_array_spiral],
        resolution_mm;
        do_split_phase = true,
        do_normalize = params_general[:do_normalize_recon],
    )
end

if params_general[:do_plot_recon]
    @info "Plotting Reconstruction"
    plotlyjs()
    plot_reconstruction(
        reco,
        1:length(selected_slice),
        resized_b0_maps,
        fig_handles = ["Original Magnitude", "Original Phase", "B0"],
        is_slice_interleaved = false,
        rotation = 90,
    )
end

@info "Successfully Completed Spiral Recon"
