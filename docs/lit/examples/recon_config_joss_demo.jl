#-----------------------------------------------------------------------------------
# # [GIRFReco.jl Example Configuration File](@id example_config)
#-----------------------------------------------------------------------------------

#=
## Introduction

This recon_config.jl file describes all general reconstruction parameters, 
as well as data locations and selections for an iterative non-Cartesian reconstruction that relies 
on an external reference scan (Cartesian) to estimate calibration maps (coil sensitivities, B0 maps).

All parameters are stored in a dictionary named `params_general`.

=#
params_general = Dict{Symbol,Any}();

# Gyromagnetic ratio, in unit of Hz
params_general[:gamma] = 42577478;

# General options for reconstruction script:
params_general[:do_load_maps] = true              # if true, reloads B0/SENSE maps instead of recalculating
params_general[:do_save_recon] = true              # if true, saves reconstruction and all auxiliary image data (maps) as NIfTI files
params_general[:do_plot_recon] = true              # if true, plots intermediate debugging and output recon figures (needs graphics, not recommended in multi-thread mode due to PyPlot)
params_general[:do_process_map_scan] = true         # if true, compute sensitivity and B0 maps from reconstructed Cartesian scan   
params_general[:do_save_processed_map_scan] = false  # save ISMRMD file of preprocessed Cartesian data (before recon)
params_general[:recon_id] = "v7"               # unique identifier for the saved result files
params_general[:do_correct_with_b0_map] = true       # whether perform off-resonance correction
params_general[:do_correct_with_girf_k1] = true    # whether perform 1st order GIRF correction
params_general[:do_correct_with_girf_k0] = true    # whether perform 1st order GIRF correction
params_general[:do_coil_compression] = false       # whether perform coil compression
params_general[:do_normalize_recon] = false       # if true, set the range of magnitude image as [0 1]

# General parameters for reconstruction:
params_general[:num_adc_samples] = 15650             # total number of ADC points BEFORE the rewinder at the end of the spiral readout. Need to check with data. # 15504 for 523, 15655 for gradient 508
params_general[:recon_size] = [200, 200, 1]          # the matrix size of the reconstructed images. Needs to specify 1 on Z dimension for 2D images
params_general[:num_recon_iterations] = 40             # number of recon iterations (for both Cartesian and Spiral recon)
params_general[:b0_map_beta] = 0.1             # for estimate_b0_maps, * `β` - Regularization parameter controlling roughness penalty (larger = smoother, default 5e-4)
params_general[:saving_scalefactor] = 1.0e9    # typical range of recon intensities is 1e-7, rescale when saving, e.g., to 0...1000 roughly for fMRI analysis
params_general[:num_virtual_coils] = 8                 # if perform coil compression, the number of coils to be compressed to
params_general[:fov_shift] = [-2, -5]               # amount of FOV shift; in unit of number of voxels in [x,y] direction (for our in-vivo dataset this is [0, -15], change as necessary)
params_general[:slice_distance_factor_percent] = 400 # Scan parameters, Additional acquisition information, e.g., slice distance etc.
params_general[:num_total_diffusion_directions] = 6          # Need to specify total diffusion directions included in the raw data
#=

## Data selector

See details in the page Advanced Usage.

Data selector is designed for calling the script repetitively through [`RunReconLoop.jl`](@__REPO_ROOT_URL__/recon/RunReconLoop.jl) 
for e.g. different diffusion directions and/or different averages. 
Note that the index of diffusion direction starts from 0 (b=0) to the total number in MDDW protocol, 
e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs. 
`boolean is_called_from_global_recon` should be set as true if this `RunReconLoop.jl` is active.
If is_called_from_global_recon is false or not defined, the data selector needs to be defined in the following code block.
=#
if !(@isdefined is_called_from_global_recon) || !is_called_from_global_recon
    global selector = Dict{Symbol,Any}()
    selector[:avg] = 1
    selector[:seg] = 1
    selector[:dif] = 0
end


#=
## Directories and File Names

We set directories for data reading in and results writing in this section. 

### Specifying Directories
=#

params_general[:project_path] = root_project_path # Root path for the project

#src # Paths (user-dependent)
#src # params_general[:data_path] = joinpath(params_general[:project_path], "data", "joss_demo", "Human", "dat")
#src # params_general[:gradients_path] = joinpath(params_general[:project_path], "data", "joss_demo", "gradients")
#src # params_general[:results_path] = joinpath(params_general[:project_path], "results", "joss_demo", "Human")
#src # params_general[:girf_path] = joinpath(params_general[:project_path], "code", "GIRFReco", "data", "GIRF", "GIRF_ISMRM2022")
#src # params_general[:recon_save_path] = joinpath(params_general[:results_path], "recon", params_general[:recon_id])

#Path to ISMRMRD files (raw k-space data) [Input]
params_general[:data_path] = params_general[:project_path]
#Path to spiral readout gradient files [Input]
params_general[:gradients_path] = joinpath(params_general[:project_path], "Gradients")
#Path to GIRF files [Input]
params_general[:girf_path] = joinpath(params_general[:project_path], "GIRF", "GIRF_ISMRM2022")
#Path to middle results (coil and B0 maps) files [Output]
params_general[:results_path] = joinpath(params_general[:project_path], "results", "phantom")
#Path to final reconstructed spiral images [Output]
params_general[:recon_save_path] = joinpath(params_general[:results_path], "recon", params_general[:recon_id]);

#=
### Specifying File Names
=#
params_general[:map_scan_filename] = "Fieldmaps/meas_MID00083_FID06181_GRE_FieldMap_DualEcho_2mm.mrd" # Cartesian dual-echo file, for coil and B0 maps calculation [Input]
params_general[:map_scan_filename_stem] = "meas_MID00083_FID06181_GRE_FieldMap_DualEcho_2mm.mrd"
# params_general[:mapTEs_ms] = [4.92,  7.38] # Two echo times, in ms
params_general[:girf_filename] = ["2021Nov_PosNeg_Gx.mat", "2021Nov_PosNeg_Gy.mat", "2021Nov_PosNeg_Gz.mat"] # Calculated GIRF for each gradient axis [Input]
params_general[:gradient_filename] = joinpath("gradients508.txt") # File name for the spiral gradient [Input]
params_general[:scan_filename] = ["Spirals/meas_MID00072_FID06170_diffSpiral_508_Intl0_b2k_4Avg.mrd"]# ISMRMRD Raw k-space data for spiral acquisition [Input]
#Remove the comment mark of the next line if you want a multi-interleave spiral reconstruction
#params_general[:scan_filename] = ["Spirals/meas_MID00072_FID06170_diffSpiral_508_Intl0_b2k_4Avg.mrd","Spirals/meas_MID00074_FID06172_diffSpiral_508_Intl1_b2k_4Avg.mrd","Spirals/meas_MID00076_FID06174_diffSpiral_508_Intl2_b2k_4Avg.mrd","Spirals/meas_MID00078_FID06176_diffSpiral_508_Intl3_b2k_4Avg.mrd"]
params_general[:scan_filename_stem] = "meas_MID00072_FID06170_diffSpiral_508_Intl0_b2k_4Avg.mrd"
params_general[:processed_map_scan_filename] = "processed_cartesian_data.h5" # file name for preprocessed data (remove oversampling, permute dimensions wrt MRIReco) [Output]
params_general[:map_save_filename] = splitext(params_general[:map_scan_filename_stem])[1] * "_reconmap.nii" # File name for reconstructed dual-echo Cartesian images [Output]
params_general[:sensitivity_save_filename] = splitext(params_general[:map_scan_filename_stem])[1] * "_sensemap.nii" # File name for calculated coil sensitivity maps [Output]
params_general[:b0_map_save_filename] = splitext(params_general[:map_scan_filename_stem])[1] * "_b0map.nii"; # File name for calculated off-resonance (B0) maps [Output]

#=
File name for the final reconstructed spiral image.
If we reconstructing multiple spiral data files (e.g. multiple interleaves) through `RunReconLoop.jl`, 
the file name for the final reconstructed image is concatenated from multiple scan file names. 
Otherwise, just append `_recon.nii` as suffix to file name.
=#
if isa(params_general[:scan_filename], AbstractVector)
    params_general[:recon_save_filename] =
        join([(x[1] * "_") for x in splitext.(params_general[:scan_filename_stem])]) *
        "dif$(selector[:dif])_" *
        "itl$(selector[:seg])_" *
        "avg$(selector[:avg])_" *
        "recon.nii"
else
    params_general[:recon_save_filename] = splitext(params_general[:scan_filename_stem])[1] * "_recon.nii"
end

#=
### Assembling Full Paths

Assembling directories and file names for final full pathes. 
These are automated operations.
=#
params_general[:gradient_fullpath] = joinpath(params_general[:gradients_path], params_general[:gradient_filename]) # Full paths of spiral readout gradients
#src # params_general[:pathload_maps] = 
#src # joinpath(params_general[:results_path], "recon", params_general[:recon_id]) # NOTE: If loaded from other recon_id, this path might differ
params_general[:girf_fullpath] = joinpath.(params_general[:girf_path], params_general[:girf_filename]) # Full paths of GIRF files
params_general[:map_scan_fullpath] = joinpath(params_general[:data_path], params_general[:map_scan_filename]) # Full path of dual-echo Cartesian data
params_general[:scan_fullpath] = joinpath.(params_general[:data_path], params_general[:scan_filename]) # Full paths of raw k-space data files of spiral acquisition
params_general[:processed_map_scan_fullpath] = joinpath(params_general[:recon_save_path], params_general[:processed_map_scan_filename]) # Full paths of pre-processed Cartesian dual-echo data [Output]
params_general[:recon_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:recon_save_filename]) # Full paths of the reconstructed spiral image [Output]
params_general[:map_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:map_save_filename]) # Full paths of reconstructed dual-echo Cartesian images [Output]
params_general[:sensitivity_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:sensitivity_save_filename]) # Full paths of calculated coil sensitivity maps [Output]
params_general[:b0_map_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:b0_map_save_filename]); # Full paths of calculated off-resonance (B0) maps [Output]

#=
## Final Steps

If the path for results writing is not existing, create it.

As the last step of configuration, copy this config file 
to the recon path for further checking and debugging purposes.
=#
if ~ispath(params_general[:recon_save_path])
    mkpath(params_general[:recon_save_path])
end

cp(@__FILE__, joinpath(params_general[:recon_save_path], "recon_config.jl"); force = true)



