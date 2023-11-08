## This ReconConfig.jl file describes all reconstruction parameters, as well as data locations and selections for an iterative non-Cartesian reconstruction that relies 
#  on an external reference scan (Cartesian) to estimate calibration maps (coil sensitivities, B0 maps)

using Dates

params_general = Dict{Symbol,Any}()

## General options for recon script
params_general[:do_load_maps] = false               # if true, reloads B0/SENSE maps instead of recalculating
params_general[:do_save_recon] = true              # if true, saves reconstruction and all auxiliary image data (maps) as NIfTI files
params_general[:do_plot_recon] = true             # if true, plots intermediate debugging and output recon figures (needs graphics, not recommended in multi-thread mode due to PyPlot)
params_general[:do_process_map_scan] = true         # if true, compute sensitivity and B0 maps from reconstructed Cartesian scan   
params_general[:do_save_processed_map_scan] = false; # save ISMRMD file of preprocessed Cartesian data (before recon)

## Reconstruction Parameters
# update time stamp for new recon, otherwise keep fixed, will create a new recon/<recon_id> directory
#params_general[:recon_id] = Dates.format(Dates.now(), "yyyy-mm-dd_HH_MM_SS") # recon ID is recon_id
# params_general[:recon_id] = "2022-10-20_09_07_07"
params_general[:recon_id] = "v6";
params_general[:do_correct_with_b0_map] = true
params_general[:do_correct_with_girf_k1] = true
params_general[:do_correct_with_girf_k0] = true

params_general[:num_virtual_coils] = 8;
params_general[:do_coil_compression] = false;
params_general[:fov_shift] = [0, -20]; # Unit: number of voxels

## Scan parameters, Additional acquisition information, e.g., slice distance etc.
params_general[:slice_distance_factor_percent] = 400 # 000

#Total number of ADC points BEFORE the rewinder at the end of the spiral readout. For gradient 508, use 15655 (out of 16084); for gradient 511, use 15445 (out of 15624).
params_general[:num_adc_samples] = 15655 # 15504
# Matrix size of the reconstructed image. For gradient 508 with all 4 interleaves, use 200 for high resolution image; otherwise consider using 112 or 84 for a lower resolution. The FOV is 220 mm for both gradients 508 and 511.
params_general[:recon_size] = [200, 200, 1] #(112, 112) #(200, 200)
params_general[:num_recon_iterations] = 20; # number of recon iterations (for both Cartesian and Spiral recon)
params_general[:b0_map_beta] = 0.1 # for estimate_b0_maps, * `β` - Regularization parameter controlling roughness penalty (larger = smoother, default 5e-4)
params_general[:do_normalize_recon] = false # set max abs to 1
params_general[:saving_scalefactor] = 1.0e9 # 1 # typical range of recon intensities is 1e-7, rescale when saving, e.g., to 0...1000 roughly for fMRI analysis

# Data selector
#  Choose diffusion direction; starting from 0 (b=0) to the total number in MDDW protocol, e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs)
# boolean is_called_from_global_recon is true, if this RunReconLoop is active
# If is_called_from_global_recon is false or not defined, the data selector needs to be defined here.
if !(@isdefined is_called_from_global_recon) || !is_called_from_global_recon
    global selector = Dict{Symbol,Any}()
    selector[:avg] = 1
    selector[:seg] = 1
    selector[:dif] = 0
end


## Data parameters (Path handling, data/results locations etc.)
# UHN work
# params_general[:project_path] = "C:\\Users\\Lars Kasper\\UHN\\Brain-TO - MRP-SPIDI - MRP-SPIDI\\SPIDI"
# Laptop home, one drive sync
# params_general[:project_path] = "C:\\Users\\kasperla\\UHN\\Brain-TO - MRP-SPIDI - MRP-SPIDI\\SPIDI"
# Gadgetron Server
# laptop home, external drive
# params_general[:data_path] = "e:\\SPIDI\\data\\SPIDI_0007\\Phantom\\rawdata"
#params_general[:project_path] = "/home/kasperl/SPIDI"
params_general[:project_path] = "/srv/data/ajaffray/TORONTO_COLLAB"

# SPIDI_0011

## Paths (user-dependent)
params_general[:data_path] = joinpath(params_general[:project_path], "data", "SPIDI_0007", "Human", "dat")
params_general[:gradients_path] = joinpath(params_general[:project_path], "data", "SPIDI_0007", "gradients")
params_general[:results_path] = joinpath(params_general[:project_path], "results", "SPIDI_0007", "Human")
params_general[:girf_path] = joinpath(params_general[:project_path], "code", "GIRFReco", "data", "GIRF", "GIRF_ISMRM2022")
params_general[:recon_save_path] = joinpath(params_general[:results_path], "recon", params_general[:recon_id])

## Files (user-dependent)

# Map scan file (Cartesian multi-echo file)
params_general[:map_scan_filename] = "field_map_132_2.h5"
params_general[:mapTEs_ms] = [4.92, 7.38]

params_general[:girf_filename] = ["2021Nov_PosNeg_Gx.mat", "2021Nov_PosNeg_Gy.mat", "2021Nov_PosNeg_Gz.mat"]

# File name for the spiral gradient
# multi-il (high-res 1.1mm) gradient file 508, single interleaf (low-res 2.6mm) gradient file 511
params_general[:gradient_filename] = joinpath("508", "gradients.txt")

# non-Cartesian (Spiral) scan file: MDDW 6
# params_general[:scan_filename]=["meas_MID00083_FID14976_diffSpiral_508_Intl0_b2k_4Avg.mrd"]
# params_general[:num_total_diffusion_directions] = 6

#  non-Cartesian (Spiral) scan file: MDDW30
# For single interleave data, use this section
# startIndexIntlv = 1 # Should always be 1 for single-interleave data.
# fname_spiralIntlv = "511_134_2.h5" # Gradient 511, b = 300, 10 diff directions
# fname_spiralIntlv = "511_136_2.h5" # Gradient 511, b = 700, 30 diff directions
# fname_spiralIntlv = "511_138_2.h5" # Gradient 511, b = 2500, 64 diff directions
# fname_spiralIntlv = "508_140_2.h5" # Gradient 508, interleave 0, b = 300, 10 diff directions
# fname_spiralIntlv = "508_142_2.h5" # Gradient 508, interleave 0, b = 700, 30 diff directions
# fname_spiralIntlv = "508_144_2.h5" # Gradient 508, interleave 0, b = 2500, 64 diff directions

# Multi-interleave data, needs all 4 file names, but will only read the corresponding one.
#fname_spiralIntlv0 = "508_124_2.h5" # Gradient 508, interleave 0, b = 2000, 6 diff directions, 4 averages
#fname_spiralIntlv1 = "508_126_2.h5" # Gradient 508, interleave 1, b = 2000, 6 diff directions, 4 averages
#fname_spiralIntlv2 = "508_128_2.h5" # Gradient 508, interleave 2, b = 2000, 6 diff directions, 4 averages
#fname_spiralIntlv3 = "508_130_2.h5" # Gradient 508, interleave 3, b = 2000, 6 diff directions, 4 averages
# params_general[:scan_filename]=["508_124_2.h5", "508_126_2.h5", "508_128_2.h5", "508_130_2.h5"]
# params_general[:num_total_diffusion_directions] = 6
# SPIDI_0007 MDDW 6
params_general[:scan_filename] = ["508_124_2.h5"]
params_general[:num_total_diffusion_directions] = 6

# NODDI
# params_general[:scan_filename]=["508_140_2.h5"]
# params_general[:num_total_diffusion_directions] = 10
# params_general[:scan_filename]=["508_142_2.h5"]
# params_general[:num_total_diffusion_directions] = 30
# params_general[:scan_filename]=["508_144_2.h5"]
# params_general[:num_total_diffusion_directions] = 64
# params_general[:scan_filename]=["508_124_2.h5"]




## Final, Automatic operations (dependent on previous sections, usually no need to change)

params_general[:gradient_fullpath] = joinpath(params_general[:gradients_path], params_general[:gradient_filename])

# NOTE: If loaded from other recon_id, this path might differ
params_general[:pathload_maps] = joinpath(params_general[:results_path], "recon", params_general[:recon_id])

# . makes join elementwise, i.e,. every file name (in array) with the same path
params_general[:girf_fullpath] = joinpath.(params_general[:girf_path], params_general[:girf_filename])

params_general[:map_scan_fullpath] = joinpath(params_general[:data_path], params_general[:map_scan_filename])
# . makes join elementwise, i.e,. every file name (in array) with the same path
params_general[:scan_fullpath] = joinpath.(params_general[:data_path], params_general[:scan_filename])

# filename for preprocessed data (remove oversampling, permute dimensions wrt MRIReco)
params_general[:processed_map_scan_filename] = "processedCartesianData.h5"
params_general[:processed_map_scan_fullpath] = joinpath(params_general[:recon_save_path], params_general[:processed_map_scan_filename])

params_general[:map_save_filename] = splitext(params_general[:map_scan_filename])[1] * "_reconmap.nii"
params_general[:sensitivity_save_filename] = splitext(params_general[:map_scan_filename])[1] * "_sensemap.nii"
params_general[:b0_map_save_filename] = splitext(params_general[:map_scan_filename])[1] * "_b0map.nii"

if isa(params_general[:scan_filename], AbstractVector)
    # for multiple files, concatenate recon name from scan file names, e.g., 508_124_2_508_126_2_508_128_2_508_130_2_recon.nii
    params_general[:recon_save_filename] =
        join([(x[1] * "_") for x in splitext.(params_general[:scan_filename])]) *
        "dif$(selector[:dif])_" *
        "itl$(selector[:seg])_" *
        "avg$(selector[:avg])_" *
        "recon.nii"
else
    # otherwise, just concat _recon.nii to file name
    params_general[:recon_save_filename] = splitext(params_general[:scan_filename])[1] * "_recon.nii"
end

params_general[:recon_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:recon_save_filename])
params_general[:map_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:map_save_filename])
params_general[:sensitivity_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:sensitivity_save_filename])
params_general[:b0_map_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:b0_map_save_filename])

if ~ispath(params_general[:recon_save_path])
    mkpath(params_general[:recon_save_path])
end

# copies this config file to the recon path for later checks of parameter functions
cp("recon/ReconConfig.jl", joinpath(params_general[:recon_save_path], "ReconConfig.jl"); force = true)

