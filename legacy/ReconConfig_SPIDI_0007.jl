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
params_general[:doCorrectWithGIRFkxyz] = true
params_general[:doCorrectWithGIRFk0] = true

params_general[:nVirtualCoils] = 8;
params_general[:doCoilCompression] = false;
params_general[:fovShift] = [0, -20]; # Unit: number of voxels

## Scan parameters, Additional acquisition information, e.g., slice distance etc.
params_general[:sliceDistanceFactor_percent] = 400 # 000

#Total number of ADC points BEFORE the rewinder at the end of the spiral readout. For gradient 508, use 15655 (out of 16084); for gradient 511, use 15445 (out of 15624).
params_general[:numADCSamples] = 15655 # 15504
# Matrix size of the reconstructed image. For gradient 508 with all 4 interleaves, use 200 for high resolution image; otherwise consider using 112 or 84 for a lower resolution. The FOV is 220 mm for both gradients 508 and 511.
params_general[:reconSize] = [200, 200, 1] #(112, 112) #(200, 200)
params_general[:nReconIterations] = 20; # number of recon iterations (for both Cartesian and Spiral recon)
params_general[:b0mapSmoothBeta] = 0.1 # for estimateB0Maps, * `β` - Regularization parameter controlling roughness penalty (larger = smoother, default 5e-4)
params_general[:doNormalizeRecon] = false # set max abs to 1
params_general[:scalingFactorSaveRecon] = 1.0e9 # 1 # typical range of recon intensities is 1e-7, rescale when saving, e.g., to 0...1000 roughly for fMRI analysis

# Data selector
#  Choose diffusion direction; starting from 0 (b=0) to the total number in MDDW protocol, e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs)
# boolean isCalledFromReconLoopGlobal is true, if this RunReconLoop is active
# If isCalledFromReconLoopGlobal is false or not defined, the data selector needs to be defined here.
if !(@isdefined isCalledFromReconLoopGlobal) || !isCalledFromReconLoopGlobal
    global selector = Dict{Symbol,Any}()
    selector[:avg] = 1
    selector[:seg] = 1
    selector[:dif] = 0
end


## Data parameters (Path handling, data/results locations etc.)
# UHN work
# params_general[:pathProject] = "C:\\Users\\Lars Kasper\\UHN\\Brain-TO - MRP-SPIDI - MRP-SPIDI\\SPIDI"
# Laptop home, one drive sync
# params_general[:pathProject] = "C:\\Users\\kasperla\\UHN\\Brain-TO - MRP-SPIDI - MRP-SPIDI\\SPIDI"
# Gadgetron Server
# laptop home, external drive
# params_general[:pathData] = "e:\\SPIDI\\data\\SPIDI_0007\\Phantom\\rawdata"
#params_general[:pathProject] = "/home/kasperl/SPIDI"
params_general[:pathProject] = "/srv/data/ajaffray/TORONTO_COLLAB"

# SPIDI_0011

## Paths (user-dependent)
params_general[:pathData] = joinpath(params_general[:pathProject], "data", "SPIDI_0007", "Human", "dat")
params_general[:pathGradients] = joinpath(params_general[:pathProject], "data", "SPIDI_0007", "gradients")
params_general[:pathResults] = joinpath(params_general[:pathProject], "results", "SPIDI_0007", "Human")
params_general[:pathGIRF] = joinpath(params_general[:pathProject], "code", "GIRFReco", "data", "GIRF", "GIRF_ISMRM2022")
params_general[:pathSaveRecon] = joinpath(params_general[:pathResults], "recon", params_general[:recon_id])

## Files (user-dependent)

# Map scan file (Cartesian multi-echo file)
params_general[:fileNameMapScan] = "field_map_132_2.h5"
params_general[:mapTEs_ms] = [4.92, 7.38]

params_general[:fileNameGIRF] = ["2021Nov_PosNeg_Gx.mat", "2021Nov_PosNeg_Gy.mat", "2021Nov_PosNeg_Gz.mat"]

# File name for the spiral gradient
# multi-il (high-res 1.1mm) gradient file 508, single interleaf (low-res 2.6mm) gradient file 511
params_general[:fileNameGradient] = joinpath("508", "gradients.txt")

# non-Cartesian (Spiral) scan file: MDDW 6
# params_general[:fileNameScan]=["meas_MID00083_FID14976_diffSpiral_508_Intl0_b2k_4Avg.mrd"]
# params_general[:nDiffusionDirections] = 6

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
# params_general[:fileNameScan]=["508_124_2.h5", "508_126_2.h5", "508_128_2.h5", "508_130_2.h5"]
# params_general[:nDiffusionDirections] = 6
# SPIDI_0007 MDDW 6
params_general[:fileNameScan] = ["508_124_2.h5"]
params_general[:nDiffusionDirections] = 6

# NODDI
# params_general[:fileNameScan]=["508_140_2.h5"]
# params_general[:nDiffusionDirections] = 10
# params_general[:fileNameScan]=["508_142_2.h5"]
# params_general[:nDiffusionDirections] = 30
# params_general[:fileNameScan]=["508_144_2.h5"]
# params_general[:nDiffusionDirections] = 64
# params_general[:fileNameScan]=["508_124_2.h5"]




## Final, Automatic operations (dependent on previous sections, usually no need to change)

params_general[:fullPathGradient] = joinpath(params_general[:pathGradients], params_general[:fileNameGradient])

# NOTE: If loaded from other recon_id, this path might differ
params_general[:pathload_maps] = joinpath(params_general[:pathResults], "recon", params_general[:recon_id])

# . makes join elementwise, i.e,. every file name (in array) with the same path
params_general[:fullPathGIRF] = joinpath.(params_general[:pathGIRF], params_general[:fileNameGIRF])

params_general[:fullPathMapScan] = joinpath(params_general[:pathData], params_general[:fileNameMapScan])
# . makes join elementwise, i.e,. every file name (in array) with the same path
params_general[:fullPathScan] = joinpath.(params_general[:pathData], params_general[:fileNameScan])

# filename for preprocessed data (remove oversampling, permute dimensions wrt MRIReco)
params_general[:fileNameProcessedMapScan] = "processedCartesianData.h5"
params_general[:fullPathProcessedMapScan] = joinpath(params_general[:pathSaveRecon], params_general[:fileNameProcessedMapScan])

params_general[:fileNamesave_mapRecon] = splitext(params_general[:fileNameMapScan])[1] * "_reconmap.nii"
params_general[:fileNameSaveSense] = splitext(params_general[:fileNameMapScan])[1] * "_sensemap.nii"
params_general[:fileNameSaveB0] = splitext(params_general[:fileNameMapScan])[1] * "_b0map.nii"

if isa(params_general[:fileNameScan], AbstractVector)
    # for multiple files, concatenate recon name from scan file names, e.g., 508_124_2_508_126_2_508_128_2_508_130_2_recon.nii
    params_general[:fileNameSaveRecon] =
        join([(x[1] * "_") for x in splitext.(params_general[:fileNameScan])]) *
        "dif$(selector[:dif])_" *
        "itl$(selector[:seg])_" *
        "avg$(selector[:avg])_" *
        "recon.nii"
else
    # otherwise, just concat _recon.nii to file name
    params_general[:fileNameSaveRecon] = splitext(params_general[:fileNameScan])[1] * "_recon.nii"
end

params_general[:fullPathSaveRecon] = joinpath(params_general[:pathSaveRecon], params_general[:fileNameSaveRecon])
params_general[:fullPathsave_mapRecon] = joinpath(params_general[:pathSaveRecon], params_general[:fileNamesave_mapRecon])
params_general[:fullPathSaveSense] = joinpath(params_general[:pathSaveRecon], params_general[:fileNameSaveSense])
params_general[:fullPathSaveB0] = joinpath(params_general[:pathSaveRecon], params_general[:fileNameSaveB0])

if ~ispath(params_general[:pathSaveRecon])
    mkpath(params_general[:pathSaveRecon])
end

# copies this config file to the recon path for later checks of parameter functions
cp("recon/ReconConfig.jl", joinpath(params_general[:pathSaveRecon], "ReconConfig.jl"); force = true)

