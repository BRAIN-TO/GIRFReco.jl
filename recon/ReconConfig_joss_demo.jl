## This ReconConfig.jl file describes all reconstruction parameters, as well as data locations and selections for an iterative non-Cartesian reconstruction that relies 
#  on an external reference scan (Cartesian) to estimate calibration maps (coil sensitivities, B0 maps)

using Dates

paramsGeneral = Dict{Symbol,Any}()

## General options for recon script
paramsGeneral[:doLoadMaps] = false               # if true, reloads B0/SENSE maps instead of recalculating
paramsGeneral[:doSaveRecon] = true              # if true, saves reconstruction and all auxiliary image data (maps) as NIfTI files
paramsGeneral[:doPlotRecon] = true             # if true, plots intermediate debugging and output recon figures (needs graphics, not recommended in multi-thread mode due to PyPlot)
paramsGeneral[:doProcessMapScan] = true         # if true, compute sensitivity and B0 maps from reconstructed Cartesian scan   
paramsGeneral[:doSaveProcessedMapScan] = false; # save ISMRMD file of preprocessed Cartesian data (before recon)

## Reconstruction Parameters
# update time stamp for new recon, otherwise keep fixed, will create a new recon/<reconId> directory
#paramsGeneral[:reconId] = Dates.format(Dates.now(), "yyyy-mm-dd_HH_MM_SS") # recon ID is reconId
# paramsGeneral[:reconId] = "2022-10-20_09_07_07"

paramsGeneral[:reconId] = "vDemo";
paramsGeneral[:doCorrectWithB0map] = true
paramsGeneral[:doCorrectWithGIRFkxyz] = true
paramsGeneral[:doCorrectWithGIRFk0] = true

paramsGeneral[:doCoilCompression] = false;
paramsGeneral[:nVirtualCoils] = 8;
paramsGeneral[:fovShift] = [0, -20]; # Unit: number of voxels in [x,y] direction

## Scan parameters, Additional acquisition information, e.g., slice distance etc.
paramsGeneral[:sliceDistanceFactor_percent] = 400 # 000

#Total number of ADC points BEFORE the rewinder at the end of the spiral readout. For gradient 508, use 15655 (out of 16084); for gradient 511, use 15445 (out of 15624).
paramsGeneral[:numADCSamples] = 15655 # 15504
# Matrix size of the reconstructed image. For gradient 508 with all 4 interleaves, use 200 for high resolution image; otherwise consider using 112 or 84 for a lower resolution. The FOV is 220 mm for both gradients 508 and 511.
paramsGeneral[:reconSize] = [200, 200,1] #(112, 112) #(200, 200)
paramsGeneral[:nReconIterations] = 20; # number of recon iterations (for both Cartesian and Spiral recon)
paramsGeneral[:b0mapSmoothBeta] = 0.1 # for estimateB0Maps, * `β` - Regularization parameter controlling roughness penalty (larger = smoother, default 5e-4)
paramsGeneral[:doNormalizeRecon] = false # set max abs to 1
paramsGeneral[:scalingFactorSaveRecon] = 1.0e9 # 1 # typical range of recon intensities is 1e-7, rescale when saving, e.g., to 0...1000 roughly for fMRI analysis

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


# Set recon file paths (tell GIRFReco where to find things)

## Data parameters (Path handling, data/results locations etc.)
# UHN work
# paramsGeneral[:pathProject] = "C:\\Users\\Lars Kasper\\UHN\\Brain-TO - MRP-SPIDI - MRP-SPIDI\\SPIDI"
# Laptop home, one drive sync
# paramsGeneral[:pathProject] = "C:\\Users\\kasperla\\UHN\\Brain-TO - MRP-SPIDI - MRP-SPIDI\\SPIDI"
# Gadgetron Server
# laptop home, external drive
# paramsGeneral[:pathData] = "e:\\SPIDI\\data\\SPIDI_0007\\Phantom\\rawdata"
#paramsGeneral[:pathProject] = "/home/kasperl/SPIDI"
paramsGeneral[:pathProject] = "/srv/data/ajaffray/TORONTO_COLLAB"

## Paths (user-dependent)
paramsGeneral[:pathData] = joinpath(paramsGeneral[:pathProject], "data", "joss_demo", "Human", "dat")
paramsGeneral[:pathGradients] = joinpath(paramsGeneral[:pathProject], "data", "joss_demo", "gradients")
paramsGeneral[:pathResults] = joinpath(paramsGeneral[:pathProject], "results", "joss_demo", "Human")
paramsGeneral[:pathGIRF] = joinpath(paramsGeneral[:pathProject], "code", "GIRFReco", "data", "GIRF", "GIRF_ISMRM2022")
paramsGeneral[:pathSaveRecon] = joinpath(paramsGeneral[:pathResults], "recon", paramsGeneral[:reconId])

## Files (user-dependent)
# Map scan file (Cartesian multi-echo file)
paramsGeneral[:fileNameMapScan] = "field_map.h5"
paramsGeneral[:mapTEs_ms] = [4.92,  7.38]

paramsGeneral[:fileNameGIRF] = ["2021Nov_PosNeg_Gx.mat", "2021Nov_PosNeg_Gy.mat", "2021Nov_PosNeg_Gz.mat"]

# File name for the spiral gradient
# multi-il (high-res 1.1mm) gradient file 508, single interleaf (low-res 2.6mm) gradient file 511
paramsGeneral[:fileNameGradient] = joinpath("gradients.txt")

# non-Cartesian (Spiral) scan file: MDDW 6
# paramsGeneral[:fileNameScan]=["meas_MID00083_FID14976_diffSpiral_508_Intl0_b2k_4Avg.mrd"]
# paramsGeneral[:nDiffusionDirections] = 6

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
# paramsGeneral[:fileNameScan]=["508_124_2.h5", "508_126_2.h5", "508_128_2.h5", "508_130_2.h5"]
# paramsGeneral[:nDiffusionDirections] = 6
# SPIDI_0007 MDDW 6
paramsGeneral[:fileNameScan]=["spiral_data.h5"]
paramsGeneral[:nDiffusionDirections] = 6

# NODDI
# paramsGeneral[:fileNameScan]=["508_140_2.h5"]
# paramsGeneral[:nDiffusionDirections] = 10
# paramsGeneral[:fileNameScan]=["508_142_2.h5"]
# paramsGeneral[:nDiffusionDirections] = 30
# paramsGeneral[:fileNameScan]=["508_144_2.h5"]
# paramsGeneral[:nDiffusionDirections] = 64
# paramsGeneral[:fileNameScan]=["508_124_2.h5"]

## Final, Automatic operations (dependent on previous sections, usually no need to change)

paramsGeneral[:fullPathGradient] = joinpath(paramsGeneral[:pathGradients], paramsGeneral[:fileNameGradient])

# NOTE: If loaded from other reconId, this path might differ
paramsGeneral[:pathLoadMaps] = joinpath(paramsGeneral[:pathResults], "recon", paramsGeneral[:reconId])

# . makes join elementwise, i.e,. every file name (in array) with the same path
paramsGeneral[:fullPathGIRF] = joinpath.(paramsGeneral[:pathGIRF], paramsGeneral[:fileNameGIRF])

paramsGeneral[:fullPathMapScan] = joinpath(paramsGeneral[:pathData], paramsGeneral[:fileNameMapScan])
# . makes join elementwise, i.e,. every file name (in array) with the same path
paramsGeneral[:fullPathScan] = joinpath.(paramsGeneral[:pathData], paramsGeneral[:fileNameScan])

# filename for preprocessed data (remove oversampling, permute dimensions wrt MRIReco)
paramsGeneral[:fileNameProcessedMapScan] = "processed_cartesian_data.h5"
paramsGeneral[:fullPathProcessedMapScan] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameProcessedMapScan])

paramsGeneral[:fileNameSaveMapRecon] = splitext(paramsGeneral[:fileNameMapScan])[1] * "_reconmap.nii"
paramsGeneral[:fileNameSaveSense] = splitext(paramsGeneral[:fileNameMapScan])[1] * "_sensemap.nii"
paramsGeneral[:fileNameSaveB0] = splitext(paramsGeneral[:fileNameMapScan])[1] * "_b0map.nii"

if isa(paramsGeneral[:fileNameScan], AbstractVector)
    # for multiple files, concatenate recon name from scan file names, e.g., 508_124_2_508_126_2_508_128_2_508_130_2_recon.nii
    paramsGeneral[:fileNameSaveRecon] = join([(x[1] * "_") for x in splitext.(paramsGeneral[:fileNameScan])]) * "dif$(selector[:dif])_" * "itl$(selector[:seg])_" * "avg$(selector[:avg])_" * "recon.nii"
else
    # otherwise, just concat _recon.nii to file name
    paramsGeneral[:fileNameSaveRecon] = splitext(paramsGeneral[:fileNameScan])[1] * "_recon.nii"
end

paramsGeneral[:fullPathSaveRecon] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameSaveRecon] )
paramsGeneral[:fullPathSaveMapRecon] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameSaveMapRecon] )
paramsGeneral[:fullPathSaveSense] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameSaveSense] )
paramsGeneral[:fullPathSaveB0] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameSaveB0] )

if ~ispath(paramsGeneral[:pathSaveRecon])
    mkpath(paramsGeneral[:pathSaveRecon])
end

# copies this config file to the recon path for later checks of parameter functions
cp("recon/ReconConfig.jl", joinpath(paramsGeneral[:pathSaveRecon], "ReconConfig.jl"); force = true)

