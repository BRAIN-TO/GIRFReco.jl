#-----------------------------------------------------------------------------------
# # [GIRFReco.jl Example Configuration File](@id example_config)
#-----------------------------------------------------------------------------------

#= This ReconConfig.jl file describes all general reconstruction parameters, 
as well as data locations and selections for an iterative non-Cartesian reconstruction that relies 
on an external reference scan (Cartesian) to estimate calibration maps (coil sensitivities, B0 maps).

All parameters are stored in a dictionary named `paramsGeneral`.

=#
paramsGeneral = Dict{Symbol,Any}()

# Gyromagnetic ratio, in unit of Hz
paramsGeneral[:gamma] = 42577478

# General options for reconstruction script
paramsGeneral[:doLoadMaps] = false              # if true, reloads B0/SENSE maps instead of recalculating
paramsGeneral[:doSaveRecon] = true              # if true, saves reconstruction and all auxiliary image data (maps) as NIfTI files
paramsGeneral[:doPlotRecon] = true              # if true, plots intermediate debugging and output recon figures (needs graphics, not recommended in multi-thread mode due to PyPlot)
paramsGeneral[:doProcessMapScan] = true         # if true, compute sensitivity and B0 maps from reconstructed Cartesian scan   
paramsGeneral[:doSaveProcessedMapScan] = false  # save ISMRMD file of preprocessed Cartesian data (before recon)
paramsGeneral[:reconId] = "vDemo"               # unique identifier for the saved result files
paramsGeneral[:doCorrectWithB0map] = true       # whether perform off-resonance correction
paramsGeneral[:doCorrectWithGIRFkxyz] = true    # whether perform 1st order GIRF correction
paramsGeneral[:doCorrectWithGIRFk0] = true      # whether perform 1st order GIRF correction
paramsGeneral[:doCoilCompression] = false       # whether perform coil compression
paramsGeneral[:nVirtualCoils] = 8;              # if perform coil compression, the number of coils to be compressed to
paramsGeneral[:fovShift] = [0, -20];            # amount of FOV shift; in unit of number of voxels in [x,y] direction
paramsGeneral[:sliceDistanceFactor_percent] = 400 # Scan parameters, Additional acquisition information, e.g., slice distance etc.

#Total number of ADC points BEFORE the rewinder at the end of the spiral readout. For gradient 508, use 15655 (out of 16084); for gradient 511, use 15445 (out of 15624).
paramsGeneral[:numADCSamples] = 15655 # 15504
# Matrix size of the reconstructed image. For gradient 508 with all 4 interleaves, use 200 for high resolution image; otherwise consider using 112 or 84 for a lower resolution. The FOV is 220 mm for both gradients 508 and 511.
paramsGeneral[:reconSize] = [200, 200,1] #(112, 112) #(200, 200)
paramsGeneral[:nReconIterations] = 20;          # number of recon iterations (for both Cartesian and Spiral recon)
paramsGeneral[:b0mapSmoothBeta] = 0.1           # for estimateB0Maps, * `β` - Regularization parameter controlling roughness penalty (larger = smoother, default 5e-4)
paramsGeneral[:doNormalizeRecon] = false        # set max abs to 1
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

# Set paths for input and output files
paramsGeneral[:pathProject] = rootProjPath # Root path for the project

#src # Paths (user-dependent)
#src # paramsGeneral[:pathData] = joinpath(paramsGeneral[:pathProject], "data", "joss_demo", "Human", "dat")
#src # paramsGeneral[:pathGradients] = joinpath(paramsGeneral[:pathProject], "data", "joss_demo", "gradients")
#src # paramsGeneral[:pathResults] = joinpath(paramsGeneral[:pathProject], "results", "joss_demo", "Human")
#src # paramsGeneral[:pathGIRF] = joinpath(paramsGeneral[:pathProject], "code", "GIRFReco", "data", "GIRF", "GIRF_ISMRM2022")
#src # paramsGeneral[:pathSaveRecon] = joinpath(paramsGeneral[:pathResults], "recon", paramsGeneral[:reconId])

#Path to ISMRMRD files (raw k-space data)
paramsGeneral[:pathData] = joinpath(paramsGeneral[:pathProject], "data", "SPIDI_0007", "Human", "dat")
#Path to spiral readout gradient files
paramsGeneral[:pathGradients] = joinpath(paramsGeneral[:pathProject], "data", "SPIDI_0007", "gradients")
#Path to GIRF files
paramsGeneral[:pathGIRF] = joinpath(paramsGeneral[:pathProject], "code", "GIRFReco", "data", "GIRF", "GIRF_ISMRM2022")

#Path to middle results (coil and B0 maps) files
paramsGeneral[:pathResults] = joinpath(paramsGeneral[:pathProject], "results", "SPIDI_0007", "Human")
#Path to final reconstructed spiral images
paramsGeneral[:pathSaveRecon] = joinpath(paramsGeneral[:pathResults], "recon", paramsGeneral[:reconId])


# Map scan file (Cartesian multi-echo file)
paramsGeneral[:fileNameMapScan] = "field_map_132_2.h5"
paramsGeneral[:mapTEs_ms] = [4.92,  7.38]

paramsGeneral[:fileNameGIRF] = ["2021Nov_PosNeg_Gx.mat", "2021Nov_PosNeg_Gy.mat", "2021Nov_PosNeg_Gz.mat"]

# File name for the spiral gradient
# multi-il (high-res 1.1mm) gradient file 508, single interleaf (low-res 2.6mm) gradient file 511
paramsGeneral[:fileNameGradient] = joinpath("508", "gradients.txt")

# SPIDI_0007 MDDW 6
paramsGeneral[:fileNameScan]=["508_124_2.h5"]
paramsGeneral[:nDiffusionDirections] = 6

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

