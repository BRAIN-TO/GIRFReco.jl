#-----------------------------------------------------------------------------------
# # [GIRFReco.jl Example Configuration File](@id example_config)
#-----------------------------------------------------------------------------------

#=
## Introduction

This ReconConfig.jl file describes all general reconstruction parameters, 
as well as data locations and selections for an iterative non-Cartesian reconstruction that relies 
on an external reference scan (Cartesian) to estimate calibration maps (coil sensitivities, B0 maps).

All parameters are stored in a dictionary named `paramsGeneral`.

=#
paramsGeneral = Dict{Symbol,Any}();

# Gyromagnetic ratio, in unit of Hz
paramsGeneral[:gamma] = 42577478;

# General options for reconstruction script:
paramsGeneral[:doLoadMaps] = false              # if true, reloads B0/SENSE maps instead of recalculating
paramsGeneral[:doSaveRecon] = true              # if true, saves reconstruction and all auxiliary image data (maps) as NIfTI files
paramsGeneral[:doPlotRecon] = true              # if true, plots intermediate debugging and output recon figures (needs graphics, not recommended in multi-thread mode due to PyPlot)
paramsGeneral[:doProcessMapScan] = true         # if true, compute sensitivity and B0 maps from reconstructed Cartesian scan   
paramsGeneral[:doSaveProcessedMapScan] = false  # save ISMRMD file of preprocessed Cartesian data (before recon)
paramsGeneral[:reconId] = "vFinal_unshifted_2_k0"               # unique identifier for the saved result files
paramsGeneral[:doCorrectWithB0map] = true       # whether perform off-resonance correction
paramsGeneral[:doCorrectWithGIRFkxyz] = true    # whether perform 1st order GIRF correction
paramsGeneral[:doCorrectWithGIRFk0] = true    # whether perform 1st order GIRF correction
paramsGeneral[:doCoilCompression] = false       # whether perform coil compression
paramsGeneral[:doNormalizeRecon] = false       # if true, set the range of magnitude image as [0 1]

# General parameters for reconstruction:
paramsGeneral[:numADCSamples] = 15650             # total number of ADC points BEFORE the rewinder at the end of the spiral readout. Need to check with data. # 15504 for 523, 15655 for gradient 508
paramsGeneral[:reconSize] = [200,200, 1]          # the matrix size of the reconstructed images. Needs to specify 1 on Z dimension for 2D images
paramsGeneral[:nReconIterations] = 40             # number of recon iterations (for both Cartesian and Spiral recon)
paramsGeneral[:b0mapSmoothBeta] = 0.1             # for estimateB0Maps, * `β` - Regularization parameter controlling roughness penalty (larger = smoother, default 5e-4)
paramsGeneral[:scalingFactorSaveRecon] = 1.0e9    # typical range of recon intensities is 1e-7, rescale when saving, e.g., to 0...1000 roughly for fMRI analysis
paramsGeneral[:nVirtualCoils] = 8                 # if perform coil compression, the number of coils to be compressed to
paramsGeneral[:fovShift] = [-2, -5]              # amount of FOV shift; in unit of number of voxels in [x,y] direction
paramsGeneral[:sliceDistanceFactor_percent] = 400 # Scan parameters, Additional acquisition information, e.g., slice distance etc.
paramsGeneral[:nDiffusionDirections] = 6          # Need to specify total diffusion directions included in the raw data
#=

## Data selector

See details in the page Advanced Usage.

Data selector is designed for calling the script repetitively through [`RunReconLoop.jl`](@__REPO_ROOT_URL__/recon/RunReconLoop.jl) 
for e.g. different diffusion directions and/or different averages. 
Note that the index of diffusion direction starts from 0 (b=0) to the total number in MDDW protocol, 
e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs. 
`boolean isCalledFromReconLoopGlobal` should be set as true if this `RunReconLoop.jl` is active.
If isCalledFromReconLoopGlobal is false or not defined, the data selector needs to be defined in the following code block.
=#
if !(@isdefined isCalledFromReconLoopGlobal) || !isCalledFromReconLoopGlobal
    global selector = Dict{Symbol,Any}()
    selector[:avg] = 1
    selector[:seg] = 1
    selector[:dif] = 0;
end


#=
## Directories and File Names

We set directories for data reading in and results writing in this section. 

### Specifying Directories
=#

paramsGeneral[:pathProject] = rootProjPath # Root path for the project

#src # Paths (user-dependent)
#src # paramsGeneral[:pathData] = joinpath(paramsGeneral[:pathProject], "data", "joss_demo", "Human", "dat")
#src # paramsGeneral[:pathGradients] = joinpath(paramsGeneral[:pathProject], "data", "joss_demo", "gradients")
#src # paramsGeneral[:pathResults] = joinpath(paramsGeneral[:pathProject], "results", "joss_demo", "Human")
#src # paramsGeneral[:pathGIRF] = joinpath(paramsGeneral[:pathProject], "code", "GIRFReco", "data", "GIRF", "GIRF_ISMRM2022")
#src # paramsGeneral[:pathSaveRecon] = joinpath(paramsGeneral[:pathResults], "recon", paramsGeneral[:reconId])

#Path to ISMRMRD files (raw k-space data) [Input]
paramsGeneral[:pathData] = paramsGeneral[:pathProject]
#Path to spiral readout gradient files [Input]
paramsGeneral[:pathGradients] = joinpath(paramsGeneral[:pathProject],"Gradients")
#Path to GIRF files [Input]
paramsGeneral[:pathGIRF] = joinpath(paramsGeneral[:pathProject], "GIRF", "GIRF_ISMRM2022")
#Path to middle results (coil and B0 maps) files [Output]
paramsGeneral[:pathResults] = joinpath(paramsGeneral[:pathProject], "results","phantom")
#Path to final reconstructed spiral images [Output]
paramsGeneral[:pathSaveRecon] = joinpath(paramsGeneral[:pathResults], "recon", paramsGeneral[:reconId]);

#=
### Specifying File Names
=#
paramsGeneral[:fileNameMapScan] = "Fieldmaps/meas_MID00083_FID06181_GRE_FieldMap_DualEcho_2mm.mrd" # Cartesian dual-echo file, for coil and B0 maps calculation [Input]
paramsGeneral[:fileNameMapStem] = "meas_MID00083_FID06181_GRE_FieldMap_DualEcho_2mm.mrd"
paramsGeneral[:mapTEs_ms] = [4.92,  7.38] # Two echo times, in ms
paramsGeneral[:fileNameGIRF] = ["2021Nov_PosNeg_Gx.mat", "2021Nov_PosNeg_Gy.mat", "2021Nov_PosNeg_Gz.mat"] # Calculated GIRF for each gradient axis [Input]
paramsGeneral[:fileNameGradient] = joinpath("gradients508.txt") # File name for the spiral gradient [Input]
paramsGeneral[:fileNameScan] = ["Spirals/meas_MID00072_FID06170_diffSpiral_508_Intl0_b2k_4Avg.mrd","Spirals/meas_MID00074_FID06172_diffSpiral_508_Intl1_b2k_4Avg.mrd","Spirals/meas_MID00076_FID06174_diffSpiral_508_Intl2_b2k_4Avg.mrd","Spirals/meas_MID00078_FID06176_diffSpiral_508_Intl3_b2k_4Avg.mrd"]# , "Spirals/508_74_2.h5","Spirals/508_76_2.h5","Spirals/508_78_2.h5"] # ISMRMRD Raw k-space data for spiral acquisition [Input]
paramsGeneral[:fileNameScanStem] = "meas_MID00072_FID06170_diffSpiral_508_Intl0_b2k_4Avg.mrd"
paramsGeneral[:fileNameProcessedMapScan] = "processed_cartesian_data.h5" # file name for preprocessed data (remove oversampling, permute dimensions wrt MRIReco) [Output]
paramsGeneral[:fileNameSaveMapRecon] = splitext(paramsGeneral[:fileNameMapStem])[1] * "_reconmap.nii" # File name for reconstructed dual-echo Cartesian images [Output]
paramsGeneral[:fileNameSaveSense] = splitext(paramsGeneral[:fileNameMapStem])[1] * "_sensemap.nii" # File name for calculated coil sensitivity maps [Output]
paramsGeneral[:fileNameSaveB0] = splitext(paramsGeneral[:fileNameMapStem])[1] * "_b0map.nii"; # File name for calculated off-resonance (B0) maps [Output]

#=
File name for the final reconstructed spiral image.
If we reconstructing multiple spiral data files (e.g. multiple interleaves) through `RunReconLoop.jl`, 
the file name for the final reconstructed image is concatenated from multiple scan file names. 
Otherwise, just append `_recon.nii` as suffix to file name.
=#
if isa(paramsGeneral[:fileNameScan], AbstractVector)
    paramsGeneral[:fileNameSaveRecon] = join([(x[1] * "_") for x in splitext.(paramsGeneral[:fileNameScanStem])]) * "dif$(selector[:dif])_" * "itl$(selector[:seg])_" * "avg$(selector[:avg])_" * "recon.nii";
else
    paramsGeneral[:fileNameSaveRecon] = splitext(paramsGeneral[:fileNameScanStem])[1] * "_recon.nii";
end

#=
### Assembling Full Paths

Assembling directories and file names for final full pathes. 
These are automated operations.
=#
paramsGeneral[:fullPathGradient] = joinpath(paramsGeneral[:pathGradients], paramsGeneral[:fileNameGradient]) # Full paths of spiral readout gradients
#src # paramsGeneral[:pathLoadMaps] = 
#src # joinpath(paramsGeneral[:pathResults], "recon", paramsGeneral[:reconId]) # NOTE: If loaded from other reconId, this path might differ
paramsGeneral[:fullPathGIRF] = joinpath.(paramsGeneral[:pathGIRF], paramsGeneral[:fileNameGIRF]) # Full paths of GIRF files
paramsGeneral[:fullPathMapScan] = joinpath(paramsGeneral[:pathData], paramsGeneral[:fileNameMapScan]) # Full path of dual-echo Cartesian data
paramsGeneral[:fullPathScan] = joinpath.(paramsGeneral[:pathData], paramsGeneral[:fileNameScan]) # Full paths of raw k-space data files of spiral acquisition
paramsGeneral[:fullPathProcessedMapScan] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameProcessedMapScan]) # Full paths of pre-processed Cartesian dual-echo data [Output]
paramsGeneral[:fullPathSaveRecon] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameSaveRecon] ) # Full paths of the reconstructed spiral image [Output]
paramsGeneral[:fullPathSaveMapRecon] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameSaveMapRecon] ) # Full paths of reconstructed dual-echo Cartesian images [Output]
paramsGeneral[:fullPathSaveSense] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameSaveSense] ) # Full paths of calculated coil sensitivity maps [Output]
paramsGeneral[:fullPathSaveB0] = joinpath(paramsGeneral[:pathSaveRecon], paramsGeneral[:fileNameSaveB0] ); # Full paths of calculated off-resonance (B0) maps [Output]

#=
## Final Steps

If the path for results writing is not existing, create it.

As the last step of configuration, copy this config file 
to the recon path for further checking and debugging purposes.
=#
if ~ispath(paramsGeneral[:pathSaveRecon])
    mkpath(paramsGeneral[:pathSaveRecon])
end

cp("recon/ReconConfig.jl", joinpath(paramsGeneral[:pathSaveRecon], "ReconConfig.jl"); force = true)

