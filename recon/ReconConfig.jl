using Dates

paramsGeneral = Dict{Symbol,Any}()

# UHN work
# paramsGeneral[:pathProject] = "C:\\Users\\Lars Kasper\\UHN\\Brain-TO - MRP-SPIDI - MRP-SPIDI\\SPIDI"
# Laptop home, one drive sync
# paramsGeneral[:pathProject] = "C:\\Users\\kasperla\\UHN\\Brain-TO - MRP-SPIDI - MRP-SPIDI\\SPIDI"
# Gadgetron Server
# laptop home, external drive
# paramsGeneral[:pathData] = "e:\\SPIDI\\data\\SPIDI_0007\\Phantom\\rawdata"

## General options for recon script
paramsGeneral[:doLoadMaps] = true
paramsGeneral[:doSaveRecon] = true
paramsGeneral[:doPlotRecon] = false

## Reconstruction Parameter
# update time stamp for new recon, otherwise keep fixed, will create a new recon/<reconId> directory
#paramsGeneral[:reconId] = Dates.format(Dates.now(), "yyyy-mm-dd_HH_MM_SS") # recon ID is reconId
# paramsGeneral[:reconId] = "2022-10-20_09_07_07"
paramsGeneral[:reconId] = "v11";
paramsGeneral[:doCorrectWithB0map] = true
paramsGeneral[:doCorrectWithGIRFkxyz] = true
paramsGeneral[:doCorrectWithGIRFk0] = true

paramsGeneral[:nVirtualCoils] = 8;
paramsGeneral[:doCoilCompression] = false;
paramsGeneral[:fovShift] = [0, -10]; # TODO: identify unit
## Scan parameters, Additional acquisition information, e.g., slice distance etc.
paramsGeneral[:sliceDistanceFactor_percent] = 000
## Total number of ADC points BEFORE the rewinder at the end of the spiral readout. For gradient 508, use 15655 (out of 16084); for gradient 511, use 15445 (out of 15624).
# numADCSamples = 15504
# numADCSamples = 15655
paramsGeneral[:numADCSamples] = 15655
# Matrix size of the reconstructed image. For gradient 508 with all 4 interleaves, use 200 for high resolution image; otherwise consider using 112 or 84 for a lower resolution. The FOV is 220 mm for both gradients 508 and 511.
paramsGeneral[:reconSize] = (200, 200) #(112, 112) #(200, 200)
paramsGeneral[:b0mapSmoothBeta] = 0.01 # for estimateB0Maps, * `β` - Regularization parameter controlling roughness penalty (larger = smoother, default 5e-4)

## Data selector
#  Choose diffusion direction; starting from 0 (b=0) to the total number in MDDW protocol, e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs)
# boolean isCalledFromReconLoopGlobal is true, if this RunReconLoop is active
if !((@isdefined isCalledFromReconLoopGlobal) && isCalledFromReconLoopGlobal)
    global selector = Dict{Symbol,Any}()
    selector[:avg] = 1;
    selector[:seg] = 1;
    selector[:dif] = 0;
end


## Path handling, data/results locations etc.
paramsGeneral[:pathProject] = "/home/kasperl/SPIDI"

# SPIDI_0011
paramsGeneral[:pathData] = joinpath("/home/wuz/spiralDiffusion/", "SPIDI_0011")
paramsGeneral[:pathGradients] = joinpath(paramsGeneral[:pathProject], "data", "SPIDI_0007", "gradients")
paramsGeneral[:pathResults] = joinpath(paramsGeneral[:pathProject], "results", "SPIDI_0011")
paramsGeneral[:pathGIRF] = joinpath(paramsGeneral[:pathProject], "code", "GIRFReco", "data", "GIRF", "GIRF_ISMRM2022")
paramsGeneral[:pathSaveRecon] = joinpath(paramsGeneral[:pathResults], "recon", paramsGeneral[:reconId])
paramsGeneral[:fileNameMapScan] = "meas_MID00075_FID14968_GRE_FieldMap_DualEcho_2mm.mrd"
paramsGeneral[:mapTEs_ms] = [4.92,  7.38]

# SPIDI_0011 
# MDDW 6
# paramsGeneral[:fileNameScan]=["meas_MID00083_FID14976_diffSpiral_508_Intl0_b2k_4Avg.mrd"]
# paramsGeneral[:nDiffusionDirections] = 6

# MDDW30
paramsGeneral[:fileNameScan]=["meas_MID00085_FID14978_diffSpiral_508_Intl0_b1000_30d.mrd"]
paramsGeneral[:nDiffusionDirections] = 30

# File name for the spiral gradient
# multi-il gradient file 508, single interleaf gradient file 511
paramsGeneral[:fullPathGradient] = joinpath(paramsGeneral[:pathGradients], "508", "gradients.txt")

if ~ispath(paramsGeneral[:pathSaveRecon])
    mkpath(paramsGeneral[:pathSaveRecon])
end

# copy this file to the recon path for later checks of parameter functions
cp("recon/ReconConfig.jl", joinpath(paramsGeneral[:pathSaveRecon], "ReconConfig.jl"); force = true)

# If loaded from other scan, this path might differ
paramsGeneral[:pathLoadMaps] = joinpath(paramsGeneral[:pathResults], "recon", paramsGeneral[:reconId])

paramsGeneral[:fileNameGIRF] = ["2021Nov_PosNeg_Gx.mat", "2021Nov_PosNeg_Gy.mat", "2021Nov_PosNeg_Gz.mat"]
# . makes join elementwise, i.e,. every file name (in array) with the same path
paramsGeneral[:fullPathGIRF] = joinpath.(paramsGeneral[:pathGIRF], paramsGeneral[:fileNameGIRF])

paramsGeneral[:fullPathMapScan] = joinpath(paramsGeneral[:pathData], paramsGeneral[:fileNameMapScan])
# . makes join elementwise, i.e,. every file name (in array) with the same path
paramsGeneral[:fullPathScan] = joinpath.(paramsGeneral[:pathData], paramsGeneral[:fileNameScan])

paramsGeneral[:doProcessMapScan] = true
paramsGeneral[:doSaveProcessedMapScan] = false;
# filename for preprocessed data (remove oversampling, permute dimensions wrt MRIReco)
paramsGeneral[:fileNameProcessedMapScan] = "processedCartesianData.h5"
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