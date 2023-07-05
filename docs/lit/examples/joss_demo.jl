#-----------------------------------------------------------------------------------
# # [GIRFReco.jl Example Script](@id example_script)
#-----------------------------------------------------------------------------------

#=
This page demonstrates an example script for using GIRFReco.jl

This page was generated from the following Julia file: [`joss_demo.jl`](@__REPO_ROOT_URL__/docs/lit/examples/joss_demo.jl)

The configuration file is [`ReconConfig_joss_demo.jl`](@__REPO_ROOT_URL__/docs/lit/examples/ReconConfig_joss_demo.jl)
=#

#=
## 1. Setup

The necessary Julia packages needed for spiral reconstruction.
=#

#Base packages for computation
using HDF5, LinearAlgebra, Dierckx, DSP, FourierTools, RegularizedLeastSquares, ImageUtils, PolygonInbounds

#Packages for figure displaying
using MosaicViews, Plots, Images

#Our developed packages
using GIRFReco, MRIGradients

#MRIReco and its sub-packages
using MRIReco, FileIO, MRIFiles, MRICoilSensitivities

#=
## 2. Configurations for reconstruction

The following file, [`ReconConfig_joss_demo.jl`](@__REPO_ROOT_URL__/docs/lit/examples/ReconConfig_joss_demo.jl),
includes general configuration for spiral reconstruction.
It is necessary to execute this file to make sure all parameters are loaded.
Sample Data that works with this script can be found at: https://doi.org/10.5281/zenodo.7779044
Please download, extract and set the rootProjPath as the top level folder (should be something like /your/path/here/data-2, I've renamed mine to SPIDI)

=#

rootProjPath = "Your/Extracted/Data/Folder" # Root path of the data extracted from Zenodo

include("ReconConfig_joss_demo.jl")


# Two user defined parameters, just for this script.
reloadSpiralData = true; # Set true if we need to reload raw data compulsively.
reloadGIRFData = true; # Set true if we need to reload GIRF data compulsively.

#=
Choose Slice ([single number] OR [1,2,31,...]）
Leave empty ([]) or remove this line to later select all slices
=#
sliceChoice = [];

#=
Choose which diffusion directions and averages to be processed. 
Diffusion direction index starts from 0 (b=0) to the total number in MDDW protocol (e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs). 
Index for average starts from 1.
=#
diffusionDirection = 0
idxAverage = 1
nTotalDiffDir = paramsGeneral[:nDiffusionDirections]

## Determine to reconstruct single-interleave data, or one interleave out of multi-interleave data.
isDataSingleIntlv = isa(paramsGeneral[:fullPathScan], String)

#=
Choose which interleave to be reconstructed. 
For multi-interleave data, the range of this value is [1:nTotalNumIntlv] (total number of interleaves)
For single-interleave data, it should always be set as 1; for multi-interleave data, the value set here will be used, indicating which interleaves to be merged and reconstructed.
=#
startIndexIntlv = 1

#===================================================
## 3. Image Reconstruction

The steps of image reconstruction starts here.

### 3.1 Calculation of B0 and Coil Sensitivity Maps

The first step in reconstruction pipeline is to calculate the off-resonance (B0) maps `b0Maps` 
and coil sensitivity maps `senseCartesian` through the Cartesian reconstruction script 
[CartesianRecon.jl](@__REPO_ROOT_URL__/recon/CartesianRecon.jl). 
Ideally this script is execute once and the calculated maps are 
saved into files, which are loaded for future usage to save calculation time. 
This is controlled by `doLoadMaps` in general parameters. 
===================================================#

if paramsGeneral[:doLoadMaps] && isfile(paramsGeneral[:fullPathSaveB0])
    @info "Loading SENSE and B0 maps from $(paramsGeneral[:fullPathSaveSense]) and $(paramsGeneral[:fullPathSaveB0])"
    b0Maps = loadMap(paramsGeneral[:fullPathSaveB0])
    
    nSlices = size(b0Maps, 3);
    sliceIndexArray = getSliceOrder(nSlices, isSliceInterleaved = true)

    b0Maps = b0Maps[:,:,invperm(sliceIndexArray)]
    senseCartesian = loadMap(paramsGeneral[:fullPathSaveSense]; doSplitPhase = true)[:,:,invperm(sliceIndexArray),:]
else
    @info "Running CartesianRecon to retrieve maps (senseCartesian and b0Maps)"
    include("../../../recon/CartesianRecon.jl")
    nSlices = size(b0Maps, 3);
end

#=
### 3.2 Preparation of Spiral Reconstruction

With off-resonance (B0) maps and coil sensitivity maps calculated, 
before the reconstruction of spiral images, there are necessary steps to prepare for 
the related data. 

#### 3.2.1 Data Selection

The first step is to select the part of spiral k-space data that we 
would like to reconstruct. This include selecting slices, diffusion directions, 
and averages that we want.

First we sort the slice index that we selected to reconstruct.
=#

if isempty(sliceChoice) || !(@isdefined sliceChoice)
    sliceChoice = collect(1:nSlices)
end

isMultiSlice = length(sliceChoice) > 1

if !isMultiSlice
    selectedSlice = sliceChoice
else
    selectedSlice = sort(vec(sliceChoice))
end

#=
Next we select the data we would like to reconstruct from the ISMRMRD file. 

The ISMRMRD data are stored in the following loops:

Slice 1, Slice 2 ... Slice N   Slice 1, Slice 2 ... Slice N     Slice 1, Slice 2 ... Slice N ... 

|______ Diff Dir 1 ______|   |______ Diff Dir 2 ______| ... |______ Diff Dir N ______| ... 

|_________________________________ Average 1 ___________________________________| ... |___ Average N___| 

Here we chose the set corresponding to the b-value = 0 images under the first average as the example.
There is a constant shift due to pre-scan data that we want to skip, which is why the data starts from `nSlices*2`.
=#
excitationList = collect(nSlices*2 + 2 : 2 : nSlices*4) .+ diffusionDirection * nSlices * 2 .+ (idxAverage - 1) * nSlices * (nTotalDiffDir + 1) * 2
sliceSelection = excitationList[selectedSlice]

#=
#### 3.2.2 Synchronization and Merging of k-space Data and Trajectory

Since the k-space data and spiral k-space trajectories are sampled under different sampling rates 
and stored in separate files, they need to be first synchronized into the frequency of k-space data 
and then merged into a single object before final spiral image reconstruction.

Here we use a dictionary `paramsSpiral` to hold the parameters for this k-space data/trajectory synchronization and merging.
=#

paramsSpiral = Dict{Symbol,Any}()
paramsSpiral[:reconSize] = Tuple(paramsGeneral[:reconSize])
paramsSpiral[:interleave] = startIndexIntlv
paramsSpiral[:numSamples] = paramsGeneral[:numADCSamples]
paramsSpiral[:delay] = 0.00000 # naive delay correction

paramsSpiral[:interleaveDataFileNames] = paramsGeneral[:fullPathScan]

paramsSpiral[:trajFilename] = paramsGeneral[:fullPathGradient]
paramsSpiral[:excitations] = sliceSelection

paramsSpiral[:doMultiInterleave] = !isDataSingleIntlv
paramsSpiral[:doOddInterleave] = false
paramsSpiral[:numInterleaves] = isDataSingleIntlv ? 1 : length(paramsSpiral[:interleaveDataFileNames]) # one interleaf per file, count files, if filenames are array of strings (not only one string)

paramsSpiral[:singleSlice] = !isMultiSlice


#=
Here we synchronize the spiral k-space data with trajectory by upsampling the trajectory. 
Subsequently, data of all the selected spiral interleaves and the corresponding trajectories 
are merged into `acqDataImaging`. 
This step is done through the function `mergeRawInterleaves`, which can be viewed in 
[Utils.jl](@__REPO_ROOT_URL__/utils/Utils.jl).

Note that we only do these steps when they have not been done yet or it's specifically required.
=#
if reloadSpiralData || !(@isdefined acqDataImaging)
    @info "Reading spiral data and merging interleaves"
    acqDataImaging = mergeRawInterleaves(paramsSpiral)
end

# for i = 1:20

#     acqDataImaging.kdata[1][:,i] .= vcat(acqDataImaging.kdata[1][4:end,i],zeros(ComplexF32,(3,1)))

# end

#=
#### 3.2.3 Correction of k-space Trajectory Using Gradient Impulse Response Function

The previously calculated GIRFs are loaded. 
The spiral trajectory is corrected by the 1st and 0th order of GIRF.

Finally we check if the k-space trajectory is normalized to the range [-0.5, 0.5].
=#

#Load GIRFs (K1)
gK1 = readGIRFFile(paramsGeneral[:fullPathGIRF][1], paramsGeneral[:fullPathGIRF][2], paramsGeneral[:fullPathGIRF][3], "GIRF_FT", false)
gAk1 = GirfApplier(gK1, paramsGeneral[:gamma])
#Load K0 GIRF
gK0 = readGIRFFile(paramsGeneral[:fullPathGIRF][1], paramsGeneral[:fullPathGIRF][2], paramsGeneral[:fullPathGIRF][3], "b0ec_FT", true)
gAk0 = GirfApplier(gK0, paramsGeneral[:gamma])

if paramsGeneral[:doCorrectWithGIRFkxyz] 
    @info "Correcting For GIRF"
    applyGIRF!(acqDataImaging, gAk1)
end

if paramsGeneral[:doCorrectWithGIRFk0]
    @info "Correcting For k₀"
    applyK0!(acqDataImaging, gAk0)
end

#Check the k-space nodes so they don't exceed frequency limits [-0.5, 0.5] (inclusive)
checkAcquisitionNodes!(acqDataImaging)

#=
#### 3.2.4 Center the Object to the Field-of-View (FOV)

If the scanned object is not in the center of the FOV, we need to shift FOV 
to place the object in the center. This is achieved through adding linear phases 
on all dimensions.
=#
shiftksp!(acqDataImaging,paramsGeneral[:fovShift])

#=
#### 3.2.5 Processing Coil Sensitivity Maps

We need to preprocess the coil sensitivity maps before reconstruction. 
This includes resizing the coil maps to the size of output encoding matrix size; 
compress the channels according to user's setting to achieve a faster reconstruction.
=#
low_freq_mask = hanning((30,30),padding=170,zerophase=true)

sensitivity = mapslices(x ->imresize(x, paramsSpiral[:reconSize][1],paramsSpiral[:reconSize][2]), senseCartesian, dims=[1,2])

#Optional: Plot the sensitivity maps of each coil on a given slice.
if paramsGeneral[:doPlotRecon]
    plotlyjs()
    plotSenseMaps(sensitivity,size(sensitivity, 4),sliceIndex = 2)
end

#Do coil compression to make recon faster
if paramsGeneral[:doCoilCompression]
    acqDataImaging, sensitivity = geometricCC_2d(acqDataImaging,sensitivity, paramsGeneral[:nVirtualCoils])
end

#=
#### 3.2.6 Processing Off-Resonance (B0) Maps

We need to resize the B0 maps to the size of output encoding matrix size.
=#
resizedB0 = mapslices(x->imresize(x,paramsSpiral[:reconSize][1],paramsSpiral[:reconSize][2]), b0Maps, dims=[1,2])

#=
#### 3.2.7 Alignment of Off-Resonance, Sensitivity, and Spiral Data

We need to make sure that the axes line up so we rotate the sensitivities and the off-resonance maps  
Depending on your geometry, this might not be necessary but it is here
=#
# resizedB0 = mapslices(x->rotl90(x),resizedB0,dims=[1,2])
# sensitivity = mapslices(x->rotl90(x),sensitivity,dims=[1,2])

#=
### 3.3 Spiral Image Reconstruction

Here we start the spiral image reconstruction.

First we need to set necessary parameters for reconstruction, 
including iterative solver's setting, coil maps, B0 maps, etc. 
These parameters are held under the dictionary `paramsRecon`.

Note that it is safer to cast B0 maps to ComplexF32 if the current version of MRIReco.jl is used.
=#

@info "Setting Reconstruction Parameters"
paramsRecon = Dict{Symbol,Any}()
paramsRecon[:reco] = "multiCoil"
paramsRecon[:reconSize] = paramsSpiral[:reconSize][1:2]
paramsRecon[:regularization] = "L2"
paramsRecon[:λ] = 1e-3
paramsRecon[:iterations] = paramsGeneral[:nReconIterations]
paramsRecon[:solver] = "cgnr"
paramsRecon[:solverInfo] = SolverInfo(ComplexF32,store_solutions=false)
paramsRecon[:senseMaps] = ComplexF32.(sensitivity[:,:,selectedSlice,:])

if paramsGeneral[:doCorrectWithB0map]
    paramsRecon[:correctionMap] = ComplexF32.(-1im.*resizedB0[:,:,selectedSlice])
end

#= 
Finally we can call reconstruction function of the package `MRIReco.jl` 
to perform final spiral image reconstruction.
=#
@info "Performing Spiral Reconstruction"
@time reco = reconstruction(acqDataImaging, paramsRecon)

GC.gc()

#=
## 4. Save and Plot the Results (Optional)

All results could be saved into NIfTI files using the `saveMap` function 
and be plotted using the `plotReconstruction` function, both located in 
the file [Utils.jl](@__REPO_ROOT_URL__/utils/Utils.jl).

=#
if paramsGeneral[:doSaveRecon] #TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    resolution_tmp = fieldOfView(acqDataImaging)[1:2]./encodingSize(acqDataImaging)
    resolution_mm = (resolution_tmp[1],resolution_tmp[2],fieldOfView(acqDataImaging)[3] *(1 + paramsGeneral[:sliceDistanceFactor_percent]/100.0)) #for 2D only, since FOV[3] is slice thickness then, but gap has to be observed

    #TODO: use slice ordering from cartesian scan directly!
    nSlices = numSlices(acqDataImaging)
    sliceIndexArray = getSliceOrder(nSlices, isSliceInterleaved = true)
    saveMap(paramsGeneral[:fullPathSaveRecon], paramsGeneral[:scalingFactorSaveRecon]*reco.data[:,:,sliceIndexArray], resolution_mm; doSplitPhase=true, doNormalize = paramsGeneral[:doNormalizeRecon])
end

if paramsGeneral[:doPlotRecon]
    @info "Plotting Reconstruction"
    plotlyjs()
    plotReconstruction(reco, 1:length(selectedSlice), resizedB0[:, :, selectedSlice], figHandles=["Original Magnitude", "Original Phase", "B0"], isSliceInterleaved=false, rotateAngle=90)
end

@info "Successfully Completed SpiralRecon"