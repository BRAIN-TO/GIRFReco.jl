using HDF5, MRIReco, LinearAlgebra, Dierckx, DSP, FourierTools, ImageBinarization, ImageEdgeDetection, MRIGradients

##
# Include tools and reader functions for running the spiral reconstruction recipe
# Note: the files are found relative of the location of the folder, not the
# environment current folder
include("../io/GradientReader.jl")
include("../utils/Utils.jl")

## ----------------------------- User-defined Variables -------------------------- ##

# All data-specific recon parameters
include("ReconConfig.jl")

## Set true if we need to reload Cartesian and/or spiral data compulsively.
reloadCartesianData = true
reloadSpiralData = true
reloadGIRFData = true

## Choose Slice (can be [single number] OR [1,2,3,...])
sliceChoice = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] # For multi-slice
# sliceChoice = [6] # For single-slice

## Choose diffusion direction; starting from 0 (b=0) to the total number in MDDW protocol, e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs)
diffusionDirection = 0
idxAverage = 4;
nDiffusionDirections = 6;

## Determine to reconstruct single-interleave data, or one interleave out of multi-interleave data.
isDataSingleIntlv = isa(paramsGeneral[:fullPathScan], String)

# Which interleave to be reconstructed. For single-interleave data, it will always be set as 1; for multi-interleave data, the value set here will be used.
# For multi-interleaved data, this value is ranging from [1:TotNumIntlv] (total number of interleaves), indicating which interleave to be reconstructed
startIndexIntlv = 1

## Total number of ADC points BEFORE the rewinder at the end of the spiral readout. For gradient 508, use 15655 (out of 16084); for gradient 511, use 15445 (out of 15624).
numADCSamples = 15504
# numADCSamples = 15655
# numADCSamples = 15445

## Gyromagnetic ratio, in unit of Hz
gamma = 42577478

## ------------------------------------------------ Calculation Starts Here ---------------------------------------------------------- ##

## Only calculate sensitivity and B0 maps when they have not been done yet, or it's specifically required.
if reloadCartesianData || !((@isdefined senseCartesian) && (@isdefined b0Maps))
    ## Executing Cartesian recon from which B0/sensitivity maps have been computed
    @info "Running CartesianRecon to retrieve maps (senseCartesian and b0Maps)"
    include("CartesianRecon_Mar2022_Human.jl")
end

totalSliceNum = size(sensitivity, 3)

## Set figures to be unlocked from the win9ow (i.e use matplotlib backend with controls)

## Spiral Reconstruction Recipe Starts Here
@info "Starting Spiral Reconstruction Pipeline"

## Default to single slice selection. Choose multi-slice only if computer is capable.
multiSlice = true

if length(sliceChoice) > 1
    multiSlice = true
end

if !multiSlice
    selectedSlice = sliceChoice
else
    selectedSlice = sort(vec(sliceChoice))
end

## The ISMRMRD File contains more than one excitation, so we choose the set corresponding to the b-value 0 images
excitationList = collect(totalSliceNum*2 + 2 : 2 : totalSliceNum*4) .+ diffusionDirection * totalSliceNum * 2 .+ (idxAverage - 1) * totalSliceNum * (nDiffusionDirections + 1) * 2 # DATASET SPECIFIC INDEXING: 15 slices, starting from profile 32
sliceSelection = excitationList[selectedSlice]

@info "Slice Chosen = $selectedSlice: \n \nExcitations Chosen = $excitationList "

# adjustmentDict is the dictionary that sets the information for correct data loading and trajectory and data synchronization
adjustmentDict = Dict{Symbol,Any}()
adjustmentDict[:reconSize] = paramsGeneral[:reconSize]
adjustmentDict[:interleave] = startIndexIntlv
adjustmentDict[:slices] = 1
adjustmentDict[:numSamples] = numADCSamples
adjustmentDict[:delay] = 0.00000 # naive delay correction

adjustmentDict[:interleaveDataFileNames] = paramsGeneral[:fullPathScan]

adjustmentDict[:trajFilename] = paramsGeneral[:fullPathGradient]
adjustmentDict[:excitations] = sliceSelection

adjustmentDict[:doMultiInterleave] = !isDataSingleIntlv
adjustmentDict[:doOddInterleave] = false
adjustmentDict[:numInterleaves] = isDataSingleIntlv ? 1 : length(adjustmentDict[:interleaveDataFileNames]) # one interleaf per file, count files, if filenames are array of strings (not only one string)

adjustmentDict[:singleSlice] = !multiSlice

# Defined recon size and parameters for data loading
@info "Using Parameters:\n\nreconSize = $(adjustmentDict[:reconSize]) \n interleave = $(adjustmentDict[:interleave]) \n slices = $(adjustmentDict[:slices]) \n coils = $(size(sensitivity, 4)) \n numSamples = $(adjustmentDict[:numSamples])\n\n"

if reloadGIRFData || !(@isdefined gK1) || !(@isdefined gAK1) || !(@isdefined gK0) || !(@isdefined gAK0)
    @info "Loading Gradient Impulse Response Functions"
    
    ## Load GIRFs (K1)
    gK1 = readGIRFFile(paramsGeneral[:fullPathGIRF][1], paramsGeneral[:fullPathGIRF][2], paramsGeneral[:fullPathGIRF][3], "GIRF_FT", false)
    gAk1 = GirfApplier(gK1, gamma)

    ## Load K₀ GIRF
    gK0 = readGIRFFile(paramsGeneral[:fullPathGIRF][1], paramsGeneral[:fullPathGIRF][2], paramsGeneral[:fullPathGIRF][3], "b0ec_FT", true)
    gAk0 = GirfApplier(gK0, gamma)
end

## Only load data when it has not been done yet, or it's specifically required.
if reloadSpiralData || !(@isdefined acqDataImaging)
    ## Convert raw to AcquisitionData

    @info "Reading spiral data and merging interleaves"
    acqDataImaging = mergeRawInterleaves(adjustmentDict)

    @info "Correcting For GIRF"
    applyGIRF!(acqDataImaging, gAk1)

    @info "Correcting For k₀"
    applyK0!(acqDataImaging, gAk0)

    ## Check the k-space nodes so they don't exceed frequency limits [-0.5, 0.5] (inclusive)
    checkAcquisitionNodes!(acqDataImaging)

end

## Sense Map loading
@info "Resizing Sense Maps"

# Resize sense maps to match encoding size of data matrix
sensitivity = mapslices(x ->imresize(x, adjustmentDict[:reconSize]), senseCartesian, dims=[1,2])

# Plot the sensitivity maps of each coil
@info "Plotting SENSE Maps"

if paramsGeneral[:doPlotRecon]
    plotSenseMaps(sensitivity,size(sensitivity, 4),sliceIndex = 10)
end


# shift FOV to middle :) 
shiftksp!(acqDataImaging,paramsGeneral[:fovShift])
#changeFOV!(acqDataImaging,[1.5,1.5])


## Do coil compression to make recon faster
if paramsGeneral[:doCoilCompression]
    acqDataImaging, sensitivity = geometricCC_2d(acqDataImaging,sensitivity, paramsGeneral[:nVirtualCoils])
end


## B0 Maps (Assumes we have a B0 map from gradient echo scan named b0)
@info "Resizing B0 Maps"
resizedB0 = mapslices(x->imresize(x,adjustmentDict[:reconSize]), b0Maps, dims=[1,2])

## Define Parameter Dictionary for use with reconstruction
# CAST TO ComplexF32 if you're using current MRIReco.jl

@info "Setting Parameters"
params = Dict{Symbol,Any}()
params[:reco] = "multiCoil"
params[:reconSize] = adjustmentDict[:reconSize]
params[:regularization] = "L2"
params[:λ] = 1e-2 # CHANGE THIS TO GET BETTER OR WORSE RECONSTRUCTION RESULTS
params[:iterations] = 20
params[:solver] = "cgnr"
params[:solverInfo] = SolverInfo(ComplexF32,store_solutions=false)
params[:senseMaps] = ComplexF32.(sensitivity[:,:,selectedSlice,:])
params[:correctionMap] = ComplexF32.(-1im.*resizedB0[:,:,selectedSlice])

## Call to reconstruction
@info "Performing Reconstruction"
@time reco = reconstruction(acqDataImaging,params)
#totalRecon = sum(abs2,reco.data,dims=5)

# save Map recon (multi-echo etc.)
if paramsGeneral[:doSaveRecon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    resolution_mm = fieldOfView(acqDataImaging)./encodingSize(acqDataImaging)
    resolution_mm[3] = fieldOfView(acqDataImaging)[3] *(1 + paramsGeneral[:sliceDistanceFactor_percent]/100.0); # for 2D only, since FOV[3] is slice thickness then, but gap has to be observed

    # TODO: use slice ordering from cartesian scan directly!
    nSlices = numSlices(acqDataImaging)
    sliceIndexArray = getSliceOrder(nSlices, isSliceInterleaved = true)
    saveMap(paramsGeneral[:fullPathSaveRecon], reco.data[:,:,sliceIndexArray], resolution_mm; doSplitPhase=true, doNormalize = true)
end

if paramsGeneral[:doPlotRecon]
    @info "Plotting Reconstruction"
    pygui(true)
    plotReconstruction(reco, 1:length(selectedSlice), resizedB0[:, :, selectedSlice], figHandles=["Original Magnitude", "Original Phase", "B0"], isSliceInterleaved=true, rotateAngle=270)
end

@info "Successfully Completed SpiralRecon"