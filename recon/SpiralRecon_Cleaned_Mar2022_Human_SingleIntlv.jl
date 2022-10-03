using HDF5, MRIReco, LinearAlgebra, Dierckx, DSP, FourierTools, ImageBinarization, ImageEdgeDetection, MRIGradients

# %%
# Include tools and reader functions for running the spiral reconstruction recipe
# Note: the files are found relative of the location of the folder, not the
# environment current folder
include("../io/GradientReader.jl")
include("../utils/Utils.jl")

## ----------------------------- User-defined Variables -------------------------- ##

## Set true if we need to reload Cartesian and/or spiral data compulsively.
reloadCartesianData = false
reloadSpiralData = true
reloadGIRFData = false

## Choose Slice (can be [single number] OR [1,2,3,...])
sliceChoice = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] # For multi-slice
# sliceChoice = [6] # For single-slice

## Matrix size of the reconstructed image. For gradient 508 with all 4 interleaves, use 200 for high resolution image; otherwise consider using 112 or 84 for a lower resolution. The FOV is 220 mm for both gradients 508 and 511.
reconSize = (112,112) #(200, 200) for gradient 508

## Choose diffusion direction; starting from 0 (b=0) to the total number in MDDW protocol, e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs)
diffusionDirection = 0

## Determine to reconstruct single-interleave data, or one interleave out of multi-interleave data.
isDataSingleIntlv = true

# Which interleave to be reconstructed. For single-interleave data, it will always be set as 1; for multi-interleave data, the value set here will be used.
# For multi-interleaved data, this value is ranging from [1:TotNumIntlv] (total number of interleaves), indicating which interleave to be reconstructed
startIndexIntlv = 3

# For single interleave data, use this section
if isDataSingleIntlv
    startIndexIntlv = 1 # Should always be 1 for single-interleave data.
    # fname_spiralIntlv = "data/Spirals/511_134_2.h5" # Gradient 511, b = 300, 10 diff directions
    # fname_spiralIntlv = "data/Spirals/511_136_2.h5" # Gradient 511, b = 700, 30 diff directions
    fname_spiralIntlv = "D:\\OneDrive - UHN\\MRP-SPIDI\\SPIDI\\data\\SPIDI_0007\\Human\\dat\\511_138_2.h5" # Gradient 511, b = 2500, 64 diff directions
    # fname_spiralIntlv = "data/Spirals/508_140_2.h5" # Gradient 508, interleave 0, b = 300, 10 diff directions
    # fname_spiralIntlv = "data/Spirals/508_142_2.h5" # Gradient 508, interleave 0, b = 700, 30 diff directions
    # fname_spiralIntlv = "data/Spirals/508_144_2.h5" # Gradient 508, interleave 0, b = 2500, 64 diff directions
else
    # Multi-interleave data, needs all 4 file names, but will only read the corresponding one.
    fname_spiralIntlv0 = "D:\\OneDrive - UHN\\MRP-SPIDI\\SPIDI\\data\\SPIDI_0007\\Human\\dat\\508_124_2.h5" # Gradient 508, interleave 0, b = 2000, 6 diff directions, 4 averages
    fname_spiralIntlv1 = "D:\\OneDrive - UHN\\MRP-SPIDI\\SPIDI\\data\\SPIDI_0007\\Human\\dat\\508_126_2.h5" # Gradient 508, interleave 1, b = 2000, 6 diff directions, 4 averages
    fname_spiralIntlv2 = "D:\\OneDrive - UHN\\MRP-SPIDI\\SPIDI\\data\\SPIDI_0007\\Human\\dat\\508_128_2.h5" # Gradient 508, interleave 2, b = 2000, 6 diff directions, 4 averages
    fname_spiralIntlv3 = "D:\\OneDrive - UHN\\MRP-SPIDI\\SPIDI\\data\\SPIDI_0007\\Human\\dat\\508_130_2.h5" # Gradient 508, interleave 3, b = 2000, 6 diff directions, 4 averages
end

## Total number of ADC points BEFORE the rewinder at the end of the spiral readout. For gradient 508, use 15655 (out of 16084); for gradient 511, use 15445 (out of 15624).
# numADCSamples = 15655
numADCSamples = 15445

## File name for the spiral gradient
# fname_gradient = "D:\\OneDrive - UHN\\MRP-SPIDI\\SPIDI\\data\\SPIDI_0007\\508\\gradients.txt" # Contains all 4 interleaves.
fname_gradient = "D:\\OneDrive - UHN\\MRP-SPIDI\\SPIDI\\data\\SPIDI_0007\\511\\gradients.txt"

fname_girfGx = "D:\\SpiralDiffusion\\DataNov2020\\GIRF\\GIRF_ISMRM2022\\2021Nov_PosNeg_Gx.mat"
fname_girfGy = "D:\\SpiralDiffusion\\DataNov2020\\GIRF\\GIRF_ISMRM2022\\2021Nov_PosNeg_Gy.mat"
fname_girfGz = "D:\\SpiralDiffusion\\DataNov2020\\GIRF\\GIRF_ISMRM2022\\2021Nov_PosNeg_Gz.mat"

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
multiSlice = false

if length(sliceChoice) > 1
    multiSlice = true
end

if !multiSlice
    selectedSlice = sliceChoice
else
    selectedSlice = sort(vec(sliceChoice))
end

## The ISMRMRD File contains more than one excitation, so we choose the set corresponding to the b-value 0 images
excitationList = collect(totalSliceNum*2 + 2 : 2 : totalSliceNum*4) .+ diffusionDirection * totalSliceNum * 2 # DATASET SPECIFIC INDEXING: 15 slices, starting from profile 32
sliceSelection = excitationList[selectedSlice]

@info "Slice Chosen = $selectedSlice: \n \nExcitations Chosen = $excitationList "

# adjustmentDict is the dictionary that sets the information for correct data loading and trajectory and data synchronization
adjustmentDict = Dict{Symbol,Any}()
adjustmentDict[:reconSize] = reconSize
adjustmentDict[:interleave] = startIndexIntlv
adjustmentDict[:slices] = 1
adjustmentDict[:numSamples] = numADCSamples
adjustmentDict[:delay] = 0.00000 # naive delay correction

if isDataSingleIntlv
    adjustmentDict[:interleaveDataFileNames] = [fname_spiralIntlv]
else
    adjustmentDict[:interleaveDataFileNames] = [fname_spiralIntlv0, fname_spiralIntlv1, fname_spiralIntlv2, fname_spiralIntlv3]
end

adjustmentDict[:trajFilename] = fname_gradient
adjustmentDict[:excitations] = sliceSelection

adjustmentDict[:doMultiInterleave] = false
adjustmentDict[:doOddInterleave] = false
adjustmentDict[:numInterleaves] = 1

adjustmentDict[:singleSlice] = !multiSlice

# Defined recon size and parameters for data loading
@info "Using Parameters:\n\nreconSize = $(adjustmentDict[:reconSize]) \n interleave = $(adjustmentDict[:interleave]) \n slices = $(adjustmentDict[:slices]) \n coils = $(size(sensitivity, 4)) \n numSamples = $(adjustmentDict[:numSamples])\n\n"

if reloadGIRFData || !(@isdefined gK1) || !(@isdefined gAK1) || !(@isdefined gK0) || !(@isdefined gAK0)
    @info "Loading Gradient Impulse Response Functions \n"
    
    ## Load GIRFs (K1)
    gK1 = readGIRFFile(fname_girfGx, fname_girfGy, fname_girfGz, "GIRF_FT", false)
    gAk1 = GirfApplier(gK1, gamma)

    ## Load K₀ GIRF
    gK0 = readGIRFFile(fname_girfGx, fname_girfGy, fname_girfGz, "b0ec_FT", true)
    gAk0 = GirfApplier(gK0, gamma)
end

## Only load data when it has not been done yet, or it's specifically required.
if reloadSpiralData || !(@isdefined acqDataImaging)
    ## Convert raw to AcquisitionData

    @info "Reading spial data and merging interleaves \n"
    acqDataImaging = mergeRawInterleaves(adjustmentDict)

    @info "Correcting For GIRF \n"
    applyGIRF!(acqDataImaging, gAk1)

    @info "Correcting For k₀ \n"
    applyK0!(acqDataImaging, gAk0)

    ## Check the k-space nodes so they don't exceed frequency limits [-0.5, 0.5] (inclusive)
    checkAcquisitionNodes!(acqDataImaging)

end

## Sense Map loading
@info "Validating Sense Maps \n"

# Resize sense maps to match encoding size of data matrix
sensitivity = mapslices(x ->imresize(x, (acqDataImaging.encodingSize[1],acqDataImaging.encodingSize[2])), senseCartesian, dims=[1,2])

# ## Plot the sensitivity maps of each coil
@info "Plotting SENSE Maps \n"
plotSenseMaps(sensitivity,size(sensitivity, 4),sliceIndex = 10)

## B0 Maps (Assumes we have a B0 map from gradient echo scan named b0)
@info "Resizing B0 Maps \n"
resizedB0 = mapslices(x->imresize(x,(acqDataImaging.encodingSize[1], acqDataImaging.encodingSize[2])), b0Maps2, dims=[1,2])

## Define Parameter Dictionary for use with reconstruction
# CAST TO ComplexF32 if you're using current MRIReco.jl

@info "Setting Parameters \n"
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
@info "Performing Reconstruction \n"
@time reco = reconstruction(acqDataImaging,params)

#totalRecon = sum(abs2,reco.data,dims=5)
@info "Plotting Reconstruction \n"
pygui(true)
plotReconstruction(reco, 1:length(selectedSlice), resizedB0[:,:,selectedSlice], figHandles = ["Original Magnitude", "Original Phase", "B0"], isSliceInterleaved = true, rotateAngle = 270)

@info "Successfully Completed SpiralRecon \n"