using PyPlot, HDF5, MRIReco, LinearAlgebra, Dierckx, DSP, FourierTools, ImageBinarization, ImageEdgeDetection, MRIGradients

# Note: the files are found relative of the location of the folder, not the
# environment current folder
include("../io/GradientReader.jl")

include("../utils/Utils.jl")

## Executing Cartesian recon from which B0/sensitivity maps have been computed
@info "Running julia_recon_cartesian to retrieve maps (senseCartesian and b0Maps)"
include("../recon/CartesianRecon.jl")


## Preparation
pygui(true)

@info "Starting Spiral Reconstruction"

## Load ISMRMRD data files (can be undersampled) THIS SHOULD BE THE ONLY SECTION NEEDED TO EDIT TO ADJUST FOR DIFFERENT SCANS
# @info "Loading Data Files"

# selectedSlice = 3

selectedSlice = vec(1:9)
# excitationList = [4]

excitationList = vec(20:2:36) # for MULTISLICE
sliceSelection = excitationList[selectedSlice]

@info "Slice Chosen = $selectedSlice: \n \nExcitations Chosen = $excitationList "

# adjustmentDict is the dictionary that sets the information for correct data loading and trajectory and data synchronization
adjustmentDict = Dict{Symbol,Any}()
adjustmentDict[:reconSize] = (200,200)
adjustmentDict[:interleave] = 1
adjustmentDict[:slices] = 1
adjustmentDict[:coils] = 5
adjustmentDict[:numSamples] = 15475
adjustmentDict[:delay] = 0.00000 # naive delay correction

# adjustmentDict[:interleaveDataFileNames] = ["data/Spirals/523_21_1_2.h5", "data/Spirals/523_23_2_2.h5", "data/Spirals/523_25_3_2.h5", "data/Spirals/523_27_4_2.h5"]

adjustmentDict[:interleaveDataFileNames] = ["data/Spirals/523_96_2.h5","data/Spirals/523_98_2.h5", "data/Spirals/523_100_2.h5", "data/Spirals/523_102_2.h5"]
adjustmentDict[:trajFilename] = "data/Gradients/gradients523.txt"
adjustmentDict[:excitations] = sliceSelection

adjustmentDict[:doMultiInterleave] = false
adjustmentDict[:doOddInterleave] = false
adjustmentDict[:numInterleaves] = 1

adjustmentDict[:singleSlice] = false

@info "Using Parameters:\n\nreconSize = $(adjustmentDict[:reconSize]) \n interleave = $(adjustmentDict[:interleave]) \n slices = $(adjustmentDict[:slices]) \n coils = $(adjustmentDict[:coils]) \n numSamples = $(adjustmentDict[:numSamples])\n\n"
# define recon size and parameters for data loading

# print(" reconSize = $(adjustmentDict[:reconSize]) \n interleave = $(adjustmentDict[:interleave]) \n slices = $(adjustmentDict[:slices]) \n coils = $(adjustmentDict[:coils]) \n numSamples = $(adjustmentDict[:numSamples])\n\n")

## Convert raw to AcquisitionData

@info "Merging interleaves and reading data \n"
acqDataImaging = mergeRawInterleaves(adjustmentDict)

@info "Loading Gradient Impulse Response Functions \n"
## Load GIRFs!
gK1 = loadGirf(1,1)
gAk1 = GirfApplier(gK1, 42577478)

gK0 = loadGirf(0,1)
gAk0 = GirfApplier(gK0, 42577478)

# @info "Correcting Coil Phase"
# calibrateAcquisitionPhase!(acqDataImaging)

@info "Correcting For GIRF \n"
applyGIRF!(acqDataImaging, gAk1)

# @info "Correcting For k₀ \n"
# applyK0!(acqDataImaging, gAk0)

checkAcquisitionNodes!(acqDataImaging)

## Sense Map Calculation

@info "Validating Sense Maps \n"
# Resize sense maps to match encoding size of data matrix
sensitivity = mapslices(x ->imresize(x, (acqDataImaging.encodingSize[1],acqDataImaging.encodingSize[2])), senseCartesian, dims=[1,2])
sensitivity = mapslices(rotl90,sensitivity,dims=[1,2])

# ## Plot the sensitivity maps of each coil
# @info "Plotting SENSE Maps \n"
plotSenseMaps(sensitivity,adjustmentDict[:coils])

## B0 Maps (Assumes we have a B0 map from gradient echo scan named b0)
@info "Validating B0 Maps \n"
resizedB0 = mapslices(x->imresize(x,(acqDataImaging.encodingSize[1], acqDataImaging.encodingSize[2])), b0Maps, dims=[1,2])

## Define Parameter Dictionary for use with reconstruction
# CAST TO ComplexF32 if you're using current MRIReco.jl

@info "Setting Parameters \n"
params = Dict{Symbol,Any}()
params[:reco] = "multiCoil"
params[:reconSize] = adjustmentDict[:reconSize]
params[:regularization] = "L2"
params[:λ] = 5.e-6
params[:iterations] = 20
params[:solver] = "cgnr"
params[:solverInfo] = SolverInfo(ComplexF32,store_solutions=false)
params[:senseMaps] = ComplexF32.(sensitivity[:,:,selectedSlice,:])
params[:correctionMap] = ComplexF32.(-1im.*resizedB0[:,:,selectedSlice])

##
@info "Performing Reconstruction \n"
reco = reconstruction(acqDataImaging,params)

## Plotting reconstruction
@info "Plotting Reconstruction \n"

# IF MULTISLICE
indexArray = [6,7,3,8,4,9]

# IF SINGLESLICE
# indexArray = 1

#totalRecon = sum(abs2,reco.data,dims=5)
plotReconstruction(reco, indexArray, resizedB0)

## Plot the image edges (feature comparison)

# img_edges₁ = detect_edges(slice1,Canny(spatial_scale = 2.6))
# img_edges₂ = detect_edges(slice2,Canny(spatial_scale = 2.7))

# imEdges = cat(img_edges₁,img_edges₂,zeros(size(img_edges₁)),dims=3)

# figure("Edge Differences")
# PyPlot.imshow(imEdges)
