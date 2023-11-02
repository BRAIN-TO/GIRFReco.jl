using HDF5, MRIReco, LinearAlgebra, Dierckx, DSP, FourierTools, ImageBinarization, ImageEdgeDetection, MRIGradients

# %%
# Include tools and reader functions for running the spiral reconstruction recipe
# Note: the files are found relative of the location of the folder, not the
# environment current folder
include("../io/GradientReader.jl")
include("../utils/Utils.jl")

## Set true if we need to reload Cartesian and/or spiral data compulsively.
reloadCartesianData = true
reloadSpiralData = true
doCoilCompression = true

## Gyromagnetic ratio, in unit of Hz
gamma = 42577478

## Only calculate sensitivity and B0 maps when they have not been done yet, or it's specifically required.
if reloadCartesianData || !((@isdefined cartesian_sensitivity) && (@isdefined b0Maps))
    ## Executing Cartesian recon from which B0/sensitivity maps have been computed
    @info "Running CartesianRecon to retrieve maps (cartesian_sensitivity and b0Maps)"
    include("CartesianRecon_Mar2022_Human.jl")
end

## Set figures to be unlocked from the win9ow (i.e use matplotlib backend with controls)
#pygui(true)

## Choose Slice (can be [single number] OR [1,2,3,4,5,6,7,8,9]
sliceChoice = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # UNCOMMENT FOR MULTISLICE
# sliceChoice = [6] # UNCOMMENT FOR SINGLESLICE (SLICES 3, 7 and 8 are good examples)
diffusionDirection = 0 # CAN BE FROM 0 (b=0) to 6 (e.g. for 6 direction MDDW, 1-6 are 6 directions)

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
# excitationList = vec(20:2:36) .+ diffusionDirection * 9 * 2 # DATASET SPECIFIC INDEXING
excitationList = vec(32:2:62) .+ diffusionDirection * 15 * 2 # DATASET SPECIFIC INDEXING: 15 slices, starting from profile 32
sliceSelection = excitationList[selectedSlice]

@info "Slice Chosen = $selectedSlice: \n \nExcitations Chosen = $excitationList "

fname_spiralIntlv1 = "data/Spirals/508_124_2.h5"
fname_spiralIntlv2 = "data/Spirals/508_126_2.h5"
fname_spiralIntlv3 = "data/Spirals/508_128_2.h5"
fname_spiralIntlv4 = "data/Spirals/508_130_2.h5"
fname_gradient = "data/Gradients/gradients508.txt"
fname_girfGx = "data/GIRF/GIRF_ISMRM2022/2021Nov_PosNeg_Gx.mat"
fname_girfGy = "data/GIRF/GIRF_ISMRM2022/2021Nov_PosNeg_Gy.mat"
fname_girfGz = "data/GIRF/GIRF_ISMRM2022/2021Nov_PosNeg_Gz.mat"

# adjustmentDict is the dictionary that sets the information for correct data loading and trajectory and data synchronization
adjustmentDict = Dict{Symbol,Any}()
adjustmentDict[:reconSize] = (200, 200)
adjustmentDict[:interleave] = 1
adjustmentDict[:slices] = 1
adjustmentDict[:coils] = 20
# adjustmentDict[:numSamples] = 16084 # Total Number of ADC event, including the period of gradient rewinder
adjustmentDict[:numSamples] = 15504 # Total Number of readout before gradient rewinder
adjustmentDict[:delay] = 0.00000 # naive delay correction

adjustmentDict[:interleaveDataFileNames] = [fname_spiralIntlv1, fname_spiralIntlv2, fname_spiralIntlv3, fname_spiralIntlv4]
adjustmentDict[:trajFilename] = fname_gradient
adjustmentDict[:excitations] = sliceSelection

adjustmentDict[:doMultiInterleave] = true
adjustmentDict[:doOddInterleave] = true
adjustmentDict[:numInterleaves] = 4

adjustmentDict[:singleSlice] = !multiSlice

# Defined recon size and parameters for data loading
@info "Using Parameters:\n\nreconSize = $(adjustmentDict[:reconSize]) \n interleave = $(adjustmentDict[:interleave]) \n slices = $(adjustmentDict[:slices]) \n coils = $(adjustmentDict[:coils]) \n numSamples = $(adjustmentDict[:numSamples])\n\n"

## Only load data when it has not been done yet, or it's specifically required.
if reloadSpiralData || !(@isdefined acqDataImaging)
    ## Convert raw to AcquisitionData

    @info "Merging interleaves and reading data \n"
    acqDataImaging = merge_raw_interleaves(adjustmentDict)

    @info "Loading Gradient Impulse Response Functions \n"
    ## Load GIRFs!
    # Tim Wu, use new read GIRF function
    #gK1 = loadGirf(1,1)
    gK1 = readGIRFFile(fname_girfGx, fname_girfGy, fname_girfGz, "GIRF_FT", false)
    gAk1 = GirfApplier(gK1, gamma)

    @info "Correcting For GIRF \n"
    apply_girf!(acqDataImaging, gAk1)

    # Load K₀ GIRF
    # Tim Wu, use new read GIRF function
    #gK0 = loadGirf(0,1)
    gK0 = readGIRFFile(fname_girfGx, fname_girfGy, fname_girfGz, "b0ec_FT", true)
    gAk0 = GirfApplier(gK0, gamma)

    @info "Correcting For k₀ \n"
    apply_k0!(acqDataImaging, gAk0)

    ## Check the k-space nodes so they don't exceed frequency limits [-0.5, 0.5] (inclusive)
    check_acquisition_nodes!(acqDataImaging)

end

## Sense Map loading
@info "Recalculating Sense Maps \n"
sensitivity = espirit(acqDataCartesian, (4, 4), 12, adjustmentDict[:reconSize], eigThresh_1 = 0.01, eigThresh_2 = 0.98)

# shift FOV to middle :) 
shiftksp!(acqDataImaging, [0, -20])
#changeFOV!(acqDataImaging,[1.5,1.5])

nvcoils = size(sensitivity, 4)

doCoilCompression = false

## Do coil compression to make recon faster
if doCoilCompression
    nvcoils = 4
    acqDataImaging, sensitivity = geometricCC_2d(acqDataImaging, sensitivity, nvcoils)
end

# ## Plot the sensitivity maps of each coil
if params_general[:do_plot_recon]
    @info "Plotting SENSE Maps"
    plot_sense_maps(sensitivity, nvcoils, sliceIndex = 10)
end

## B0 Maps (Assumes we have a B0 map from gradient echo scan named b0)
@info "Resizing B0 Maps"
resizedB0 = mapslices(x -> imresize(x, adjustmentDict[:reconSize]), b0Maps, dims = [1, 2])

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
params[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false)
params[:senseMaps] = ComplexF32.(sensitivity[:, :, selectedSlice, :])
params[:correctionMap] = ComplexF32.(-1im .* resizedB0[:, :, selectedSlice])

## Call to reconstruction
@info "Performing Reconstruction \n"
@time reco = reconstruction(acqDataImaging, params)

#totalRecon = sum(abs2,reco.data,dims=5)
@info "Plotting Reconstruction \n"
plot_reconstruction(cartesianReco, 1:length(selectedSlice), resizedB0[:, :, selectedSlice], isSliceInterleaved = true, rotateAngle = 270)

## Plot the image edges (feature comparison)

# img_edges₁ = detect_edges(slice1,Canny(spatial_scale = 2.6))
# img_edges₂ = detect_edges(slice2,Canny(spatial_scale = 2.7))

# imEdges = cat(img_edges₁,img_edges₂,zeros(size(img_edges₁)),dims=3)

# figure("Edge Differences")
# imshow(imEdges)

@info "Successfully Completed SpiralRecon \n"
