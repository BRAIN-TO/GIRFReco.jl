using HDF5,
    MRIReco,
    LinearAlgebra,
    Dierckx,
    DSP,
    FourierTools,
    ImageBinarization,
    ImageEdgeDetection,
    MRIGradients,
    FileIO,
    MRIFiles,
    MRICoilSensitivities,
    RegularizedLeastSquares,
    GIRFReco,
    MosaicViews,
    Plots,
    Images

# All data-specific recon parameters
include("ReconConfig_SPIDI_0007.jl")

## ----------------------------- User-defined Variables -------------------------- ##

## Set true if we need to reload raw data compulsively.
reloadSpiralData = true
reloadGIRFData = true

# Choose Slice (can be [single number] OR [1,2,3,...])
# Leave empty ([]) or remove this line to later select all slices
sliceChoice = []; # TODO: read from ISMRMRD itself

## Gyromagnetic ratio, in unit of Hz
gamma = 42577478

## Choose diffusion direction; starting from 0 (b=0) to the total number in MDDW protocol, e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs)
diffusionDirection = selector[:dif]
idxAverage = selector[:avg]
nDiffusionDirections = params_general[:nDiffusionDirections] # TODO: Read from ISMRMRD itself

# Which interleave to be reconstructed. For single-interleave data, it will always be set as 1; for multi-interleave data, the value set here will be used.
# For multi-interleaved data, this value is ranging from [1:TotNumIntlv] (total number of interleaves), indicating which interleave to be reconstructed
startIndexIntlv = selector[:seg]

## Determine to reconstruct single-interleave data, or one interleave out of multi-interleave data.
isDataSingleIntlv = isa(params_general[:fullPathScan], String)


## ------------------------------------------------ Calculation Starts Here ---------------------------------------------------------- ##

if params_general[:do_load_maps] && isfile(params_general[:fullPathSaveB0]) # # TODO ask for sense map (but split in magn/phase)
    @info "Loading SENSE and B0 maps from $(params_general[:fullPathSaveSense]) and $(params_general[:fullPathSaveB0])"
    # load maps, permute slice, sice files have geometric slice order
    b0Maps = load_map(params_general[:fullPathSaveB0])

    nSlices = size(b0Maps, 3)
    sliceIndexArray = get_slice_order(nSlices, isSliceInterleaved = true)

    b0Maps = b0Maps[:, :, invperm(sliceIndexArray)]
    cartesian_sensitivity = load_map(params_general[:fullPathSaveSense]; doSplitPhase = true)[:, :, invperm(sliceIndexArray), :]
else
    ## Only calculate sensitivity and B0 maps when they have not been done yet, or it's specifically required.
    ## Executing Cartesian recon from which B0/sensitivity maps have been computed
    @info "Running CartesianRecon to retrieve maps (cartesian_sensitivity and b0Maps)"
    include("CartesianRecon.jl")
    nSlices = size(b0Maps, 3)
end

## Spiral Reconstruction Recipe Starts Here
@info "Starting Spiral Reconstruction Pipeline"

if isempty(sliceChoice) || !(@isdefined sliceChoice)
    sliceChoice = collect(1:nSlices)
end

isMultiSlice = length(sliceChoice) > 1

if !isMultiSlice
    selectedSlice = sliceChoice
else
    selectedSlice = sort(vec(sliceChoice))
end

## The ISMRMRD File contains more than one excitation, so we choose the set corresponding to the b-value 0 images
excitationList = collect(nSlices*2+2:2:nSlices*4) .+ diffusionDirection * nSlices * 2 .+ (idxAverage - 1) * nSlices * (nDiffusionDirections + 1) * 2 # DATASET SPECIFIC INDEXING: 15 slices, starting from profile 32
sliceSelection = excitationList[selectedSlice]

@info "Slice Chosen = $selectedSlice: \n \nExcitations Chosen = $excitationList "

# adjustmentDict is the dictionary that sets the information for correct data loading and trajectory and data synchronization
adjustmentDict = Dict{Symbol,Any}()
adjustmentDict[:reconSize] = Tuple(params_general[:reconSize])
adjustmentDict[:interleave] = startIndexIntlv
adjustmentDict[:numSamples] = params_general[:numADCSamples]
adjustmentDict[:delay] = 0.00000 # naive delay correction

adjustmentDict[:interleaveDataFileNames] = params_general[:fullPathScan]

adjustmentDict[:trajFilename] = params_general[:fullPathGradient]
adjustmentDict[:excitations] = sliceSelection

adjustmentDict[:doMultiInterleave] = !isDataSingleIntlv
adjustmentDict[:doOddInterleave] = false
adjustmentDict[:numInterleaves] = isDataSingleIntlv ? 1 : length(adjustmentDict[:interleaveDataFileNames]) # one interleaf per file, count files, if filenames are array of strings (not only one string)

adjustmentDict[:singleSlice] = !isMultiSlice

# Defined recon size and parameters for data loading
@info "Using Parameters:\n\nreconSize = $(adjustmentDict[:reconSize]) \n interleave = $(adjustmentDict[:interleave]) \n slices = $(sliceChoice) \n coils = $(size(cartesian_sensitivity, 4)) \n numSamples = $(adjustmentDict[:numSamples])\n\n"

if reloadGIRFData || !(@isdefined gK1) || !(@isdefined gAK1) || !(@isdefined gK0) || !(@isdefined gAK0)
    @info "Loading Gradient Impulse Response Functions"

    ## Load GIRFs (K1)
    gK1 = readGIRFFile(params_general[:fullPathGIRF][1], params_general[:fullPathGIRF][2], params_general[:fullPathGIRF][3], "GIRF_FT", false)
    gAk1 = GirfApplier(gK1, gamma)

    ## Load K₀ GIRF
    gK0 = readGIRFFile(params_general[:fullPathGIRF][1], params_general[:fullPathGIRF][2], params_general[:fullPathGIRF][3], "b0ec_FT", true)
    gAk0 = GirfApplier(gK0, gamma)
end

## Only load data when it has not been done yet, or it's specifically required.
if reloadSpiralData || !(@isdefined acqDataImaging)
    ## Convert raw to AcquisitionData

    @info "Reading spiral data and merging interleaves"
    acqDataImaging = merge_raw_interleaves(adjustmentDict)

    if params_general[:doCorrectWithGIRFkxyz]
        @info "Correcting For GIRF"
        apply_girf!(acqDataImaging, gAk1)
    end

    if params_general[:doCorrectWithGIRFk0]
        @info "Correcting For k₀"
        apply_k0!(acqDataImaging, gAk0)
    end

    ## Check the k-space nodes so they don't exceed frequency limits [-0.5, 0.5] (inclusive)
    check_acquisition_nodes!(acqDataImaging)

end

## Sense Map loading
@info "Resizing Sense Maps"

# Resize sense maps to match encoding size of data matrix
sensitivity = mapslices(x -> imresize(x, adjustmentDict[:reconSize][1], adjustmentDict[:reconSize][2]), cartesian_sensitivity, dims = [1, 2])

# Plot the sensitivity maps of each coil
@info "Plotting SENSE Maps"

if params_general[:do_plot_recon]
    plot_sense_maps(sensitivity, size(sensitivity, 4), sliceIndex = 10)
end


# shift FOV to middle :) 
shiftksp!(acqDataImaging, params_general[:fovShift])
# changeFOV!(acqDataImaging,[0.99, 0.99])


## Do coil compression to make recon faster
if params_general[:doCoilCompression]
    acqDataImaging, sensitivity = geometricCC_2d(acqDataImaging, sensitivity, params_general[:nVirtualCoils])
end


## B0 Maps (Assumes we have a B0 map from gradient echo scan named b0)
@info "Resizing B0 Maps"
resizedB0 = mapslices(x -> imresize(x, adjustmentDict[:reconSize][1], adjustmentDict[:reconSize][2]), b0Maps, dims = [1, 2])

## Define Parameter Dictionary for use with reconstruction
# CAST TO ComplexF32 if you're using current MRIReco.jl

@info "Setting Parameters"
params = Dict{Symbol,Any}()
params[:reco] = "multiCoil"
params[:reconSize] = adjustmentDict[:reconSize][1:2]
params[:regularization] = "L2"
params[:λ] = 1e-2 # CHANGE THIS TO GET BETTER OR WORSE RECONSTRUCTION RESULTS
params[:iterations] = params_general[:nReconIterations]
params[:solver] = "cgnr"
params[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false)
params[:senseMaps] = ComplexF32.(sensitivity[:, :, selectedSlice, :])

if params_general[:do_correct_with_b0_map]
    params[:correctionMap] = ComplexF32.(-1im .* resizedB0[:, :, selectedSlice])
end

## Call to reconstruction
@info "Performing Reconstruction"
@time reco = reconstruction(acqDataImaging, params)
#totalRecon = sum(abs2,reco.data,dims=5)

# save Map recon (multi-echo etc.)
if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    resolution_tmp = fieldOfView(acqDataImaging)[1:2] ./ encodingSize(acqDataImaging)
    resolution_mm = (resolution_tmp[1], resolution_tmp[2], fieldOfView(acqDataImaging)[3] * (1 + params_general[:sliceDistanceFactor_percent] / 100.0)) # for 2D only, since FOV[3] is slice thickness then, but gap has to be observed

    # TODO: use slice ordering from cartesian scan directly!
    nSlices = numSlices(acqDataImaging)
    sliceIndexArray = get_slice_order(nSlices, isSliceInterleaved = true)
    save_map(
        params_general[:fullPathSaveRecon],
        params_general[:scalingFactorSaveRecon] * reco.data[:, :, sliceIndexArray],
        resolution_mm;
        doSplitPhase = true,
        doNormalize = params_general[:doNormalizeRecon],
    )
end

if params_general[:do_plot_recon]
    @info "Plotting Reconstruction"
    #pygui(true)
    plot_reconstruction(
        reco,
        1:length(selectedSlice),
        resizedB0[:, :, selectedSlice],
        figHandles = ["Original Magnitude", "Original Phase", "B0"],
        isSliceInterleaved = true,
        rotateAngle = 270,
    )
end

@info "Successfully Completed SpiralRecon"
