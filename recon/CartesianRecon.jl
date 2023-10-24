using HDF5, MRIReco, LinearAlgebra, DSP, FourierTools, ROMEO, MRIGradients, MRIFiles, MRIFieldmaps

include("../utils/fieldMapEstimator.jl")

## Dictionary of frequently changed parameters
# include("ReconConfig.jl")

## Load data files

@info "Loading Data Files"

b0FileName = paramsGeneral[:fullPathMapScan];

# filename for preprocessed data (remove oversampling, permute dimensions wrt MRIReco)
processedFileName = paramsGeneral[:fullPathProcessedMapScan]

if paramsGeneral[:doProcessMapScan]

    # Set the data file name (Change this for your own system)
    dataFileCartesian = ISMRMRDFile(b0FileName)

    # read in the raw data from the ISMRMRD file into a RawAcquisitionData object
    r = RawAcquisitionData(dataFileCartesian)

    # does not change anything...
    # r.params["reconFOV"] = [230, 230, 2]

    # Preprocess Data and save!
    preprocessCartesianData(r::RawAcquisitionData, true; fname = processedFileName)

end

# Load preprocessed data!
dataFileNew = ISMRMRDFile(processedFileName)


rawDataNew = RawAcquisitionData(dataFileNew)
acqDataCartesian = AcquisitionData(rawDataNew, estimateProfileCenter = true)

# Define coils and slices
nCoils = size(acqDataCartesian.kdata[1], 2)
nSlices = numSlices(acqDataCartesian)

# Define TEs 
# Echo times for field map raw data, in ms
TE1 = rawDataNew.params["TE"][1]
TE2 = rawDataNew.params["TE"][2]

sliceIndexArray = getSliceOrder(nSlices, isSliceInterleaved = true)
# shift FOV to middle :) 
#TODO: in MRIReco v0.7, try: correctOffset(acqDataCartesian, [0 -20 0])
shiftksp!(acqDataCartesian, paramsGeneral[:fovShift]) # amount of FOV shift; in unit of number of voxels in [x,y] direction

#changeFOV!(acqDataCartesian,[1.5,1.5])

## Don't have to recalculate sense maps for both scans but possibly it could make a
#  difference in Diffusion scans

@info "Calculating Sense Maps"

# from reading docstring, min thresholding
senseCartesian = espirit(acqDataCartesian, (6, 6), 30, eigThresh_1 = 0.01, eigThresh_2 = 0.9)

# from Lars (7T spirals)
#senseCartesian = espirit(acqDataCartesian,(6,6),30,eigThresh_1=0.02, eigThresh_2=0.98)
# from Alexander (Phantom?)
#senseCartesian = espirit(acqDataCartesian,(4,4),12,eigThresh_1=0.01, eigThresh_2=0.98)

# normalize for consistency with saving/loading and better ranges of reconstruction values
senseCartesian /= maximum(abs.(senseCartesian))
sensitivity = senseCartesian

res_x = fieldOfView(acqDataCartesian)[1] ./ size(sensitivity)[1]
res_y = fieldOfView(acqDataCartesian)[2] ./ size(sensitivity)[2]
res_z = fieldOfView(acqDataCartesian)[3] .* (1 + paramsGeneral[:sliceDistanceFactor_percent] ./ 100.0) # for 2D only, since FOV[3] is slice thickness then, but gap has to be observed
resolution_mm = (res_x, res_y, res_z)

# save SENSE maps
if paramsGeneral[:doSaveRecon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    # TODO: use correct slice order everywhere, e.g., when saving/loading maps for spiral recon
    saveMap(paramsGeneral[:fullPathSaveSense], sensitivity[:, :, sliceIndexArray, :], resolution_mm; doSplitPhase = true)
end

## Parameter dictionary definition for reconstruction

@info "Setting Parameters"
paramsCartesian = Dict{Symbol,Any}() # instantiate dictionary
paramsCartesian[:reco] = "multiCoil" # choose multicoil reconstruction
paramsCartesian[:reconSize] = (acqDataCartesian.encodingSize[1], acqDataCartesian.encodingSize[2]) # set recon size to be the same as encoded size
paramsCartesian[:regularization] = "L2" # choose regularization for the recon algorithm
paramsCartesian[:λ] = 1.e-2 # recon parameter (there may be more of these, need to dig into code to identify them for solvers other than cgnr)
paramsCartesian[:iterations] = paramsGeneral[:nReconIterations] # number of CG iterations
paramsCartesian[:solver] = "cgnr" # inverse problem solver method
paramsCartesian[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false) # turn on store solutions if you want to see the reconstruction convergence (uses more memory)
paramsCartesian[:senseMaps] = ComplexF32.(sensitivity) # set sensitivity map array


## Call the reconstruction function

@info "Performing Reconstruction"
@time cartesianReco = reconstruction(acqDataCartesian, paramsCartesian)

# save Map recon (multi-echo etc.)
if paramsGeneral[:doSaveRecon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    saveMap(paramsGeneral[:fullPathSaveMapRecon], cartesianReco.data[:, :, sliceIndexArray, :, :, :], resolution_mm; doSplitPhase = true)
end

## Calculate B0 maps from the acquired images (if two TEs)
@info "Calculating B0 Maps"
slices = 1:length(sliceIndexArray)
b0Maps = zeros(200, 200, 15)
b0Method = "2D_2008" # Can be "Simple","2D_2008" or "3D_2020" (How do we incorporate this into the recon demo?)

if b0Method == "2D_2008"

    b0Maps = estimateB0Maps(cartesianReco.data, slices, TE1, TE2, true; β = paramsGeneral[:b0mapSmoothBeta], reltol = 1e-4)

elseif b0Method == "3D_2020"

    niter = 100 # usually good enough

    for i = 1:15

        (b0Maps[:, :, i], times, out) = b0map(
            reshape(cartesianReco.data[:, :, i, :, 1, 1], (200, 200, 1, 2)) ./ maximum(abs.(cartesianReco.data)),
            [TE1 / 1000, TE2 / 1000];
            order = 2,
            l2b = -6,
            gamma_type = :PR,
            niter = niter,
            precon = :diag,
            track = false,
        )

    end

    b0Maps = 2 * pi * b0Maps

end

# save B0 map
if paramsGeneral[:doSaveRecon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    saveMap(paramsGeneral[:fullPathSaveB0], b0Maps[:, :, sliceIndexArray], resolution_mm; doNormalize = false) # no normalization, we want absolute values for offres maps
end


if paramsGeneral[:doPlotRecon]
    @info "Plotting Cartesian Results (Sensitivity Maps and B0 Maps)"
    # plotSenseMaps(sensitivity,nCoils)
    plotlyjs()
    plotReconstruction(cartesianReco[:, :, :, 1], 1:size(cartesianReco, 3), b0Maps, isSliceInterleaved = false, rotateAngle = 180)
end

# cleanup unused file
if !paramsGeneral[:doSaveProcessedMapScan] && paramsGeneral[:doProcessMapScan]
    rm(processedFileName)
end

@info "Successfully Completed CartesianReconstruction"

