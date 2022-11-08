using HDF5, MRIReco, LinearAlgebra, DSP, FourierTools, ROMEO, MRIGradients

include("../utils/Utils.jl")
include("../utils/fieldMapEstimator.jl")

## Dictionary of frequently changed parameters
include("ReconConfig.jl")

## Load data files

# Echo times for field map raw data, in ms
TE1 = 4.92
TE2 = 7.38 

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
acqDataCartesian= AcquisitionData(rawDataNew, estimateProfileCenter=true)

# Define coils and slices
nCoils = size(acqDataCartesian.kdata[1],2)
nSlices = numSlices(acqDataCartesian)


sliceIndexArray = getSliceOrder(nSlices, isSliceInterleaved = true)
# shift FOV to middle :) 
#TODO: in MRIReco v0.7, try: correctOffset(acqDataCartesian, [0 -20 0])
shiftksp!(acqDataCartesian, paramsGeneral[:fovShift])
#changeFOV!(acqDataCartesian,[1.5,1.5])

## Don't have to recalculate sense maps for both scans but possibly it could make a
#  difference in Diffusion scans

@info "Calculating Sense Maps"

# from reading docstring, min thresholding
senseCartesian = espirit(acqDataCartesian,(6,6),30,eigThresh_1=0.01, eigThresh_2=0.9)

# from Lars (7T spirals)
#senseCartesian = espirit(acqDataCartesian,(6,6),30,eigThresh_1=0.02, eigThresh_2=0.98)
# from Alexander (Phantom?)
#senseCartesian = espirit(acqDataCartesian,(4,4),12,eigThresh_1=0.01, eigThresh_2=0.98)

# normalize for consistency with saving/loading and better ranges of reconstruction values
senseCartesian /= maximum(abs.(senseCartesian))
sensitivity = senseCartesian 

resolution_mm = fieldOfView(acqDataCartesian)./size(sensitivity)[1:3]
resolution_mm[3] = fieldOfView(acqDataCartesian)[3] *(1 + paramsGeneral[:sliceDistanceFactor_percent]/100.0); # for 2D only, since FOV[3] is slice thickness then, but gap has to be observed


# save SENSE maps
if paramsGeneral[:doSaveRecon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    # TODO: use correct slice order everywhere, e.g., when saving/loading maps for spiral recon
    saveMap(paramsGeneral[:fullPathSaveSense], sensitivity[:,:,sliceIndexArray,:], resolution_mm; doSplitPhase=true)
end

plotSenseMaps(sensitivity,nCoils)

#acqDataCartesian.traj[1].cartesian = false
#acqDataCartesian.traj[2].cartesian = false
## Parameter dictionary definition for reconstruction

@info "Setting Parameters"
paramsCartesian = Dict{Symbol,Any}() # instantiate dictionary
paramsCartesian[:reco] = "multiCoil" # choose multicoil reconstruction

# TODO: make recon size and FOV variable!
paramsCartesian[:reconSize] = (acqDataCartesian.encodingSize[1],acqDataCartesian.encodingSize[2]) # set recon size to be the same as encoded size
paramsCartesian[:regularization] = "L2" # choose regularization for the recon algorithm
paramsCartesian[:λ] = 1.e-2 # recon parameter (there may be more of these, need to dig into code to identify them for solvers other than cgnr)
paramsCartesian[:iterations] = 20 # number of CG iterations
paramsCartesian[:solver] = "cgnr" # inverse problem solver method
paramsCartesian[:solverInfo] = SolverInfo(ComplexF32,store_solutions=false) # turn on store solutions if you want to see the reconstruction convergence (uses more memory)
paramsCartesian[:senseMaps] = ComplexF32.(sensitivity) # set sensitivity map array


## Call the reconstruction function

@info "Performing Reconstruction"
@time cartesianReco = reconstruction(acqDataCartesian,paramsCartesian)

# save Map recon (multi-echo etc.)
if paramsGeneral[:doSaveRecon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    saveMap(paramsGeneral[:fullPathSaveMapRecon], cartesianReco.data[:,:,sliceIndexArray,:,:,:], resolution_mm; doSplitPhase=true)
end

## Calculate B0 maps from the acquired images (if two TEs)

slices = 1:length(sliceIndexArray)

@info "Calculating B0 Maps"
# b0Maps = calculateB0Maps(cartesianReco.data,slices, TE1, TE2)
b0Maps = estimateB0Maps(cartesianReco.data,slices,TE1,TE2,true; β = 0.01, reltol = 1e-4)

# save B0 map
if paramsGeneral[:doSaveRecon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    saveMap(paramsGeneral[:fullPathSaveB0], b0Maps[:,:,sliceIndexArray], resolution_mm; doNormalize = false) # no normalization, we want absolute values for offres maps
end


if paramsGeneral[:doPlotRecon]
    @info "Plotting Cartesian Results (Sensitivity Maps and B0 Maps)"
    pygui(true) # Leave this code till we need plotting.
    # plotSenseMaps(sensitivity,nCoils)
    plotReconstruction(cartesianReco[:,:,:,1], 1:size(cartesianReco,3), b0Maps, isSliceInterleaved = true, rotateAngle = 270)
end

# cleanup unused file
if !paramsGeneral[:doSaveProcessedMapScan] && paramsGeneral[:doProcessMapScan]
    rm(processedFileName)
end

@info "Successfully Completed CartesianReconstruction"

