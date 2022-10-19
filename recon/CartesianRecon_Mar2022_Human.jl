using HDF5, MRIReco, LinearAlgebra, DSP, FourierTools, ROMEO, MRIGradients, NIfTI

include("../utils/Utils.jl")
include("../utils/fieldMapEstimator.jl")

## function to calculate the B0 maps from the two images with different echo times
# TODO have the b0 map calculation be capable of handling variable echo times
function calculateB0Maps(imData,slices,echoTime1,echoTime2)

    # b0Maps = mapslices(x -> rotl90(x),ROMEO.unwrap(angle.(imData[:,:,slices,2,1].*conj(imData[:,:,slices,1,1]))),dims=(1,2))./((7.38-4.92)/1000)
    b0Maps = mapslices(x -> x, ROMEO.unwrap(angle.(imData[:,:,slices,2,1].*conj(imData[:,:,slices,1,1]))),dims=(1,2))./((echoTime2-echoTime1)/1000)

end

## Dictionary of frequently changed parameters
include("ReconConfig.jl")

## Load data files

makeMaps = true
saveMaps = true

# Echo times for field map raw data, in ms
TE1 = 4.92
TE2 = 7.38 

@info "Loading Data Files"

b0FileName = paramsGeneral[:fullPathMultiEcho];

# filename for preprocessed data (remove oversampling, permute dimensions wrt MRIReco)
processedFileName = paramsGeneral[:fullPathProcessedCartesian] 

if makeMaps

    # Set the data file name (Change this for your own system)
    dataFileCartesian = ISMRMRDFile(b0FileName)

    # read in the raw data from the ISMRMRD file into a RawAcquisitionData object
    r = RawAcquisitionData(dataFileCartesian)

    # Preprocess Data and save!
    preprocessCartesianData(r::RawAcquisitionData, saveMaps; fname = processedFileName)

end

# Load preprocessed data!
dataFileNew = ISMRMRDFile(processedFileName)
rawDataNew = RawAcquisitionData(dataFileNew)
acqDataCartesian= AcquisitionData(rawDataNew, estimateProfileCenter=true)

# Define coils and slices
nCoils = size(acqDataCartesian.kdata[1],2)
nSlices = numSlices(acqDataCartesian)


sliceIndexArray = getSliceOrder(nSlices, isSliceInterleaved = true)



## Don't have to recalculate sense maps for both scans but possibly it could make a
#  difference in Diffusion scans

@info "Calculating Sense Maps"
senseCartesian = espirit(acqDataCartesian,(6,6),24,eigThresh_1=0.02, eigThresh_2=0.95)
sensitivity = senseCartesian


# save SENSE maps
if paramsGeneral[:doSaveRecon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    resolution_mm[1:2] = fieldOfView(acqDataCartesian)[1:2]./size(sensitivity)[1:2]
    resolution_mm[3] = fieldOfView(acqDataCartesian)[3]; # for 2D only, since FOV[3] is slice thickness then

    # TODO: use correct slice order everywhere, e.g., when saving/loading maps for spiral recon
    saveMap(paramsGeneral[:fullPathSaveSense], sensitivity[:,:,sliceIndexArray,:], resolution_mm; doSplitPhase=true)
end

## Parameter dictionary definition for reconstruction

@info "Setting Parameters"
paramsCartesian = Dict{Symbol,Any}() # instantiate dictionary
paramsCartesian[:reco] = "multiCoil" # choose multicoil reconstruction
paramsCartesian[:reconSize] = (acqDataCartesian.encodingSize[1],acqDataCartesian.encodingSize[2]) # set recon size to be the same as encoded size
paramsCartesian[:regularization] = "L2" # choose regularization for the recon algorithm
paramsCartesian[:λ] = 1.e-2 # recon parameter (there may be more of these, need to dig into code to identify them for solvers other than cgnr)
paramsCartesian[:iterations] = 20 # number of CG iterations
paramsCartesian[:solver] = "cgnr" # inverse problem solver method
paramsCartesian[:solverInfo] = SolverInfo(ComplexF32,store_solutions=false) # turn on store solutions if you want to see the reconstruction convergence (uses more memory)
paramsCartesian[:senseMaps] = ComplexF32.(sensitivity) # set sensitivity map array
# paramsCartesian[:correctionMap] = ComplexF32.(-1im.*b0Maps)




## Call the reconstruction function

@info "Performing Reconstruction"
@time cartesianReco = reconstruction(acqDataCartesian,paramsCartesian)

## Calculate B0 maps from the acquired images (if two TEs)

slices = 1:length(indexArray)

@info "Calculating B0 Maps"
# b0Maps = calculateB0Maps(cartesianReco.data,slices, TE1, TE2)
b0Maps2 = estimateB0Maps(cartesianReco.data,slices,TE1,TE2,0.00001,true)

# save B0 map
if paramsGeneral[:doSaveRecon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    resolution_mm[1:2] = fieldOfView(acqDataCartesian)[1:2]./size(b0Maps2)[1:2]
    resolution_mm[3] = fieldOfView(acqDataCartesian)[3]; # for 2D only, since FOV[3] is slice thickness then
    saveMap(paramsGeneral[:fullPathSaveB0], b0Maps2[:,:,sliceIndexArray], resolution_mm)
end


@info "Plotting Cartesian Results (Sensitivity Maps and B0 Maps)"
pygui(true) # Leave this code till we need plotting.
# plotSenseMaps(sensitivity,nCoils)
plotReconstruction(cartesianReco[:,:,:,1], 1:size(cartesianReco,3), b0Maps2, isSliceInterleaved = true, rotateAngle = 270)

@info "Successfully Completed CartesianReconstruction"
