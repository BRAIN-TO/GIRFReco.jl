using HDF5, MRIReco, LinearAlgebra, DSP, FourierTools, ROMEO, MRIGradients

include("../utils/Utils.jl")
include("../utils/fieldMapEstimator.jl")
include("../utils/shiftksp.jl")

## function to calculate the B0 maps from the two images with different echo times
# TODO have the b0 map calculation be capable of handling variable echo times
function calculateB0Maps(imData,slices,echoTime1,echoTime2)

    # b0Maps = mapslices(x -> rotl90(x),ROMEO.unwrap(angle.(imData[:,:,slices,2,1].*conj(imData[:,:,slices,1,1]))),dims=(1,2))./((7.38-4.92)/1000)
    b0Maps = mapslices(x -> x, ROMEO.unwrap(angle.(imData[:,:,slices,2,1].*conj(imData[:,:,slices,1,1]))),dims=(1,2))./((echoTime2-echoTime1)/1000)

end

## Load data files
makeMaps = true
saveMaps = true

# Echo times for field map raw data, in ms
TE1 = 4.92
TE2 = 7.38 

reconSize = (200,200)

@info "Loading Data Files"

b0FileName = "data/Fieldmaps/field_map_132_2.h5"
processedFileName = "data/Fieldmaps/processedCartesianData.h5" # filename for preprocessed data 

if makeMaps

    # Set the data file name (Change this for your own system)
    dataFileCartesian = ISMRMRDFile(b0FileName)

    # read in the raw data from the ISMRMRD file into a RawAcquisitionData object
    r = RawAcquisitionData(dataFileCartesian)

    # Preprocess Data and save!
    preprocessCartesianData(r::RawAcquisitionData, saveMaps; fname = processedFileName)

end

# removeOversampling!(r)

# Load preprocessed data!
dataFileNew = ISMRMRDFile(processedFileName)
rawDataNew = RawAcquisitionData(dataFileNew)
acqDataCartesian= AcquisitionData(rawDataNew, estimateProfileCenter=true)

# Define coils and slices
nCoils = size(acqDataCartesian.kdata[1],2)
nSlices = numSlices(acqDataCartesian)

# shift FOV to middle :) 
shiftksp!(acqDataCartesian,[0,-20])
#changeFOV!(acqDataCartesian,[1.5,1.5])

## Don't have to recalculate sense maps for both scans but possibly it could make a
#  difference in Diffusion scans

@info "Calculating Sense Maps"
senseCartesian = espirit(acqDataCartesian,(4,4),12,eigThresh_1=0.01, eigThresh_2=0.98)
sensitivity = senseCartesian

plotSenseMaps(sensitivity,nCoils)

acqDataCartesian.traj[1].cartesian = false
acqDataCartesian.traj[2].cartesian = false
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
## Defining array mapping from acquisition number to slice number (indexArray[slice = 1:9] = [acquisitionNumbers])

# indexArray = [5,1,6,2,7,3,8,4,9] # for 9 slice phantom
indexArray = [8,1,9,2,10,3,11,4,12,5,13,6,14,7,15] # for 15 slice phantom
#indexArray = 1 # for 1 slice phantom

## Call the reconstruction function

@info "Performing Reconstruction"
@time cartesianReco = reconstruction(acqDataCartesian,paramsCartesian)

## Calculate B0 maps from the acquired images (if two TEs)

slices = 1:length(indexArray)

@info "Calculating B0 Maps"
b0Maps = calculateB0Maps(cartesianReco.data,slices, TE1, TE2)
@time b0Maps2 = estimateB0Maps(cartesianReco.data,slices,TE1,TE2,true; β = 0.5, reltol = 1e-4)