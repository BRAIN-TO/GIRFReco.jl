using Flux, HDF5, MRIReco, LinearAlgebra, DSP, FourierTools, ROMEO, MRIGradients

include("../utils/Utils.jl")
include("../utils/fieldMapEstimator.jl")

## function to calculate the B0 maps from the two images with different echo times
# TODO have the b0 map calculation be capable of handling variable echo times
function calculateB0Maps(imData,slices)

    b0Maps = mapslices(x -> rotl90(x),ROMEO.unwrap(angle.(imData[:,:,slices,2,1].*conj(imData[:,:,slices,1,1]))),dims=(1,2))./((7.38-4.92)/1000)

end

## Load data files

makeMaps = true
saveMaps = true

reconSize = (64,64)

@info "Loading Data Files"

if makeMaps

    # Set the data file name (Change this for your own system)
    dataFileCartesian = ISMRMRDFile("data/Fieldmaps/fieldMap_105_gre_2.h5")

    # read in the raw data from the ISMRMRD file into a RawAcquisitionData object
    r = RawAcquisitionData(dataFileCartesian)

    # Set filename for preprocessed data 
    fname = "data/testFile.h5"

    # Preprocess Data and save!
    preprocessCartesianData(r::RawAcquisitionData, saveMaps; fname)

end

# removeOversampling!(r)

# Load preprocessed data!
dataFileNew = ISMRMRDFile("data/testFile.h5")
rawDataNew = RawAcquisitionData(dataFileNew)
acqDataCartesian= AcquisitionData(rawDataNew, estimateProfileCenter=true)

# Define coils and slices
nCoils = size(acqDataCartesian.kdata[1],2)
nSlices = numSlices(acqDataCartesian)

## Don't have to recalculate sense maps for both scans but possibly it could make a
#  difference in Diffusion scans

@info "Calculating Sense Maps"
senseCartesian = espirit(acqDataCartesian,(6,6),24,eigThresh_1=0.02, eigThresh_2=0.95)
sensitivity = senseCartesian

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

indexArray = [5,1,6,2,7,3,8,4,9] # for 9 slice phantom
#indexArray = 1 # for 1 slice phantom

## Call the reconstruction function

@info "Performing Reconstruction"
cartesianReco = reconstruction(acqDataCartesian,paramsCartesian)

## Calculate B0 maps from the acquired images (if two TEs)

slices = 1:length(indexArray)

@info "Calculating B0 Maps"
b0Maps = calculateB0Maps(cartesianReco.data,slices)
b0Maps2 = estimateB0Maps(cartesianReco.data,slices,4.92,7.38,0.00001,false)

@info "Successfully Completed CartesianReconstruction"
