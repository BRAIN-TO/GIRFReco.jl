using PyPlot, HDF5, MRIReco, LinearAlgebra, DSP, FourierTools, ROMEO

include("../utils/utils.jl")

## function to calculate the B0 maps from the two images with different echo times
# TODO have the b0 map calculation be capable of handling variable echo times
function calculateB0Maps(imData,slices)

    b0Maps = mapslices(x -> rotl90(x),ROMEO.unwrap(angle.(imData[:,:,slices,2,1].*conj(imData[:,:,slices,1,1]))),dims=(1,2))./((7.38-4.92)/1000)

end

## function to crop the B0 map FOV down to the same aspect ratio as that of the spiral scans (our data is 2x oversampled in the readout direction so the default FOV is 2x the size it needs to be in that direction)
function cropB0Maps(b0_l)


   b0 = mapslices(x->x[:,33:96],b0_l,dims=(1,2))

end

## Load data files

reconSize = (400,200)

@info "Loading Data Files"
# Set the data file name (Change this for your own system)
dataFileCartesian = ISMRMRDFile("data/Fieldmaps/fieldMap_30_2.h5")

# read in the raw data from the ISMRMRD file into a RawAcquisitionData object
r = RawAcquisitionData(dataFileCartesian)

# Set filename for preprocessed data 
fname = "data/testFile.h5"

# Preprocess Data and save!
preprocessCartesianData(r::RawAcquisitionData, fname)

# Load preprocessed data!
dataFileNew = ISMRMRDFile("/home/ajaffray/Documents/testFile.h5")
rawDataNew = RawAcquisitionData(dataFileNew)
acqDataCartesian= AcquisitionData(rawDataNew, estimateProfileCenter=true)

# Define coils and slices
nCoils = size(acqDataCartesian.kdata[1],2)
nSlices = numSlices(acqDataCartesian)

## Don't have to recalculate sense maps for both scans but possibly it could make a
#  difference in Diffusion scans

@info "Calculating Sense Maps"
@time senseCartesian = espirit(acqDataCartesian,(6,6),30,eigThresh_1=0.05, eigThresh_2=0.98)
sensitivity = senseCartesian

## Resize sense maps to match encoding size of data matrix
sensitivity = imresize(senseCartesian,(acqDataCartesian.encodingSize[1],acqDataCartesian.encodingSize[2],nSlices,nCoils))
#sensitivity = mapslices(rotl90,sensitivity,dims=[1,2])

## Visualization

@info "Plotting Sense Maps"
plotSenseMaps(sensitivity,nCoils)

## Parameter dictionary definition for reconstruction

@info "Setting Parameters"
paramsCartesian = Dict{Symbol,Any}() # instantiate dictionary
paramsCartesian[:reco] = "multiCoil" # choose multicoil reconstruction
paramsCartesian[:reconSize] = (acqDataCartesian.encodingSize[1],acqDataCartesian.encodingSize[2]) # set recon size to be the same as encoded size
paramsCartesian[:regularization] = "L2" # choose regularization for the recon algorithm
paramsCartesian[:λ] = 1.e-2 # recon parameter (there may be more of these, need to dig into code to identify them for solvers other than cgnr)
paramsCartesian[:iterations] = 20 # number of CG iterations
paramsCartesian[:solver] = "cgnr" # inverse problem solver method
paramsCartesian[:solverInfo] = SolverInfo(ComplexF64,store_solutions=false) # turn on store solutions if you want to see the reconstruction convergence (uses more memory)
paramsCartesian[:senseMaps] = sensitivity # set sensitivity map array

## Defining array mapping from acquisition number to slice number (indexArray[slice = 1:9] = [acquisitionNumbers])

#indexArray = [5,1,6,2,7,3,8,4,9] # for 9 slice phantom
indexArray = 1 # for 1 slice phantom

## Call the reconstruction function

@info "Performing Reconstruction"
cartesianReco = reconstruction(acqDataCartesian,paramsCartesian)
#@time reco2 = reconstruction(acqData2,params)

## If not multi-coil, perform L2 combination of coil images
# if multi-coil, this doesn't do anything to the data

reconstructed = sum(abs2,cartesianReco.data,dims=5)

## Calculate B0 maps from the acquired images (if two TEs)

slices = 1:length(indexArray)
b0Maps = calculateB0Maps(cartesianReco.data,slices)
b0 = cropB0Maps(b0Maps)

## Resize fourier data and replot (Only for comparison purposes with hi-res spiral images)

#resampledRecon = Array{ComplexF64,5}(undef,(reconSize[1], reconSize[2],1,2,1))

resampledRecon = mapslices(x->rotl90(FourierTools.resample(x,reconSize))[:,101:300], cartesianReco.data; dims=[1,2])
resampledB0 = mapslices(x->FourierTools.resample(x,(200,200)), b0;dims=[1,2])

## Plotting Reconstruction

@info "Plotting Reconstruction"
plotReconstruction(resampledRecon,indexArray,b0)
