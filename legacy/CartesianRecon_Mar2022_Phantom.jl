using HDF5, MRIReco, LinearAlgebra, DSP, FourierTools, ROMEO, MRIGradients

include("../utils/Utils.jl")

## function to calculate the B0 maps from the two images with different echo times
# TODO have the b0 map calculation be capable of handling variable echo times
function calculate_b0_maps(imData, slices, echoTime1, echoTime2)

    # b0Maps = mapslices(x -> rotl90(x),ROMEO.unwrap(angle.(imData[:,:,slices,2,1].*conj(imData[:,:,slices,1,1]))),dims=(1,2))./((7.38-4.92)/1000)
    b0Maps =
        mapslices(x -> x, ROMEO.unwrap(angle.(imData[:, :, slices, 2, 1] .* conj(imData[:, :, slices, 1, 1]))), dims = (1, 2)) ./
        ((echoTime2 - echoTime1) / 1000)

end

## Load data files

makeMaps = true
save_maps = true

# Echo times for field map raw data, in ms
TE1 = 4.92
TE2 = 7.38

@info "Loading Data Files"

b0FileName = "D:\\OneDrive - UHN\\MRP-SPIDI\\SPIDI\\data\\SPIDI_0007\\Phantom\\dat\\field_map_83_2.h5"
processedFileName = "D:\\OneDrive - UHN\\MRP-SPIDI\\SPIDI\\data\\SPIDI_0007\\Phantom\\dat\\processedCartesianData.h5" # filename for preprocessed data 

if makeMaps

    # Set the data file name (Change this for your own system)
    dataFileCartesian = ISMRMRDFile(b0FileName)

    # read in the raw data from the ISMRMRD file into a RawAcquisitionData object
    r = RawAcquisitionData(dataFileCartesian)

    # Preprocess Data and save!
    preprocess_cartesian_data(r::RawAcquisitionData, save_maps; fname = processedFileName)

end

# remove_oversampling!(r)

# Load preprocessed data!
dataFileNew = ISMRMRDFile(processedFileName)
rawDataNew = RawAcquisitionData(dataFileNew)
acqDataCartesian = AcquisitionData(rawDataNew, estimateProfileCenter = true)

# Define coils and slices
nCoils = size(acqDataCartesian.kdata[1], 2)
nSlices = numSlices(acqDataCartesian)

## Don't have to recalculate sense maps for both scans but possibly it could make a
#  difference in Diffusion scans

@info "Calculating Sense Maps"
cartesian_sensitivity = espirit(acqDataCartesian, (4, 4), 10, eigThresh_1 = 0.01, eigThresh_2 = 0.98)
sensitivity = cartesian_sensitivity

plot_sense_maps(sensitivity, nCoils)

## Parameter dictionary definition for reconstruction

@info "Setting Parameters"
params_cartesian = Dict{Symbol,Any}() # instantiate dictionary
params_cartesian[:reco] = "multiCoil" # choose multicoil reconstruction
params_cartesian[:reconSize] = (acqDataCartesian.encodingSize[1], acqDataCartesian.encodingSize[2]) # set recon size to be the same as encoded size
params_cartesian[:regularization] = "L2" # choose regularization for the recon algorithm
params_cartesian[:λ] = 1.e-2 # recon parameter (there may be more of these, need to dig into code to identify them for solvers other than cgnr)
params_cartesian[:iterations] = 20 # number of CG iterations
params_cartesian[:solver] = "cgnr" # inverse problem solver method
params_cartesian[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false) # turn on store solutions if you want to see the reconstruction convergence (uses more memory)
params_cartesian[:senseMaps] = ComplexF32.(sensitivity) # set sensitivity map array
# params_cartesian[:correctionMap] = ComplexF32.(-1im.*b0Maps)
## Defining array mapping from acquisition number to slice number (indexArray[slice = 1:9] = [acquisitionNumbers])

# indexArray = [5,1,6,2,7,3,8,4,9] # for 9 slice phantom
indexArray = [8, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15] # for 15 slice phantom
#indexArray = 1 # for 1 slice phantom

## Call the reconstruction function

@info "Performing Reconstruction"
cartesianReco = reconstruction(acqDataCartesian, params_cartesian)

## Calculate B0 maps from the acquired images (if two TEs)

slices = 1:length(indexArray)

@info "Calculating B0 Maps"
# b0Maps = calculate_b0_maps(cartesianReco.data,slices, TE1, TE2)
b0Maps = estimateB0Maps(cartesianReco.data, slices, TE1, TE2, 0.00001, true)

@info "Plotting Cartesian Results (Sensitivity Maps and B0 Maps) \n"
pygui(true) # Leave this code till we need plotting.
# plot_sense_maps(sensitivity,nCoils)
plot_reconstruction(cartesianReco[:, :, :, 1], 1:size(cartesianReco, 3), b0Maps, isSliceInterleaved = true, rotateAngle = 270)

@info "Successfully Completed CartesianReconstruction"
