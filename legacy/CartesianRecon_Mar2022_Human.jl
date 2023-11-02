using HDF5, MRIReco, LinearAlgebra, DSP, FourierTools, ROMEO, MRIGradients

include("../utils/Utils.jl")
include("../utils/FieldMapEstimator.jl")

## Dictionary of frequently changed parameters
include("ReconConfig.jl")

## Load data files

# Echo times for field map raw data, in ms
TE1 = params_general[:mapTEs_ms][1]
TE2 = params_general[:mapTEs_ms][2]

@info "Loading Data Files"

b0FileName = params_general[:fullPathMapScan];

# filename for preprocessed data (remove oversampling, permute dimensions wrt MRIReco)
processedFileName = params_general[:fullPathProcessedMapScan]

if params_general[:do_process_map_scan]

    # Set the data file name (Change this for your own system)
    dataFileCartesian = ISMRMRDFile(b0FileName)

    # read in the raw data from the ISMRMRD file into a RawAcquisitionData object
    r = RawAcquisitionData(dataFileCartesian)

    # does not change anything...
    # r.params["reconFOV"] = [230, 230, 2]

    # Preprocess Data and save!
    preprocess_cartesian_data(r::RawAcquisitionData, true; fname = processedFileName)

end

# Load preprocessed data!
dataFileNew = ISMRMRDFile(processedFileName)


rawDataNew = RawAcquisitionData(dataFileNew)
acqDataCartesian = AcquisitionData(rawDataNew, estimateProfileCenter = true)

# Define coils and slices
nCoils = size(acqDataCartesian.kdata[1], 2)
nSlices = numSlices(acqDataCartesian)


sliceIndexArray = get_slice_order(nSlices, isSliceInterleaved = true)
# shift FOV to middle :) 
#TODO: in MRIReco v0.7, try: correctOffset(acqDataCartesian, [0 -20 0])
shiftksp!(acqDataCartesian, params_general[:fovShift])
#changeFOV!(acqDataCartesian,[1.5,1.5])

## Don't have to recalculate sense maps for both scans but possibly it could make a
#  difference in Diffusion scans

@info "Calculating Sense Maps"

# from reading docstring, min thresholding
cartesian_sensitivity = espirit(acqDataCartesian, (6, 6), 30, eigThresh_1 = 0.01, eigThresh_2 = 0.9)

# from Lars (7T spirals)
#cartesian_sensitivity = espirit(acqDataCartesian,(6,6),30,eigThresh_1=0.02, eigThresh_2=0.98)
# from Alexander (Phantom?)
#cartesian_sensitivity = espirit(acqDataCartesian,(4,4),12,eigThresh_1=0.01, eigThresh_2=0.98)

# normalize for consistency with saving/loading and better ranges of reconstruction values
cartesian_sensitivity /= maximum(abs.(cartesian_sensitivity))
sensitivity = cartesian_sensitivity

resolution_mm = fieldOfView(acqDataCartesian) ./ size(sensitivity)[1:3]
resolution_mm[3] = fieldOfView(acqDataCartesian)[3] * (1 + params_general[:sliceDistanceFactor_percent] / 100.0); # for 2D only, since FOV[3] is slice thickness then, but gap has to be observed


# save SENSE maps
if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    # TODO: use correct slice order everywhere, e.g., when saving/loading maps for spiral recon
    save_map(params_general[:fullPathSaveSense], sensitivity[:, :, sliceIndexArray, :], resolution_mm; doSplitPhase = true)
end

plot_sense_maps(sensitivity, nCoils)

#acqDataCartesian.traj[1].cartesian = false
#acqDataCartesian.traj[2].cartesian = false
## Parameter dictionary definition for reconstruction

@info "Setting Parameters"
params_cartesian = Dict{Symbol,Any}() # instantiate dictionary
params_cartesian[:reco] = "multiCoil" # choose multicoil reconstruction

# TODO: make recon size and FOV variable!
params_cartesian[:reconSize] = (acqDataCartesian.encodingSize[1], acqDataCartesian.encodingSize[2]) # set recon size to be the same as encoded size
params_cartesian[:regularization] = "L2" # choose regularization for the recon algorithm
params_cartesian[:λ] = 1.e-2 # recon parameter (there may be more of these, need to dig into code to identify them for solvers other than cgnr)
params_cartesian[:iterations] = params_general[:nReconIterations] # number of CG iterations
params_cartesian[:solver] = "cgnr" # inverse problem solver method
params_cartesian[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false) # turn on store solutions if you want to see the reconstruction convergence (uses more memory)
params_cartesian[:senseMaps] = ComplexF32.(sensitivity) # set sensitivity map array


## Call the reconstruction function

@info "Performing Reconstruction"
@time cartesianReco = reconstruction(acqDataCartesian, params_cartesian)

# save Map recon (multi-echo etc.)
if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    save_map(params_general[:fullPathsave_mapRecon], cartesianReco.data[:, :, sliceIndexArray, :, :, :], resolution_mm; doSplitPhase = true)
end

## Calculate B0 maps from the acquired images (if two TEs)

slices = 1:length(sliceIndexArray)

@info "Calculating B0 Maps"
# b0Maps = calculate_b0_maps(cartesianReco.data,slices, TE1, TE2)
b0Maps = estimateB0Maps(cartesianReco.data, slices, TE1, TE2, true; β = params_general[:b0mapSmoothBeta], reltol = 1e-4)

# save B0 map
if params_general[:do_save_recon] # TODO: include elements to save as tuple, e.g., ["b0", "sense", "recon"], same for load
    save_map(params_general[:fullPathSaveB0], b0Maps[:, :, sliceIndexArray], resolution_mm; doNormalize = false) # no normalization, we want absolute values for offres maps
end


if params_general[:do_plot_recon]
    @info "Plotting Cartesian Results (Sensitivity Maps and B0 Maps)"
    pygui(true) # Leave this code till we need plotting.
    # plot_sense_maps(sensitivity,nCoils)
    plot_reconstruction(cartesianReco[:, :, :, 1], 1:size(cartesianReco, 3), b0Maps, isSliceInterleaved = true, rotateAngle = 270)
end

# cleanup unused file
if !params_general[:do_save_processed_map_scan] && params_general[:do_process_map_scan]
    rm(processedFileName)
end

@info "Successfully Completed CartesianReconstruction"

