using PyPlot, HDF5, MRIReco, LinearAlgebra, DSP, Images, FourierTools, ROMEO

## Plotting Functions
function plotSenseMaps(sense,n_channels)

    # Create figure and plot the sensitivity maps for each coil composed as a collage

    # Magnitude maps
    figure("Sensitivity Map Magnitude"); clf(); for ch in 1:n_channels; subplot(8,4,ch); imshow((abs.(sense[:,:,1,ch]))); end;
    subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
    gcf()

    # Phase maps
    figure("Sensitivity Map Phase"); clf(); for ch in 1:n_channels; subplot(8,4,ch); imshow(ROMEO.unwrap(angle.(sense[:,:,1,ch]))); end;
    subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
    gcf()

end

## General plotting function for the reconstruction
function plotReconstruction(images, slices, b0)

    # Slice ordering check (show the correct order of slices as the images are read in not in geometrically sequential order)
    indexArray = slices

    # Plot magnitude images (normalize)
    figure("Magnitude Images")
    absData = abs.(images[:,:,slices,1,1])./maximum(abs.(images[:,:,slices,1,1]))
    absMosaic = mosaicview(absData, nrow=Int(floor(sqrt(length(slices)))),npad=5,rowmajor=true, fillvalue=0)

    imshow(absMosaic,cmap="gray")
    colorbar()

    gcf().suptitle("|Images|")

    # Plot phase images
    figure("Phase Images")

    phaseData = mapslices(x ->ROMEO.unwrap(x),angle.(images[:,:,slices,1,1]),dims=(1,2))
    phaseMosaic = mosaicview(phaseData,nrow=Int(floor(sqrt(length(slices)))),npad=5,rowmajor=true, fillvalue=0)

    imshow(phaseMosaic, cmap="inferno",vmax = 3*pi,vmin=-3*pi)
    colorbar()

    gcf().suptitle("∠Images")

    # Plot B0 maps
    figure("B₀ Map Images")

    b0Mosaic = mosaicview(b0[:,:,slices],nrow=Int(floor(sqrt(length(slices)))),npad = 5, rowmajor=true, fillvalue=0)

    imshow(b0Mosaic, cmap="inferno",vmax=500,vmin=-500)
    colorbar()

    gcf().suptitle("B₀ Maps")

end

## function prints out all profiles in the acquisition to check consistency with ISMRMRD file
function checkProfiles(rawData)

    numProfiles2 = 128 # Set to the number of profiles that you would like to see

    for l = 1:numProfiles2
        figure("Profile $l")
        subplot(2,1,1)
        plot(abs.(rawData.profiles[l].data[:,1]))
        subplot(2,1,2)
        plot(angle.(rawData.profiles[l].data[:,1]))
    end

end

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
dataFileCartesian = ISMRMRDFile("GIRFReco/data/Fieldmaps/fieldMap_30_2.h5")

# read in the raw data from the ISMRMRD file into a RawAcquisitionData object
rawDataCartesian = RawAcquisitionData(dataFileCartesian)

# Convert rawAcquisitionData object to an AcquisitionData object (these can be reconstructed)
acqDataCartesian = AcquisitionData(rawDataCartesian,estimateProfileCenter=true)

## Properly arrange data from the converted siemens file

@info "Correcting conversion errors"
acqDataCartesian.fov = acqDataCartesian.fov/1000 # Units... MRIReco.jl takes units of [m] and siemens gives in [mm]

# Fix the FOV (can be set incorrectly)
acqDataCartesian.fov[1] = 0.22

# Need to permute the dimensions of kdata to match the convention of MRIReco.jl
permutedims(acqDataCartesian.kdata,[3,2,1])

# Define coils and slices
nCoils = size(acqDataCartesian.kdata[1],2)
nSlices = numSlices(acqDataCartesian)

## Don't have to recalculate sense maps for both scans but possibly it could make a
#  difference in Diffusion scans

@info "Calculating Sense Maps"
@time senseCartesian = espirit(acqDataCartesian,(6,6),30,eigThresh_1=0.05, eigThresh_2=0.98)
sensitivity = senseCartesian;

## Resize sense maps to match encoding size of data matrix
sensitivity = imresize(senseCartesian,(acqDataCartesian.encodingSize[1],acqDataCartesian.encodingSize[2],nSlices,nCoils));
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
