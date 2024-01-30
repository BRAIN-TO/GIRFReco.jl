
imsize = 64

#normal image
N = imsize ÷ 2
img = shepp_logan(N)
msk = zeros(N, N)
msk[findall(x -> x != 0, img)] .= 1

#larger image
img2 = shepp_logan(imsize)
msk2 = zeros(imsize, imsize)
msk2[findall(x -> x != 0, img2)] .= 1

# coil sensitivites for normal image size
smaps = birdcageSensitivity(N, 8, 1.5)
snorm = sqrt.(sum(abs.(smaps) .^ 2, dims = 4))
for i = 1:8
    smaps[:, :, 1, i] .= msk .* smaps[:, :, 1, i] ./ snorm[:, :, 1, 1]
end

# coil sensitivites for larger image size
smaps2 = birdcageSensitivity(imsize, 8, 1.5)
snorm2 = sqrt.(sum(abs.(smaps2) .^ 2, dims = 4))
for i = 1:8
    smaps2[:, :, 1, i] .= msk2 .* smaps2[:, :, 1, i] ./ snorm2[:, :, 1, 1]
end

# simulation for normal image size
params = Dict{Symbol,Any}()
params[:simulation] = "fast"
params[:trajName] = "Cartesian"
params[:numProfiles] = floor(Int64, N)
params[:numSamplingPerProfile] = N
params[:senseMaps] = smaps

acqData = simulation(img, params)
acqData = MRIReco.sample_kspace(acqData, 2.0, "poisson", calsize = 15)

ksize = (6, 6) # kernel size
ncalib = 15 # number of calibration lines
eigThresh_1 = 0.02  # threshold for picking singular vectors of calibration matrix
eigThresh_2 = 0.95  # threshold for eigen vector decomposition in image space

emaps = espirit(acqData, ksize, ncalib, (imsize, imsize), eigThresh_1 = eigThresh_1, eigThresh_2 = eigThresh_2)

# simulation for larger image size
params = Dict{Symbol,Any}()
params[:simulation] = "fast"
params[:trajName] = "Cartesian"
params[:numProfiles] = floor(Int64, imsize)
params[:numSamplingPerProfile] = imsize
params[:senseMaps] = smaps2

acqData2 = simulation(img2, params)
acqData2 = MRIReco.sample_kspace(acqData2, 2.0, "poisson", calsize = 15)
emaps2 = espirit(acqData2, ksize, ncalib, (imsize, imsize), eigThresh_1 = eigThresh_1, eigThresh_2 = eigThresh_2)

# evaluate error only on the support of smaps
for i = 1:8
    emaps2[:, :, 1, i] = msk2 .* emaps2[:, :, 1, i]
    emaps[:, :, 1, i] = msk2 .* emaps[:, :, 1, i]
end

err = norm(vec(emaps2) .- vec(emaps)) / norm(vec(emaps))

function getCalib(acqData, imsize)

    T = Float64

    nxAcq, nyAcq = acqData.encodingSize[1:2]
    nx, ny = imsize
    match_acq_size = all(imsize .== acqData.encodingSize[1:2])

    #  Force maps to be at least same size as acqData.encodingSize.
    if (nxAcq >= nx || nyAcq >= ny)
        nx, ny = acqData.encodingSize[1:2]
    end

    numChan, numSl = numChannels(acqData), numSlices(acqData)
    calib = zeros(Complex{T}, ncalib, ncalib, numChan)

    idx = match_acq_size ? acqData.subsampleIndices[1] : findIndices((nx, ny), (nxAcq, nyAcq))[acqData.subsampleIndices[1]]

    # form zeropadded array with kspace data
    kdata = zeros(Complex{T}, nx * ny, numChan)
    for coil = 1:numChan
        kdata[idx, coil] .= kData(acqData, 1, coil, 1)
    end
    kdata = reshape(kdata, nx, ny, numChan)

    calib = MRIReco.crop(kdata, (ncalib, ncalib, numChan))

    return calib

end

