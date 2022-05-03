using DelimitedFiles, MRIReco, PyPlot, Dierckx, MAT, DSP

export read_gradient_txt_file

function read_gradient_txt_file(fileName, reconSize, delay)
    
    gradientData = readdlm(fileName,'\n')

    ## Read in the header data of the gradient text file (lines 1 to 21)
    dataDict = Dict{Symbol,Any}()
    dataDict[:versionNr] = gradientData[1]
    dataDict[:numSamples] = gradientData[2]
    dataDict[:dwellTime] = gradientData[3] # [seconds]
    dataDict[:samplesPerInterleave] = gradientData[4]
    dataDict[:numInterleaves] = gradientData[5]
    dataDict[:numDims] = gradientData[6]
    dataDict[:timeToCenterKSpace] = gradientData[7] # [seconds]
    dataDict[:acqDuration] = gradientData[8]
    dataDict[:samplesPerAcq] = gradientData[9]
    dataDict[:numAcquisitions] = gradientData[10]
    dataDict[:acqTR] = gradientData[11]
    dataDict[:gradientAcqStartDelay] = gradientData[12]
    dataDict[:echoTimeShiftSamples] = gradientData[13]
    dataDict[:FOV] = gradientData[14:16] # [m]
    dataDict[:voxelDims] = gradientData[17:19] # [m]
    dataDict[:gradientStrengthFactor] = gradientData[20] # [mT/m]
    dataDict[:isBinary] = gradientData[21]
    dataDict[:gamma] = 42577.478 # [Hz/mT] CAN CHANGE
    dataDict[:fieldStrength] = 3 # [T] CAN CHANGE

    #print(dataDict)

    ## reading and data scaling of gradient data
    dataDict[:gradData] = gradientData[22:end]
    interleaveGradArray = dataDict[:gradientStrengthFactor]*reshape(dataDict[:gradData],dataDict[:samplesPerInterleave],dataDict[:numInterleaves],dataDict[:numDims]) #[mT/m]

    plannedTimes = dataDict[:dwellTime].*(0:(dataDict[:samplesPerInterleave]-1))
    delayedTimes = plannedTimes .- delay .- dataDict[:dwellTime]./2 # seconds (dwellTime/2 compensates for integration)

    interleaveGradArrayFlexible = Array{Float64,3}(undef,size(interleaveGradArray))

    #print(size(interleaveGradArrayFlexible))

    ## Loop over all of the unique excitation trajectories and create an interpolant of the gradient
    for dim in 1:dataDict[:numDims]

        for l in 1:dataDict[:numInterleaves]

            #print((dim,l),"\n")

            sp = Spline1D(plannedTimes,interleaveGradArray[:,l,dim],w=ones(length(plannedTimes)), k=1, bc="zero", s=0.0)

            # evaluate the interpolant at the sampling times of the kspace data
            interleaveGradArrayFlexible[:,l,dim] = sp(delayedTimes)

            #print(interleaveGradArrayFlexible[:,l,dim][end],"\n")

        end

    end

    ## cumulative summation and numerical integration of the gradient data, resulting in the kspace trajectory
    kSpaceTrajArrayFlexible = dataDict[:gamma]*dataDict[:dwellTime]*cumsum(interleaveGradArrayFlexible,dims=1) # [rad/m]

    ## Conversion to the trajectory scaling convention in MRIReco.jl
    #  Currently only 2d Trajectories
    convertedKSpaceArrayFlexible = kSpaceTrajArrayFlexible
    convertedKSpaceArrayFlexible[:,:,1] *= dataDict[:FOV][1] ./ reconSize[1]
    convertedKSpaceArrayFlexible[:,:,2] *= dataDict[:FOV][2] ./ reconSize[2]

    ## Construction of the trajectory object ##

    ## Reshaping of the array to the format expected by the Trajectory constructor in MRIReco.jl
    # - dim 1 = kspace dimension
    # - dim 2 = kspace position (with interleaves/profiles arranged consecutively)
    permutedTrajectory = permutedims(reshape(convertedKSpaceArrayFlexible,dataDict[:samplesPerInterleave]*dataDict[:numInterleaves],dataDict[:numDims]),[2,1])

    ## Construction of the trajectory
    # - Note: timing vectors are automatically generated - seems to be consistent with the dwell time
    trajectoryObject = Trajectory(permutedTrajectory,dataDict[:numInterleaves],dataDict[:samplesPerInterleave],TE=dataDict[:echoTimeShiftSamples],AQ=dataDict[:acqDuration], numSlices=9, cartesian=false,circular=false)

    return trajectoryObject

end

function testGradReader()

    ## Testing
    gradFile = "data/Gradients/gradients523.txt"

    ##
    kSpaceTrajectory_2 = read_gradient_txt_file(gradFile,(200,200),0.00000)

    ##
    pulledTrajectory21 = kspaceNodes(kSpaceTrajectory_2)[1,:]
    pulledTrajectory22 = kspaceNodes(kSpaceTrajectory_2)[2,:]

    #
    fig = figure(234, figsize=(10,10))
    ax = fig.gca()
    ax.scatter(pulledTrajectory21,pulledTrajectory22, label="Nominal")
    xlabel("kx")
    ylabel("ky")
    title("K-space Center")
    xlim((-0.05,0.05))
    ylim((-0.05,0.05))
    legend()

end