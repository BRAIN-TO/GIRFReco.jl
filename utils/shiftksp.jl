function shiftksp!(acqData,shift)

    numSl = numSlices(acqData)
    numRep, numContr = numRepetitions(acqData), numContrasts(acqData)

    smat = prod(exp.(1im .* acqData.traj[1].nodes[:,acqData.subsampleIndices[1]] .* shift .* 2 .* pi),dims=1)

    for slice = 1:numSl
        for contr = 1:numContr
            for rep = 1:numRep                    
                acqData.kdata[contr,slice,rep] = acqData.kdata[contr,slice,rep] .* smat'                    
            end
        end
    end

end

function changeFOV!(acqData,factor)
    # assume vectorized nodes

    numTr = length(acqData.traj)

    for tr = 1:numTr                    
        acqData.traj[tr].nodes .= acqData.traj[tr].nodes .* factor
                    
    end


end
    