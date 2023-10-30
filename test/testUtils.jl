function test_getSliceOrder()

end

function test_syncTrajAndData!()

end

function test_do_k0_correction!()

end

function test_adjustHeader!()

end

function test_checkAcquisitionNodes!()

end

function test_validateSiemensMRD!()

end

function test_validateAcqData!()

end

function test_preprocessCartesianData()

end

function test_removeOversampling!()

end

function test_mergeRawInterleaves()

end

function test_applyGIRF!()

end

function test_applyK0!()

end

function test_saveMap()

end

function test_loadMap()

end

function test_shift_kspace!()

end

function testUtils(N=32)
    @testset "Utilities" begin
        
        test_adjustHeader!()
        test_applyGIRF!()
        test_applyK0!()
        test_checkAcquisitionNodes!()
        test_do_k0_correction!()
        test_getSliceOrder()
        test_loadMap()
        test_mergeRawInterleaves()
        test_preprocessCartesianData()
        test_removeOversampling!()
        test_saveMap()
        test_shift_kspace!()
        test_syncTrajAndData!()
        test_validateAcqData!()
        test_validateSiemensMRD!()

    end
end

testUtils()