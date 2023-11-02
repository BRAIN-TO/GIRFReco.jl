function test_get_slice_order()

end

function test_sync_traj_and_data!()

end

function test_do_k0_correction!()

end

function test_adjust_header!()

end

function test_check_acquisition_nodes!()

end

function test_validate_siemens_mrd!()

end

function test_validate_acq_data!()

end

function test_preprocess_cartesian_data()

end

function test_remove_oversampling!()

end

function test_merge_raw_interleaves()

end

function test_apply_girf!()

end

function test_apply_k0!()

end

function test_save_map()

end

function test_load_map()

end

function test_shift_kspace!()

end

function testUtils(N=32)
    @testset "Utilities" begin
        
        test_adjust_header!()
        test_apply_girf!()
        test_apply_k0!()
        test_check_acquisition_nodes!()
        test_do_k0_correction!()
        test_get_slice_order()
        test_load_map()
        test_merge_raw_interleaves()
        test_preprocess_cartesian_data()
        test_remove_oversampling!()
        test_save_map()
        test_shift_kspace!()
        test_sync_traj_and_data!()
        test_validate_acq_data!()
        test_validate_siemens_mrd!()

    end
end

testUtils()