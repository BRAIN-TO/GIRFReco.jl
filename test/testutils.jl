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

function test_save_and_load_map()
    N = 256
    I = shepp_logan(N)

    # simulation parameters
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Radial"
    params[:numProfiles] = floor(Int64, pi/2*N)
    params[:numSamplingPerProfile] = 2*N

    # do simulation
    acqData = simulation(I, params)

    # reco parameters
    params = Dict{Symbol, Any}()
    params[:reco] = "direct"
    params[:reconSize] = (N,N)
    params[:alpha] = 1.75
    params[:m] = 4.0
    params[:K] = 28
    Ireco = reconstruction(acqData, params)

    res_x = fieldOfView(acqData)[1] ./ encodingSize(acqData)[1]
    res_y = fieldOfView(acqData)[2] ./ encodingSize(acqData)[2]
    res_z = fieldOfView(acqData)[3]
    resolution_mm =(res_x, res_y, res_z)

    # Complex image save & load
    save_map(joinpath(tmpdir, "testComplexIO.nii"), Ireco, resolution_mm; do_split_phase = true, do_normalize = false)
    IrecoLoaded = load_map(joinpath(tmpdir, "testComplexIO.nii"), do_split_phase = true )
    @test (norm(vec(Ireco)-vec(IrecoLoaded))/norm(vec(Ireco))) < 1e-5

    # Magnitude image save & load
    save_map( joinpath(tmpdir, "testMagIO.nii"), abs.(Ireco), resolution_mm; do_split_phase = false, do_normalize = false)
    IrecoLoadedMag = load_map(joinpath(tmpdir, "testMagIO.nii"), do_split_phase = false )
    @test (norm(vec(abs.(Ireco))-vec(IrecoLoadedMag))/norm(vec(abs.(Ireco)))) < 1e-5
end

function test_shift_kspace!()

end

function test_utils(N=32)
    @testset "Utilities" begin
        
        test_adjust_header!()
        test_apply_girf!()
        test_apply_k0!()
        test_check_acquisition_nodes!()
        test_do_k0_correction!()
        test_get_slice_order()
        test_merge_raw_interleaves()
        test_preprocess_cartesian_data()
        test_remove_oversampling!()
        test_save_and_load_map()
        test_shift_kspace!()
        test_sync_traj_and_data!()
        test_validate_acq_data!()
        test_validate_siemens_mrd!()

    end
end

test_utils()
