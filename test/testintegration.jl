function test_merge_raw_interleaves()

    N = 32
    I_slice = shepp_logan(N)
    I = I_slice .* ones(N,N,10)

    # simulation parameters
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Spiral"
    params[:numProfiles] = 4
    params[:numSamplingPerProfile] = 15650

    acqData = simulation(I, params)
    acqData.encodingSize = (N,N)

    rawData = RawAcquisitionData(acqData)
    rawData.params["encodedSize"] = [N,N,1]
    
    # Set up raw data splitting
    rawData_il1 = deepcopy(rawData)
    rawData_il2 = deepcopy(rawData)
    rawData_il3 = deepcopy(rawData)
    rawData_il4 = deepcopy(rawData)

    deleteat!(rawData_il1.profiles,1:4:40)
    deleteat!(rawData_il2.profiles,2:4:40)
    deleteat!(rawData_il3.profiles,3:4:40)
    deleteat!(rawData_il4.profiles,4:4:40)
    
    fout_il1 = ISMRMRDFile("$tmpdir/temp_il1")
    fout_il2 = ISMRMRDFile("$tmpdir/temp_il2")
    fout_il3 = ISMRMRDFile("$tmpdir/temp_il3")
    fout_il4 = ISMRMRDFile("$tmpdir/temp_il4")

    save(fout_il1,rawData_il1)
    save(fout_il2,rawData_il2)
    save(fout_il3,rawData_il3)
    save(fout_il4,rawData_il4)

    filenames_interleaves = ["$tmpdir/temp_il1","$tmpdir/temp_il2","$tmpdir/temp_il3","$tmpdir/temp_il4"]
    is_single_interleave = false

    params_spiral = Dict{Symbol,Any}()
    params_spiral[:recon_size] = (N,N,1)
    params_spiral[:interleave] = 1
    params_spiral[:num_samples] = 15650
    params_spiral[:delay] = 0.00000 # naive delay correction
    params_spiral[:interleave_data_filenames] = filenames_interleaves
    params_spiral[:traj_filename] = "data/gradients508.txt"
    params_spiral[:excitations] = [1,2,3,4,5,6,7,8,9,10]
    params_spiral[:do_multi_interleave] = !is_single_interleave
    params_spiral[:do_odd_interleave] = false
    params_spiral[:num_interleaves] = is_single_interleave ? 1 : length(params_spiral[:interleave_data_filenames]) # one interleaf per file, count files, if filenames are array of strings (not only one string)
    params_spiral[:single_slice] = false
    
    acq_out = merge_raw_interleaves(params_spiral,false)

    @test size(acq_out.traj[1].nodes,2) == 62600
    @test √(sum(abs2,acq_out.traj[1].nodes[:,62600] - [-1.1819346;-2.884264])) < 1e-6 # known gt

end

function test_apply_girf!(acq_data,girf_applier)

    acq_data_ref = deepcopy(acq_data)
    apply_girf!(acq_data,girf_applier)

    # these should be pretty unique checks against ground-truth
    @test  mean(abs2.(fft(acq_data.traj[1].nodes[1,:]))./abs2.(fft(acq_data_ref.traj[1].nodes[1,:]))) - 0.98 < 1e-4

end

function test_apply_k0!(acq_data,girf_applier)

    acq_data_ref = deepcopy(acq_data)
    apply_k0!(acq_data,girf_applier)

    # these should be pretty unique checks against ground-truth
    @test  mean(angle.(acq_data_ref.kdata[1]) .- angle.(acq_data.kdata[1])) - 0.00118 < 1e-5
end

function test_integration()

    # load generated AcquisitionData objects and GirfApplier objects to run the test set
    acq_data = load_object("data/acq_data.jld2")
    girf_applier_k1 = load_object("data/girf_app_k1.jld2")
    girf_applier_k0 = load_object("data/girf_app_k0.jld2")

    @testset "Integration" begin
        test_merge_raw_interleaves()
        test_apply_girf!(acq_data,girf_applier_k1)
        test_apply_k0!(acq_data,girf_applier_k0)
    end

end

test_integration()