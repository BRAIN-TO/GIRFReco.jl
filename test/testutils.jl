function test_get_slice_order()

    slice_positions = [-7, -3, 1, 5, 9, -9, -5, -1, 3, 7]

    N = 32
    I_slice = shepp_logan(N)
    I = I_slice .* ones(N,N,10)

    # simulation parameters
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Radial"
    params[:numProfiles] = floor(Int64, pi/2*N)
    params[:numSamplingPerProfile] = 2*N

    acqData = simulation(I, params)
    rawData = RawAcquisitionData(acqData)

    tempNum = 0
    for ii = 1:params[:numProfiles]:500
        tempNum = tempNum + 1
        for jj = 1:params[:numProfiles]

            rawData.profiles[ii + jj - 1].head.position = (0.0,0.0,slice_positions[tempNum])

        end

    end

    @test get_slice_order(rawData, 10,1,50) == [6,1,7,2,8,3,9,4,10,5]

end

function test_sync_traj_and_data!()

    # General strategy for this test: Get traj from simulated object, interpolate it onto different grid then reinterpolate and verify that things still work
    traj_input = trajectory(Float64,"Spiral",10,400)

    N = 32
    I_slice = shepp_logan(N)
    I = I_slice .* ones(N,N,10)

    # simulation parameters
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Spiral"
    params[:numProfiles] = 10
    params[:numSamplingPerProfile] = 2000

    acqData = simulation(I, params)
    rawData = RawAcquisitionData(acqData)
    rawData_old = deepcopy(rawData)

    sampTimes = sync_traj_and_data!(rawData,traj_input,2000,1)

    @test mean(√,abs2.(rawData.profiles[1].traj[:] .- rawData_old.profiles[1].traj[:])) < 1e-2

end

function test_calculate_b0_maps()

    N = 32
    T = ComplexF32
    nCh = 4
    slices = 5
    nEchos = 2
    TE = 7.0

    medata = ones(ComplexF32,N,N,slices,2,nCh)
    medata[:,:,:,1,:] .*= exp(-1im .* 2*pi * 0.01)
    medata[:,:,:,2,:] .*= exp(-1im .* 2*pi * 0.02)

    b0_maps = calculate_b0_maps(medata,1:slices,2.0,4.0)

    @test median(b0_maps) + 2*pi*0.01 ./ 2e-3 < 1e-5 

end

function test_do_k0_correction!()

    N = 32
    I = shepp_logan(N)

    # simulation parameters
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Radial"
    params[:numProfiles] = floor(Int64, pi/2*N)
    params[:numSamplingPerProfile] = 2*N

    acqData = simulation(I, params)
    rawData = RawAcquisitionData(acqData)
    rawData.params["encodedSize"] = [rawData.params["encodedSize"][1],rawData.params["encodedSize"][2],1]
    rawData_orig = deepcopy(rawData)
    phase_mod = pi .* 0.1 .* ones(ComplexF32, length(rawData.profiles[1].data),2)
    do_k0_correction!(rawData,phase_mod,1)
    acqData2 = AcquisitionData(rawData)

    @test median(angle.(acqData2.kdata[1]) .- angle.(acqData.kdata[1])) .- pi .* 0.1 < 1e-4

end

function test_adjust_header!()

    N = 32
    I = shepp_logan(N)

    # simulation parameters
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Radial"
    params[:numProfiles] = floor(Int64, pi/2*N)
    params[:numSamplingPerProfile] = 2*N

    # do simulation
    acqData = simulation(I, params)

    rawData = RawAcquisitionData(acqData)

    adjust_header!(rawData,[N,N,1],2*N,1,true)

    @test rawData.profiles[1].head.discard_post == 0
    @test rawData.profiles[1].head.discard_pre == 0

end

function test_check_acquisition_nodes!()

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

    acqData.traj[1].nodes .*= 1.5

    check_acquisition_nodes!(acqData)

    @test abs.(maximum(acqData.traj[1].nodes,dims = [1,2]))[1] .< 0.51

end

function test_validate_siemens_mrd!()

    N = 32
    I = shepp_logan(N)

    # simulation parameters
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Radial"
    params[:numProfiles] = floor(Int64, pi/2*N)
    params[:numSamplingPerProfile] = 2*N

    # do simulation
    acqData = simulation(I, params)

    rawData = RawAcquisitionData(acqData)

    rawData.params["encodedFOV"] = 800

    validate_siemens_mrd!(rawData)

    @test rawData.params["encodedFOV"] == 0.8

end

function test_validate_acq_data!()

    N = 128
    T = ComplexF32
    nCh = 4
    nEchos = 2
    TE = 7.0
    
    x = T.(shepp_logan(N))

    rmap = 0.05*abs.(x)
    TEnum = Float64.(collect(TE:TE:TE*nEchos))

    coilsens = T.(birdcageSensitivity(N, nCh, 4.))
    params = Dict{Symbol,Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N
    params[:r2map] = rmap
    params[:T_echo] = TEnum
    params[:seqName] = "ME"
    params[:refocusingAngles] = Float64.(repeat([pi], length(TEnum)))
    params[:senseMaps] = coilsens

    acqData = simulation(real(x), params)

    validate_acq_data!(acqData)

    @test size(acqData.kdata) == (2,1,1)

end

function test_preprocess_cartesian_data()

    N = 32
    I = shepp_logan(N)

    # simulation parameters
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N

    # do simulation
    acqData = simulation(I, params)
    acqData.traj[1].nodes .*= 1.5

    rawData = RawAcquisitionData(acqData)
    rawData.params["encodedSize"] = [rawData.params["encodedSize"][1],rawData.params["encodedSize"][2],1]

    acqData2 = preprocess_cartesian_data(rawData,false)

    @test length(acqData2.kdata[1]) == length(acqData.kdata[1])÷2
    @test maximum(abs.(acqData2.traj[1].nodes),dims=[1,2])[1] <= 0.5 

end

function test_remove_oversampling!()

    N = 32
    I = shepp_logan(N)

    # simulation parameters
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N

    # do simulation
    acqData = simulation(I, params)

    rawData = RawAcquisitionData(acqData)

    remove_oversampling!(rawData)

    @test rawData.params["encodedSize"] == [N÷2,N]
    @test length(rawData.profiles[1].data) == N÷2

end

function test_save_and_load_map()
    N = 128
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
    N = 32
    I = shepp_logan(N)
    I[16,16] = 500

    # simulation parameters
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Radial"
    params[:numProfiles] = floor(Int64, pi/2*N)
    params[:numSamplingPerProfile] = 2*N

    acqData = simulation(I, params)
    
    # reco parameters
    params = Dict{Symbol, Any}()
    params[:reco] = "direct"
    params[:reconSize] = (N,N)
    params[:alpha] = 1.75
    params[:m] = 4.0
    params[:K] = 28
    Ireco = reconstruction(acqData, params)

    shift_kspace!(acqData,[5,5])

    Ireco_shifted = reconstruction(acqData,params)

    @test (abs.(Ireco[16,16]) > 300) & (abs.(Ireco[21,21]) < 300)
    @test (abs.(Ireco_shifted[16,16]) < 300) & (abs.(Ireco_shifted[21,21]) > 300)
    @test abs.(Ireco_shifted[21,21]) .- abs.(Ireco[16,16]) < 1e1
end

function test_utils(N=32)
    @testset "Utilities" begin
        
        test_adjust_header!()
        test_calculate_b0_maps()
        test_check_acquisition_nodes!()
        test_do_k0_correction!()
        test_get_slice_order()
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
