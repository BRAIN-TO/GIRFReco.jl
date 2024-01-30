function test_merge_raw_interleaves()

end

function test_apply_girf!(acq_data,girf_applier)

    acq_data_ref = deepcopy(acq_data)
    apply_girf!(acq_data,girf_applier)

    @test  mean(abs2.(fft(acq_data.traj[1].nodes[1,:]))./abs2.(fft(acq_data_ref.traj[1].nodes[1,:]))) - 0.98 < 1e-4

end

function test_apply_k0!(acq_data,girf_applier)

end

function test_integration()

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