function test_merge_raw_interleaves()

end

function test_apply_girf!()

end

function test_apply_k0!()

end

function test_integration()
    @testset "Integration" begin
        test_merge_raw_interleaves()
        test_apply_girf!()
        test_apply_k0!()
    end
end

test_integration()