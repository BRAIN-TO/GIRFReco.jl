function test_gradient_reader()

    grad_filename = "data/gradients508.txt"
    read_gradient_text_file(grad_filename,[200,200],0.0)
    @test true

end

function test_io()
    @testset "I/O" begin
        
        test_gradient_reader()

    end
end

test_io()