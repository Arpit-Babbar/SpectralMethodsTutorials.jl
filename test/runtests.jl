using SpectralMethodsTutorials
using SpectralMethodsTutorials: test_ch3, test_dir
using Test

include(test_dir() * "/ch6.jl")

@testset "SpectralMethodsTutorials.jl" begin
    test_ch3()
    test_ch6()
end
