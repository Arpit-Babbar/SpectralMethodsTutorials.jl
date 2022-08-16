using SpectralMethodsTutorials
using SpectralMethodsTutorials: test_dir
using Test

include(test_dir() * "/ch6.jl")
include(test_dir() * "/ch3.jl")

@testset "SpectralMethodsTutorials.jl" begin
    test_ch3()
    test_ch6()
end
