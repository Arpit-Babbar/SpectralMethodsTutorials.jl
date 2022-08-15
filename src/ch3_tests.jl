using Test
(
using SpectralMethodsTutorials: barycentric_weights,
                                lagrange_interpolation,
                                lagrange_derivative,
                                differentiation_matrix
)

function test_weights()
   x = [1.0, 2.0, 3.0]
   @test barycentric_weights(x) ≈ [0.5, -1.0, 0.5]
end

function test_interpolation()
   xg = [1.0, 2.0]
   f_values = [2.0, 4.0]
   wg = barycentric_weights(xg)
   @test lagrange_interpolation(1.5, xg, f_values, wg) ≈ 3.0
end

function test_diff1()
   xg = [1.0, 2.0, 3.0]
   xg_ = [1.1, 2.1, 3.1]
   f_vals = 2.0*xg
   wg = barycentric_weights(xg)
   p(x) = lagrange_derivative(x, xg, f_vals, wg)
   @test p.(xg) ≈ [2.0,2.0,2.0]
   @test p.(xg_) ≈ [2.0,2.0,2.0]
end

function test_diff2()
   xg = [1.0, 2.0, 3.0, 4.0]
   xg_ = [1.3, 2.3, 3.3]
   f_vals = xg.^3
   wg = barycentric_weights(xg)
   p(x) = lagrange_derivative(x, xg, f_vals, wg)
   @test p.(xg) ≈ 3.0*xg.^2
   @test p.(xg_) ≈ 3.0*xg_.^2
end

function test_diff_mat1()
   xg = [1.0, 2.0, 3.0]
   D = differentiation_matrix(1, xg)
   f_vals = 2.0*xg
   @test D*f_vals ≈ [2.0,2.0,2.0]
end

function test_diff_mat2()
   xg = [1.0, 2.0, 3.0, 4.0]
   D = differentiation_matrix(1, xg)
   f_vals = xg.^3
   @test D*f_vals ≈ 3.0*xg.^2
end

function test_diff_mat3()
   xg = [1.0, 2.0, 3.0, 4.0]
   D = differentiation_matrix(2, xg)
   f_vals = xg.^3
   @test D*f_vals ≈ 6.0*xg
end

function test_diff_mat4()
   xg = [1.0, 2.0, 3.0, 4.0]
   D = differentiation_matrix(3, xg)
   f_vals = xg.^3
   nd = length(xg)
   @test D*f_vals ≈ 6.0*ones(nd)
end

function test_ch3()
   test_weights()
   test_interpolation()
   test_diff1()
   test_diff2()
   test_diff_mat1()
   test_diff_mat2()
   test_diff_mat3()
   test_diff_mat4()
end

@testset begin
   test_ch3()
end