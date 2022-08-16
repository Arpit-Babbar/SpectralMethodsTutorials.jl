using Test
(
using SpectralMethodsTutorials: barycentric_weights,
                                lagrange_interpolation,
                                lagrange_derivative,
                                differentiation_matrix,
                                chebyshev_lobatto
)
using ChebyshevApprox
using OffsetArrays

osv(v) = OffsetArray(v, OffsetArrays.Origin(0))
osm(m) = OffsetArray(m, OffsetArrays.Origin(0,0))

function mul(A::AbstractArray, b::AbstractArray)
   nd = length(b)
   N = nd - 1
   c = OffsetArray(zeros(nd), OffsetArrays.Origin(0))
   for j in 0:N, k in 0:N
      c[j] += A[j,k]*b[k]
   end
   return c
end

function test_weights()
   x = osv([1.0, 2.0, 3.0])
   @test barycentric_weights(x) ≈ osv([0.5, -1.0, 0.5])
end

function test_interpolation()
   xg = osv([1.0, 2.0])
   f_values = osv([2.0, 4.0])
   wg = barycentric_weights(xg)
   @test lagrange_interpolation(1.5, xg, f_values, wg) ≈ 3.0
end

function test_diff1()
   xg = osv([1.0, 2.0, 3.0])
   xg_ = osv([1.1, 2.1, 3.1])
   f_vals = 2.0*xg
   wg = barycentric_weights(xg)
   p(x) = lagrange_derivative(x, xg, f_vals, wg)
   @test p.(xg) ≈ osv([2.0,2.0,2.0])
   @test p.(xg_) ≈ osv([2.0,2.0,2.0])
end

function test_diff2()
   xg = osv([1.0, 2.0, 3.0, 4.0])
   xg_ = osv([1.3, 2.3, 3.3])
   f_vals = xg.^3
   wg = barycentric_weights(xg)
   p(x) = lagrange_derivative(x, xg, f_vals, wg)
   @test p.(xg) ≈ 3.0*xg.^2
   @test p.(xg_) ≈ 3.0*xg_.^2
end

function test_diff_mat1()
   xg = osv([1.0, 2.0, 3.0])
   D = differentiation_matrix(1, xg)
   f_vals = 2.0*xg
   @test mul(D, f_vals) ≈ osv([2.0,2.0,2.0])
end

function test_diff_mat2()
   xg = osv([1.0, 2.0, 3.0, 4.0])
   D = differentiation_matrix(1, xg)
   f_vals = xg.^3
   @test mul(D, f_vals) ≈ 3.0*xg.^2
end

function test_diff_mat3()
   xg = osv([1.0, 2.0, 3.0, 4.0])
   D = differentiation_matrix(2, xg)
   f_vals = xg.^3
   @test mul(D,f_vals) ≈ 6.0*xg
end

function test_diff_mat4()
   xg = osv([1.0, 2.0, 3.0, 4.0])
   D = differentiation_matrix(3, xg)
   f_vals = xg.^3
   nd = length(xg)
   @test mul(D,f_vals) ≈ 6.0*osv(ones(nd))
end

function test_chebyshev()
   @test chebyshev_lobatto(degree = 3)[1] ≈ osv(chebyshev_extrema(4))
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
   test_chebyshev()
end
