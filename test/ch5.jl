using Test
using SpectralMethodsTutorials: solve_potential

function test_solve_potential()
   test_result = true
   # zero test
   p, max_error= solve_potential(5, 5;
                                 source = (x,y) -> zero(eltype(x)),
                                 mask = (true, true, true, true)
                                 )
   @test max_error < 1e-14

   # one test
   p, max_error = solve_potential(5, 5;
                                  source = (x,y) -> zero(eltype(x)),
                                  boundary_condition = (x,y) -> one(eltype(x)),
                                  mask = (false, false, false, false)
                                 )
   @test max_error < 1e-14

   exact_solution = (x,y) -> cospi(2.0*x) * sinpi(2.0*y)
   source = (x,y) -> -8.0 * pi^2 * exact_solution(x,y)

   p, max_error = solve_potential(10, 10,
                                  boundary_condition = exact_solution,
                                  source = source,
                                  mask = (false, false, false, false))
   expected_error = 0.0014297320973608585
   @test abs(max_error - expected_error) < 1e-12

   return test_result
end

function test_ch5()
   test_solve_potential()
end


