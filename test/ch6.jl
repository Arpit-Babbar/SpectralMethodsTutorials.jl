using Test
using SpectralMethodsTutorials
(
using SpectralMethodsTutorials: construct_transfinite_quad_M1,
                                test_object_with_data,
                                Nodal2DStorage,
                                MappedGeometry
)

using FastGaussQuadrature

function test_mapped_geometry(test = true; )
   N, M, sol_points, curves_resolution_degree = 10, 10, gausslobatto, 8
   nodal_storage = Nodal2DStorage(N, M, sol_points = sol_points)
   quad_map = construct_transfinite_quad_M1(resolution = curves_resolution_degree)
   geometry = MappedGeometry(nodal_storage, quad_map)
   # return geometry
   @test test_object_with_data(test, geometry, "mg_m1.json")
end

function test_ch6()
   test_mapped_geometry()
end
