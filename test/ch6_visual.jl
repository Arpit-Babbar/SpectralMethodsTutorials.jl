using Plots
using UnPack
using SpectralMethodsTutorials
using OffsetArrays
(
using SpectralMethodsTutorials: Nodal2DStorage,
                                construct_transfinite_quad_M1,
                                MappedGeometry
)
using FastGaussQuadrature

osv(v) = OffsetArray(v, OffsetArrays.Origin(0))
osm(m) = OffsetArray(m, OffsetArrays.Origin(0,0))

function test_quad_map()
   quad_map_id = QuadMap( (0.,0. ), (1.,0.), (1.,1.), (0.,1.)) # map to unit square
   quad_map = QuadMap( (-0.7,-1.3 ), (1.5,-2.), (0.0,2.4), (-1.,1.))
   lin     = LinRange(-1.0,1.0,100);
   domain = [ (a,b) for a in lin, b in lin ]
   image_square = quad_map_id.(domain)
   image = quad_map.(domain)
   image_square_x = [x[1] for x in image_square]
   image_square_y = [x[2] for x in image_square]
   p1 = scatter(image_square_x, image_square_y, color = "black", label = false)
   image_x = [x[1] for x in image]
   image_y = [x[2] for x in image]
   p2 = scatter(image_x, image_y, color = "black", label = false)
   return p1, p2
end

function draw_quad!(p, quad_map, xl, xr, yl, yr)
   draw_quad!(p, quad_map( (xl, yl) )
               , quad_map( (xr, yl) )
               , quad_map( (xr, yr) )
               , quad_map( (xl, yr) ))
end

function draw_quad!(p, X1, X2, X3, X4)
   x1, y1 = X1
   x2, y2 = X2
   x3, y3 = X3
   x4, y4 = X4
   p!(x,y) = plot!(p, x, y, color = :black, label = false)
   p!([x1,x2],[y1,y2])
   p!([x2,x3],[y2,y3])
   p!([x3,x4],[y3,y4])
   p!([x4,x1],[y4,y1])
end



function test_transfinite_quad(ncells; curves_resolution_degree = 8)
   # TODO - Replace lin by Chebyshev points
   quad_map = construct_transfinite_quad_M1(resolution = curves_resolution_degree)
   p = plot()
   for j in 1:ncells, i in 1:ncells
      xl, xr = -1.0 + (i-1)*2.0/ncells, -1.0 + i*2.0/ncells
      yl, yr = -1.0 + (j-1)*2.0/ncells, -1.0 + j*2.0/ncells
      draw_quad!(p, quad_map, xl, xr, yl, yr)
   end
   return p
end

function draw_mapped_geometry_interior!(p, geometry)
   @unpack x, y, N, M = geometry
   for j in 0:M-1, i in 0:M-1
      draw_quad!(p, (x[i,j],y[i,j]), (x[i+1,j],y[i+1,j]),
                    (x[i+1,j+1],y[i+1,j+1]), (x[i,j+1],y[i,j+1]))
   end
end

function draw_mapped_geometry_boundaries(p, geometry)
   @unpack x, y = geometry
   @unpack N, M = geometry
   b_plot(X, Y) = plot!(p, X, Y, color = :red, label = false)
   for j in 0:M-1
      b_plot([x[0,j], x[0,j+1]], [y[0,j], y[0,j+1]])
      b_plot([x[end,j], x[end,j+1]], [y[end,j], y[end,j+1]])
   end
   for i in 0:N-1
      b_plot([x[i,0], x[i+1,0]], [y[i,0], y[i+1,0]])
      b_plot([x[i,end], x[i+1,end]], [y[i,end], y[i+1,end]])
   end
end

function draw_mapped_geometry_normals(p, geometry; frequency = 3)
   @unpack x, y, N, M = geometry
   n_plot(X, normal) = plot!(p, [X[1], X[1]+normal[1]], [X[2],X[2]+normal[2]],
                             color = :darkblue, label = false, arrow = true )
   for j in 0:M
      if j % frequency == 0
         n_plot( [x[0,j], y[0,j]],
               0.5*collect(geometry.normal[j,4]) )
         n_plot( [x[end,j], y[end,j]],
               0.5*collect(geometry.normal[j,2]) )
      end
   end
   for i in 0:N
      if i % frequency == 0
         n_plot( [x[i,0], y[i,0]],
               0.5*collect(geometry.normal[i,1]) )
         n_plot( [x[i,end], y[i,end]],
               0.5*collect(geometry.normal[i,3]) )
      end
   end
end

function test_mapped_geometry_visual(;N = 10, M = 10,
                                      sol_points = gausslobatto,
                                      curves_resolution_degree = 8,
                                      normal_frequency = 3)
   nodal_storage = Nodal2DStorage(N, M, sol_points = sol_points)
   quad_map = construct_transfinite_quad_M1(resolution = curves_resolution_degree)
   geometry = MappedGeometry(nodal_storage, quad_map)

   p = plot()
   draw_mapped_geometry_interior!(p, geometry)
   draw_mapped_geometry_boundaries(p, geometry)
   draw_mapped_geometry_normals(p, geometry, frequency = normal_frequency)

   return p
end
