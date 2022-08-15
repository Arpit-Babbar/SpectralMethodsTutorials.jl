using UnPack

struct QuadMap
   X1::Tuple{Float64,Float64}
   X2::Tuple{Float64,Float64}
   X3::Tuple{Float64,Float64}
   X4::Tuple{Float64,Float64}
end

function (quad_map::QuadMap)(ξ,η)
   @unpack X1, X2, X3, X4 = quad_map
   x,y = (
            X1[i]*(1.0-ξ)*(1.0-η) # map (-1,-1) to X1
          + X2[i]*(ξ+1.0)*(1.0-η) # map (+1,-1) to X2
          + X3[i]*(ξ+1.0)*(1.0+η) # map (+1,+1) to X3
          + X4[i]*(1.0-ξ)*(η+1.0) # map (-1,+1) to X4
            for i in 1:2
         )
   return 0.25*x,0.25*y
end

function (quad_map::QuadMap)(Ξ)
   quad_map(Ξ[1],Ξ[2])
end

struct CurveInterpolant # TODO - Communicate array size to compiler?
   N::Float64
   xg::Vector{Float64} # single parameter values at which curve is known
   coords_x::Vector{Float64} # Known x-locations of curve where parameter values are xg
   coords_y::Vector{Float64} # Known y-locations of curve where paramater values are xg
   wg::Vector{Float64}  # Barycentric weights
end

function CurveInterpolant(N, xg, coords)
   wg = barycentric_weights(xg)
   coords_x, coords_y = coords
   return CurveInterpolant(N, xg, coords_x, coords_y, wg)
end

function CurveInterpolant(N, xg, coords_x, coords_y)
   wg = barycentric_weights(xg)
   return CurveInterpolant(N, xg, coords_x, coords_y, wg)
end

function eval_at(curve::CurveInterpolant, s)
   @unpack xg, coords_x, coords_y, wg = curve
   x = lagrange_interpolation(s, xg, coords_x, wg)
   y = lagrange_interpolation(s, xg, coords_y, wg)
   return x,y
end

function diff_at(curve::CurveInterpolant, s)
   @unpack xg, coords_x, coords_y, wg = curve
   x = lagrange_derivative(s, xg, coords_x, wg)
   y = lagrange_derivative(s, xg, coords_y, wg)
   return x, y
end

struct TransfiniteQuadMap{G1<:CurveInterpolant, G2<:CurveInterpolant,
                          G3<:CurveInterpolant, G4<:CurveInterpolant}
   g1::G1
   g2::G2
   g3::G3
   g4::G4
end

# Constructor when curves are given as functions with nodes
function TransfiniteQuadMap(nodes::AbstractVector,
                            g1::Function, g2::Function,
                            g3::Function, g4::Function)
   resolution = length(nodes)
   g1 = CurveInterpolant(resolution, nodes, g1(nodes))
   g2 = CurveInterpolant(resolution, nodes, g2(nodes))
   g3 = CurveInterpolant(resolution, nodes, g3(nodes))
   g4 = CurveInterpolant(resolution, nodes, g4(nodes))
   return TransfiniteQuadMap(g1, g2, g3, g4)
end

function (transfinite_quad_map::TransfiniteQuadMap)(ξ,η)
   @unpack g1, g2, g3, g4 = transfinite_quad_map

   x1, y1 = eval_at(g1, -1.0)
   x2, y2 = eval_at(g1,  1.0)
   x3, y3 = eval_at(g3,  1.0)
   x4, y4 = eval_at(g3, -1.0)

   X1, Y1 = eval_at(g1, ξ)
   X2, Y2 = eval_at(g2, η)
   X3, Y3 = eval_at(g3, ξ)
   X4, Y4 = eval_at(g4, η)
   x  = 0.5*( (1.0-ξ)*X4 + (1.0+ξ)*X2 + (1.0-η)*X1 + (1.0+η)*X3 )
   x += -0.25 * (  (1.0-ξ)*( (1.0-η)*x1 + (1.0+η)*x4 )
                 + (1.0+ξ)*( (1.0-η)*x2 + (1.0+η)*x3 ))

   y  = 0.5*( (1.0-ξ)*Y4 + (1.0+ξ)*Y2 + (1.0-η)*Y1 + (1.0+η)*Y3 )
   y += -0.25 * (  (1.0-ξ)*( (1.0-η)*y1 + (1.0+η)*y4 )
                 + (1.0+ξ)*( (1.0-η)*y2 + (1.0+η)*y3 ))

   return x, y
end

function (transfinite_quad_map::TransfiniteQuadMap)(Ξ)
   return transfinite_quad_map(Ξ[1],Ξ[2])
end

function construct_transfinite_quad_M1(;resolution = 8)
   # We create polynomial approximation of the boundary curves at nodes
   nodes, _ = chebyshev_lobatto(degree = resolution)
   g1(t) = (2.0 .+ t, zero(t))
   g2(t) = (3.0 * cos.( 0.25 * pi * (t .+ 1.0) ),
            3.0 * sin.( 0.25 * pi * (t .+ 1.0) ))
   g3(t) = (zero(t), 2.0 .+ t)
   g4(t) = (cos.(0.25 * pi * (t .+ 1.0)),
            sin.(0.25 * pi * (t .+ 1.0)))
   return TransfiniteQuadMap(nodes, g1, g2, g3, g4)
end

# Computation of the metric terms on a curve-bounded quadrilateral
function transfinite_quad_metrics(quad_map::TransfiniteQuadMap, ξ, η)
   @unpack g1, g2, g3, g4 = quad_map

   x1, y1 = eval_at(g1, -1.0)
   x2, y2 = eval_at(g1,  1.0)
   x3, y3 = eval_at(g3,  1.0)
   x4, y4 = eval_at(g3, -1.0)

   X1, Y1 = eval_at(g1, ξ)
   X2, Y2 = eval_at(g2, η)
   X3, Y3 = eval_at(g3, ξ)
   X4, Y4 = eval_at(g4, η)

   X1p, Y1p = diff_at(g1, ξ)
   X2p, Y2p = diff_at(g2, η)
   X3p, Y3p = diff_at(g3, ξ)
   X4p, Y4p = diff_at(g4, η)

   Xξ  =  0.5  * ( X2 - X4 + (1.0-η)*X1p + (1.0+η)*X3p )
   Xξ += -0.25 * ( (1.0-η)*(x2-x1) + (1.0+η)*(x3-x4) )

   Yξ  =  0.5  * ( Y2 - Y4 + (1.0-η)*Y1p + (1.0+η)*Y3p )
   Yξ += -0.25 * ( (1.0-η)*(y2-y1) + (1.0+η)*(y3-y4) )

   Xη  =  0.5 * ( (1.0-ξ)*X4p + (1.0+ξ)*X2p + X3 - X1 )
   Xη += -0.25 * ( (1.0-ξ)*(x4-x1) + (1.0+ξ)*(x3-x2) )

   Yη  =  0.5 * ( (1.0-ξ)*Y4p + (1.0+ξ)*Y2p + Y3 - Y1 )
   Yη += -0.25 * ( (1.0-ξ)*(y4-y1) + (1.0+ξ)*(y3-y2) )

   return Xξ, Xη, Yξ, Yη
end

# Computation of metric terms on a straight sided quadrilateral
function quad_map_metric(quad_map::QuadMap, ξ, η)
   @unpack X1, X2, X3, X4 = quad_map

   Xξ = 0.25 * ( (1.0-η)*(X2-X1) + (1.0+η)*(X3-X4) )
   Yξ = 0.25 * ( (1.0-η)*(Y2-Y1) + (1.0+η)*(Y3-Y4) )
   Xη = 0.25 * ( (1.0-ξ)*(X4-X1) + (1.0+ξ)*(X3-X2) )
   Yη = 0.25 * ( (1.0-ξ)*(Y4-Y1) + (1.0+ξ)*(Y3-Y2) )

   return Xξ, Xη, Yξ, Yη
end

struct MappedGeometry{RealType}
   # Total points along x,y directions
   N::Int64
   M::Int64

   # Node locations
   x::Array{RealType, 2}
   y::Array{RealType, 2}

   # Boundary node locations
   xb::Array{RealType, 2} # second index to mark curve from 1:4
   yb::Array{RealType, 2} # second index to mark curve from 1:4

   # Metric terms
   Xξ::Array{RealType, 2}
   Xη::Array{RealType, 2}
   Yξ::Array{RealType, 2}
   Yη::Array{RealType, 2}

   J::Array{RealType, 2} # Jacobian

   normal::Array{Tuple{RealType, RealType}, 2} # Boundary normals for each of 4 sides

   scal::Array{RealType, 2} # Scaling factor for each of 4 sides
end

# Constructor
function MappedGeometry(spA::Nodal2DStorage,
                        transfinite_quad_map::TransfiniteQuadMap)
   @unpack N, M, ξ, η, Dξ, Dη = spA
   Nd, Md = N+1, M+1 # Number of nodes in each direction
   RealType = eltype(ξ)
   transfinite_metrics(ξ_, η_) = transfinite_quad_metrics(transfinite_quad_map, ξ_, η_)

   arr = () -> Array{RealType}(undef, Nd, Md)
   x, y, Xξ, Xη, Yξ, Yη, J = ( arr() for _ in 1:7)

   arr = () -> Array{RealType}(undef, max(Nd, Md), 4)
   xb, yb, scal = ( arr() for _ in 1:3 )

   normal = Array{Tuple{RealType,RealType}}(undef , max(Nd,Md), 4)

   # Inner nodes, metrics, jacobian
   for j in 1:Md, i in 1:Nd
      x[i,j], y[i,j] = transfinite_quad_map(ξ[i], η[j])
      Xξ_, Xη_, Yξ_, Yη_ = transfinite_metrics(ξ[i], η[j])
      J[i,j] = Xξ_*Yη_ - Xη_*Yξ_
      Xξ[i,j], Xη[i,j], Yξ[i,j], Yη[i,j] = Xξ_, Xη_, Yξ_, Yη_
   end

   # Boundaries
   # ____3____
   # |       |
   # 4       2
   # |       |
   # ____1____

   # Boundary mapped by vertical edges of reference cell
   for j in 1:Md
      # Right vertical edge
      xb[j,2], yb[j,2] = transfinite_quad_map(1.0, η[j])
      Xξ_, Xη_, Yξ_, Yη_ = transfinite_metrics(1.0, η[j])
      J_ = Xξ_*Yη_ - Xη_*Yξ_
      scal[j,2] = sqrt(Yη_^2 + Xη_^2)
      normal[j,2] = (sign(J_)*  Yη_ /scal[j,2],
                     sign(J_)*(-Xη_)/scal[j,2])

      # Left vertical edge
      xb[j,4], yb[j,4] = transfinite_quad_map(-1.0, η[j])
      Xξ_, Xη_, Yξ_, Yη_ = transfinite_metrics(-1.0, η[j])
      J_ = Xξ_*Yη_ - Xη_*Yξ_
      scal[j,4] = sqrt(Yη_^2 + Xη_^2)
      normal[j,4] = ( -sign(J_)*  Yη_ /scal[j,4],
                      -sign(J_)*(-Xη_)/scal[j,4])
   end

   # Boundary mapped by horizontal edges of reference cell
   for i in 1:Nd
      # Bottom horizontal edge
      xb[i,1], yb[i,1] = transfinite_quad_map(ξ[i], -1.0)
      Xξ_, Xη_, Yξ_, Yη_ = transfinite_metrics(ξ[i], -1.0)
      J_ = Xξ_*Yη_ - Xη_*Yξ_
      scal[i,1] = sqrt(Yξ_^2 + Xξ_^2)
      normal[i,1] = ( -sign(J_)*(-Yξ_) / scal[i,1] ,
                      -sign(J_)*  Xξ_  / scal[i,1] )

      # Top horizontal edge
      xb[i,3], yb[i,3] = transfinite_quad_map(ξ[i], 1.0)
      Xξ_, Xη_, Yξ_, Yη_ = transfinite_metrics(ξ[i], 1.0)
      J_ = Xξ_*Yη_ - Xη_*Yξ_
      scal[i,3] = sqrt(Yξ_^2 + Xξ_^2)
      normal[i,3] = ( sign(J_)*(-Yξ_)/scal[i,3],
                      sign(J_)*  Xξ_ /scal[i,3] )
   end

   MappedGeometry(N, M, x, y, xb, yb, Xξ, Xη, Yξ, Yη, J, normal, scal)

end


