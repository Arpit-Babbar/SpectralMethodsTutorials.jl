using UnPack
using StaticArrays
using FastGaussQuadrature
using OffsetArrays
using LinearSolve
using Plots

# Vector and Matrix type are parametrized so user can choose static
# if size is small
struct Nodal2DStorage{RealType<:Real, AbstractVectorNd, AbstractVectorMd,
                      AbstractMatrixTypeNd, AbstractMatrixTypeMd}
   N::RealType
   M::RealType
   ξ::AbstractVectorNd # N+1 Gauss-Legendre-Lobatto (GLL) points in [-1,1]
   η::AbstractVectorMd # M+1 Gauss-Legendre-Lobatto (GLL) points in [-1,1]
   Dξ::AbstractMatrixTypeNd # N+1 point GLL differentiation matrix
   Dη::AbstractMatrixTypeMd # M+1 point GLL differentiation matrix
   Dξ2::AbstractMatrixTypeNd # N+1 point GLL 2nd derivative differentiation matrix
   Dη2::AbstractMatrixTypeMd # M+1 point GLL 2nd derivative differentiation matrix
end

# Constructor
function Nodal2DStorage(N::Int, M::Int; sol_points=gausslobatto)
   Nd, Md = N + 1, M + 1
   ξ, _ = sol_points(Nd)
   η, _ = sol_points(Md)
   osv(v) = OffsetArray(v, OffsetArrays.Origin(0))
   ξ, η = osv(SVector{Nd}(ξ)), osv(SVector{Md}(η))
   Dξ = differentiation_matrix(1, ξ)
   Dη = differentiation_matrix(1, η)
   D2ξ = differentiation_matrix(2, ξ)
   D2η = differentiation_matrix(2, η)
   # TODO - Remove the static type, these are big arrays
   Nodal2DStorage( N, M, ξ, η,
                   Dξ  , Dη, D2ξ , D2η )
end

# TODO - Rename to NodalPotentialSquare?
struct NodalPotential{RealType, BoundaryCondition}
   nodal_storage::Nodal2DStorage
   Φ::OffsetArray{RealType, 2, Array{RealType, 2}}
   s::OffsetArray{RealType, 2, Array{RealType, 2}}
   mask::Tuple{Bool, Bool, Bool, Bool}
   boundary_condition::BoundaryCondition
end

function NodalPotential(N::Int, M::Int;
                        source::Function = (x,y) -> zero(eltype(x)),
                        mask = (true, true, true, true),
                        boundary_condition::Function = (x,y) -> zero(eltype(x)))
   nodal_storage = Nodal2DStorage(N, M)
   @unpack ξ, η = nodal_storage
   osm(m) = OffsetArray(m, OffsetArrays.Origin(0,0))
   Φ, s = osm(zeros(N+1, M+1)), osm(zeros(N+1, M+1))

   # fill source values
   for j in 0:M, i in 0:N
      s[i,j] = source(ξ[i], η[j])
   end

   return NodalPotential(nodal_storage, Φ, s, mask, boundary_condition)
end

function mask_sides!(nodal_potential::NodalPotential)
   @unpack mask, Φ, nodal_storage, boundary_condition  = nodal_potential
   @unpack ξ, η = nodal_storage
   reset(u) = fill!(u, zero(eltype(u)))
   if mask[1] == true
      @views reset(Φ[:,0])
   else
      bc_bottom = x -> boundary_condition(x,η[0])
      Φ[:,0] .= bc_bottom.(ξ)
   end
   if mask[2] == true
      @views reset(Φ[end,:])
   else
      bc_right = y -> boundary_condition(ξ[end],y)
      Φ[end,:] .= bc_right.(η)
   end
   if mask[3] == true
      @views reset(Φ[:,end])
   else
      bc_top = x -> boundary_condition(x, η[end])
      Φ[:,end] .= bc_top.(ξ)
   end
   if mask[4] == true
      @views reset(Φ[0,:])
   else
      bc_left = y -> boundary_condition(ξ[0], y)
      Φ[0,:] .= bc_left.(η)
   end
end

# TODO - Isn't this assuming Dirichlet bc?
function collocate_rhs(nodal_potential::NodalPotential)
   @unpack s, Φ, nodal_storage = nodal_potential
   @unpack N, M, Dξ2, Dη2 = nodal_storage
   L = (N-1)*(M-1)
   rhs = zeros(L)
   for j in 1:M-1, i in 1:N-1
      n = flatten_indices(i,j,N)
      rhs[n] = (s[i,j] - Dξ2[i,0]*Φ[0,j] - Dξ2[i,N]*Φ[N,j]
                       - Dη2[j,0]*Φ[i,0] - Dη2[j,M]*Φ[i,M])
   end
   return rhs
end

function collocate_laplace_matrix(nodal_potential::NodalPotential)
   @unpack Φ, nodal_storage = nodal_potential
   @unpack N, M, Dξ2, Dη2 = nodal_storage
   L = (N-1)*(M-1)
   A = zeros(L, L)

   for j in 1:M-1, i in 1:N-1
      n = flatten_indices(i,j,N)
      for k in 1:N-1
         m = flatten_indices(k,j,N)
         A[n,m] = Dξ2[i,k]
      end

      for k in 1:M-1
         m = flatten_indices(i,k,N)
         A[n,m] += Dη2[j,k]
      end
   end

   return A
end

function solve_potential(N, M;
                         source::Function = (x,y) -> zero(eltype(x)),
                         boundary_condition::Function = (x,y) -> zero(eltype(x)),
                         mask = (true, true, true, true))
   nodal_potential = NodalPotential(N, M, source = source, mask = mask,
                                    boundary_condition = boundary_condition)

   # Construct the matrix problem
   mask_sides!(nodal_potential)
   rhs = collocate_rhs(nodal_potential)
   A = collocate_laplace_matrix(nodal_potential)

   # solve the matrix problem
   prob = LinearProblem(A, rhs)
   linsolve = init(prob)
   sol = LinearSolve.solve(linsolve)
   u = reshape(sol.u, (N-1, M-1))

   # Compute error
   @unpack nodal_storage = nodal_potential
   @unpack N, M, ξ, η = nodal_storage
   max_error = 0.0
   for j in 1:M-1, i in 1:N-1
      max_error = max(abs(boundary_condition(ξ[i], η[j]) -  u[i,j] ), max_error)
   end

   p = @views surface(ξ[1:N-1], η[1:M-1], u, xlabel = "x", ylabel = "y")

   return p, max_error
end

# Functions to be added for iterative solver
# LaplcianOnTheSquare, MatrixAction
