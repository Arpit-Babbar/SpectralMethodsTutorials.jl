using UnPack
using StaticArrays
using FastGaussQuadrature


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
function Nodal2DStorage(N::Real, M::Real; sol_points=gausslobatto)
   Nd, Md = N + 1, M + 1
   ξ, _ = sol_points(Nd)
   η, _ = sol_points(Md)
   Dξ = differentiation_matrix(1, ξ)
   Dη = differentiation_matrix(1, η)
   D2ξ = differentiation_matrix(2, ξ)
   D2η = differentiation_matrix(2, η)
   # TODO - Remove the static type, these are big arrays
   Nodal2DStorage( N, M, SVector{Nd}(ξ), SVector{Md}(η),
                   SMatrix{Nd,Nd}(Dξ), SMatrix{Md,Md}(Dη),
                   SMatrix{Nd,Nd}(D2ξ), SMatrix{Md,Md}(D2η) )
end

