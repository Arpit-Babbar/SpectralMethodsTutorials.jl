using StaticArrays
using OffsetArrays

osv(v) = OffsetArray(v, OffsetArrays.Origin(0))
osm(m) = OffsetArray(m, OffsetArrays.Origin(0,0))

function barycentric_weights(x)
   nd = length(x)
   N = nd - 1
   w = osv(ones(eltype(x), nd))
   for j in 1:N
      for k in 0:j-1
         w[k] *= x[k] - x[j]
         w[j] *= x[j] - x[k]
      end
   end
   for j in 0:N
      w[j] = 1.0/w[j]
   end
   return w
end

function lagrange_interpolation(x, xg, f_vals, wg)
   N = length(xg) - 1
   numerator = denominator = zero(eltype(x))
   for j in 0:N
      if x ≈ xg[j]
         return f_vals[j]
      end
      t = wg[j]/(x-xg[j])
      numerator += f_vals[j] * t
      denominator += t
   end
   return numerator/denominator
end

function lagrange_derivative(x, xg, f_vals, wg)
   N = length(xg) - 1
   at_node = false
   numerator = zero(eltype(x))
   local i,p,denominator
   for j in 0:N
      if x ≈ xg[j]
         at_node = true
         i = j
         p = f_vals[j]
         denominator = -wg[j]
         break
      end
   end
   if at_node
      for j in 0:N
         if j != i
            numerator += wg[j] * (p-f_vals[j])/(x-xg[j])
         end
      end
   else
      denominator = zero(eltype(x))
      p = lagrange_interpolation(x, xg, f_vals, wg)
      for j in 0:N
         t = wg[j]/(x-xg[j])
         numerator += t * (p-f_vals[j])/(x-xg[j])
         denominator += t
      end
   end
   return numerator/denominator
end

function differentiation_matrix(xg::AbstractVector)
   nd = length(xg)
   N = nd - 1
   D = osm(zeros(eltype(xg), nd, nd))
   w = barycentric_weights(xg)
   for j in 0:N, i in 0:N
      if j != i
         D[i,j]  = w[j]/ (w[i] * (xg[i] - xg[j]))
         D[i,i] -= D[i,j]
      end
   end
   Dm = OffsetArray(SMatrix{nd,nd}(D), OffsetArrays.Origin(0))
   return Dm
end

function differentiation_matrix(m::Int, x::AbstractVector)
   nd = length(x)
   N = nd - 1
   w = barycentric_weights(x)
   D = osm(zeros(eltype(x), nd, nd))
   D .= differentiation_matrix(x)
   if m == 1
      return D
   end
   for k in 2:m
      for i in 0:N
         Dii = D[i,i]
         D[i,i] = zero(eltype(D))
         for j in 0:N
            if j != i
               D[i,j]  = k/(x[i]-x[j]) * ( w[j]/w[i] * Dii - D[i,j] )
               D[i,i] -= D[i,j]
            end
         end
      end
   end
   return osm(SMatrix{nd,nd}(D))
end

function chebyshev_lobatto(;degree::Int)
   nd = degree + 1
   x, w = osv(zeros(nd)), osv(zeros(nd))
   for j in 0:degree
      x[j] = -cospi(j/degree)
      w[j] = pi/degree
   end
   w[1]   *= 0.5
   w[end] *= 0.5

   return osv(SVector{nd}(x)), osv(SVector{nd}(w))
end
