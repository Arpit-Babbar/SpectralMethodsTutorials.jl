using StaticArrays

function barycentric_weights(x)
   nd = length(x)
   w = ones(eltype(x), nd)
   for j in 2:nd
      for k in 1:j-1
         w[k] *= x[k] - x[j]
         w[j] *= x[j] - x[k]
      end
   end
   for j in 1:nd
      w[j] = 1.0/w[j]
   end
   return w
end

function lagrange_interpolation(x, xg, f_vals, wg)
   nd = length(xg)
   numerator = denominator = zero(eltype(x))
   for j in 1:nd
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
   nd = length(xg)
   at_node = false
   numerator = zero(eltype(x))
   local i,p,denominator
   for j in 1:nd
      if x ≈ xg[j]
         at_node = true
         i = j
         p = f_vals[j]
         denominator = -wg[j]
         break
      end
   end
   if at_node
      for j in 1:nd
         if j != i
            numerator += wg[j] * (p-f_vals[j])/(x-xg[j])
         end
      end
   else
      denominator = zero(eltype(x))
      p = lagrange_interpolation(x, xg, f_vals, wg)
      for j in 1:nd
         t = wg[j]/(x-xg[j])
         numerator += t * (p-f_vals[j])/(x-xg[j])
         denominator += t
      end
   end
   return numerator/denominator
end

function differentiation_matrix(xg::AbstractVector)
   nd = length(xg)
   D = zeros(eltype(xg), nd, nd)
   w = barycentric_weights(xg)
   for j in 1:nd, i in 1:nd
      if j != i
         D[i,j]  = w[j]/ (w[i] * (xg[i] - xg[j]))
         D[i,i] -= D[i,j]
      end
   end
   return SMatrix{nd,nd}(D)
end

function differentiation_matrix(m::Int64, x::AbstractVector)
   nd = length(x)
   w = barycentric_weights(x)
   D = zeros(eltype(x), nd, nd)
   D .= differentiation_matrix(x)
   if m == 1
      return SMatrix{nd,nd}(D)
   end
   for k in 2:m
      for i in 1:nd
         Dii = D[i,i]
         D[i,i] = zero(eltype(D))
         for j in 1:nd
            if j != i
               D[i,j]  = k/(x[i]-x[j]) * ( w[j]/w[i] * Dii - D[i,j] )
               D[i,i] -= D[i,j]
            end
         end
      end
   end
   return SMatrix{nd,nd}(D)
end

# TODO - Where do they come from??
function chebyshev_lobatto(;degree)
   nd = degree + 1
   x, w = zeros(nd), zeros(nd)
   for j in 0:degree
      x[j+1] = -cospi(j/degree)
      w[j+1] = pi/degree
   end
   w[1]   *= 0.5
   w[end] *= 0.5
   return x, w
end
