module bspline
export N, knot_parameters, global_basis_number, INC, IEN, element_number

# The following function evaluates a b-spline function of order 'p' and local index 'i' corresponding
# to the knot vector Ξ. The function is evaluated at the point ξ. Clearly, the arguments 'i', 'p' and 'ξ'
# must be compatible with the knot vector.

function N(i::Int64, p::Int64, Ξ::Array, ξ::Float64)
  if (p < 0)
    print("p < 0")
    return 0
  else
    if p == 0
      if(ξ >= Ξ[i] && ξ< Ξ[i+1])
        return 1
      else
        return 0
      end

    else
      first_term_numerator = (ξ - Ξ[i])*N(i,p-1,Ξ,ξ)
      first_term_denominator = (Ξ[i+p] - Ξ[i])
      second_term_numerator = (Ξ[i+p+1] - ξ)*N(i+1,p-1,Ξ,ξ)
      second_term_denominator = (Ξ[i+p+1] - Ξ[i+1])

      if (first_term_numerator == 0. && first_term_denominator == 0.)
        first_term = 0.
      elseif (first_term_denominator == 0.)
        print("divide by zero error")
        first_term = first_term_numerator/first_term_denominator
      else
        first_term = first_term_numerator/first_term_denominator
      end

      if (second_term_numerator == 0. && second_term_denominator == 0.)
        second_term = 0.
      elseif (second_term_denominator == 0.)
        print("divide by zero error")
        second_term = second_term_numerator/second_term_denominator
      else
        second_term = second_term_numerator/second_term_denominator
      end

      return (first_term + second_term)
    end
  end
end

# The following function returns the primary parameters associated with a knot.
# The function assumes that the first knot is an open knot, and therefore has multiplicity
# of (p+1). The function returns a tuple of integers [n,p]. Here, p is highest order basis
# that can be obtained from the given knot vector. 'n' is the corresponding number of basis
# vectors of order 'p'.

function knot_parameters(Ξ::Array)
  p = 0
  for i in 1:(length(Ξ) - 1)
    if Ξ[i] < Ξ[i+1]
      p = i - 1
      break
    end
  end

  n = length(Ξ) - 1 - p
  return n,p
end

# The following function returns the number of basis functions of order 'p' that are associated with
# the knot vector Ξ.

function knot_parameters(Ξ::Array, p::Int64)
  n1,p1 = knot_parameters(Ξ)
  if p > p1
    print("Order p must be compatible with knot vector Ξ!")
    n = NaN
  else
    n = n1 + p1 - p
  end
  return n
end

# The following function returns the global basis number corresponding to a tensor product NURBS
# construction from the knot vectors Ξ and H. By the ordering (Ξ,H) it is implied that Ξ corresponds to
# "direction 1" and H corresponds to "direction 2". The pair (i,j) corresponds to the NURBS co-ordinate
# of the vertex. For instance, the pair (3,2) corresponds to the vertex of the intersection of ξ₃ and η₂.

function global_basis_number ( Ξ::Array, H::Array, i::Int64, j::Int64, order::Int64)
  n1 = knot_parameters(Ξ, order)
  n2 = knot_parameters(H, order)
  if i > n1 || j > n2 || i < 1 || j < 1
    print("Invalid NURBS co-ordinate! Co-ordinate is out of bounds")
    return NaN
  else
    A = n1*(j-1)+i
  end
  return A
end

# The following function returns the INC (INurbsCoordinate) array corresponding to the knot vectors
# Ξ and H for a given order of bspline basis. The INC array is such that, given a global basis function number
# (refer to the function global_basis_number) and a co-ordinate direction (i.e. 1 or 2), the INC array returns
# the local basis numbering in the specified direction

function INC ( Ξ::Array, H::Array, order::Int64 )
  n1 = knot_parameters(Ξ, order)
  n2 = knot_parameters(H, order)
  A_max = global_basis_number( Ξ, H, n1, n2, order)
  INC_array = zeros(2,A_max)
  for j = 1:n2
    for i = 1:n1
      INC_array[1,n1*(j-1)+i] = i
      INC_array[2,n1*(j-1)+i] = j
    end
  end
  return INC_array
end

# The following function returns the element number of the element that is located in the parametric
# domain Ωᵉ = [ξ(i),ξ(i+1)] x [η(j),η(j+1)]

function element_number( Ξ::Array, H::Array, i::Int64, j::Int64, order::Int64 )
  n1 = knot_parameters(Ξ, order)
  n2 = knot_parameters(H, order)
  if i < order + 1 || i > n1 || j < order + 1 || j > n2
    print("Invalid NURBS co-ordinate for element. Ensure p+1 <= i <= n ")
    return NaN
  else
    e = (j - order - 1)*(n1 - order) + (i - order)
    return e
  end
end

# The following function returns the IEN array corresponding to the knot vectors
# Ξ and H for basis functions of prescribed order. Given an element number 'e' and local basis
# function number 'b', the entry IEN[b,e] is the global number of this basis function.

function IEN ( Ξ::Array, H::Array, order::Int64 )
  n1 = knot_parameters(Ξ, order)
  n2 = knot_parameters(H, order)
  A,elmt = 1,0
  num_of_local_basis_functions = (order+1)*(order+1)
  num_of_elements = element_number(Ξ,H,n1,n2,order)
  IEN_array = zeros(num_of_elements,num_of_local_basis_functions)
  for j = 1:n2
    for i = 1:n1
      if i >= order+1 && j >= order+1
        elmt = elmt+1
        for jloc = 0:order
          for iloc = 0:order
            B = A - jloc*n1 - iloc
            b = jloc*(order+1) + iloc + 1
            IEN_array[elmt,b] = B
          end
        end
      end
      A = A+1
    end
  end
  return IEN_array
end

# Code beyond this point is still under testing and development
# Shape function derivative


function dN_dξ ( i::Int64, p::Int64 , Ξ::Array , x1::Float64 )
  first_term_denominator = (Ξ[i+p] - Ξ[i])
  second_term_denominator = (Ξ[i+p+1] - Ξ[i+1])
  if (first_term_denominator == 0.)
    first_term = 0.
  else
    first_term = 1.*p/(first_term_denominator)*spline_basis(i,p,Ξ,x1)
  end
  if (second_term_denominator == 0.)
    second_term = 0.
  else
    second_term = 1.*p/(second_term_denominator)*spline_basis(i+1,p-1,Ξ,x1)
  end
  return (first_term - second_term)
end

function parametric_coordinates ( KV_Xi, ξ̃ , ni )
  ξ = [ ((KV_Xi[i, ni+1] - KV_Xi[i, ni])*ξ̃[i] + (KV_Xi[i, ni+1] + KV_Xi[i, ni]))/2. for i = 1:3 ...]
  return ξ
end

#bspline module ends here
end
#bspline module ends here
