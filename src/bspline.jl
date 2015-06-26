module bspline
export N, knot_parameters

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

function knot_parameters(Ξ::Array)
  p = 0
  for i in 1:(length(Ξ) - 1)
    if Ξ[i] < Ξ[i+1]
      p = i - 1
      break
    end
  end

  n = length(Ξ) - p
  return n,p
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
