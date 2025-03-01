module QuantumGroups
using AbstractAlgebra.Generic: LaurentPolyWrap
import ..Oscar: add!, addmul!, div!, mul!, neg!, sub!, submul!, one
import ..Oscar:
  coefficient_ring, leading_coefficient, leading_exponent_vector, leading_monomial, trailing_coefficient
  
  include("exports.jl")
  include("QuantumGroup.jl")
end

include("exports.jl")
