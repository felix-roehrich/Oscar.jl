module QuantumGroups
using ..Oscar

import Base: deepcopy_internal

using AbstractAlgebra.PrettyPrinting
using AbstractAlgebra.Generic: LaurentPolyWrap

import ..Oscar: add!, addmul!, div!, mul!, neg!, sub!, submul!, one, one!, zero, zero!
import ..Oscar:
  coefficient_ring, gen, gens, leading_coefficient, leading_exponent_vector, leading_monomial, trailing_coefficient
  
  include("exports.jl")
  include("QuantumGroup.jl")
end

using .QuantumGroups
include("exports.jl")
