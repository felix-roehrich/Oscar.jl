module QuantumGroups
using ..Oscar

import Base: deepcopy_internal

using AbstractAlgebra.PrettyPrinting
using AbstractAlgebra.Generic: FracFieldElem, LaurentPolyWrap

import ..Oscar: add!, addmul!, div!, isone, iszero, mul!, neg!, one, one!, sub!, submul!, zero, zero!
import ..Oscar:
  coefficient_ring, gen, gens, leading_coefficient, leading_exponent_vector, leading_monomial, ngens, parent, trailing_coefficient

include("exports.jl")
include("QuantumGroup.jl")
include("CanonicalBasis.jl")
end

using .QuantumGroups
include("exports.jl")
