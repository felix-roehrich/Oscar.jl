module QuantumGroups
using ..Oscar

import Base: deepcopy_internal

using AbstractAlgebra.PrettyPrinting
using AbstractAlgebra.Generic: FracFieldElem, LaurentPolyWrap

import ..Oscar: add!, addmul!, div!, isone, iszero, mul!, neg!, one, one!, sub!, submul!, zero, zero!
import ..Oscar:
  coeff, coefficient_ring, expressify, gen, gens, leading_coefficient, leading_exponent_vector, leading_monomial, ngens, parent, set_coeff!, trailing_coefficient

include("exports.jl")
include("Types.jl")

include("CanonicalBasis.jl")
include("Homomorphism.jl")
include("QuantumGroup.jl")

end

using .QuantumGroups
include("exports.jl")
