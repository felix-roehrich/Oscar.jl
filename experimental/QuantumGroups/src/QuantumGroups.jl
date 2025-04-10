module QuantumGroups
using ..Oscar

import Base: deepcopy_internal

using AbstractAlgebra.PrettyPrinting
using AbstractAlgebra.Generic: FracFieldElem, LaurentPolyWrap

import Base: factorial

import ..Oscar:
  add!, addmul!, div!, isone, iszero, mul!, neg!, one, one!, sub!, submul!, zero, zero!
import ..Oscar:
  coeff, coefficient_ring, expressify, gen, gens, leading_coefficient,
  leading_exponent_vector, leading_monomial, length, ngens, parent, setcoeff!,
  trailing_coefficient
import ..Oscar: weyl_group, simple_module

include("exports.jl")
include("Types.jl")

include("CanonicalBasis.jl")
include("Homomorphism.jl")
include("Module.jl")
include("QuantumGroup.jl")

end

using .QuantumGroups
include("exports.jl")
