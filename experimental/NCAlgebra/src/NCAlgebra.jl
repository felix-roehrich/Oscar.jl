module NCAlgebra

using ..Oscar

using Base: deepcopy_internal

import AbstractAlgebra.Generic: pow!, MPolyRing, MPoly, LaurentPolyWrap

using AbstractAlgebra.PrettyPrinting: Lowercase, LowercaseOff, pretty

import ..Oscar: coeff, coefficient_ring, leading_exponent_vector,
  leading_term, gen, gens, ngens, tail, zeros

import ..Oscar: add!, mul!, one, one!, zero, zero!, is_one, is_zero

include("PBWAlgebra.jl")

export PBWAlg, PBWAlgebraElem
end

using .NCAlgebra

export PBWAlg, PBWAlgElem
