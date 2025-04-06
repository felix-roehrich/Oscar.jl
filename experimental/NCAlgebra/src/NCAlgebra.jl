module NCAlgebra

using ..Oscar

using Base: deepcopy_internal

import AbstractAlgebra.Generic: pow!, MPolyRing, MPoly, LaurentPolyWrap, s_polynomial

using AbstractAlgebra.PrettyPrinting: Lowercase, LowercaseOff, pretty

import ..Oscar: coeff, coefficient_ring, elem_type, leading_exponent_vector,
  leading_term, length, gen, gens, ngens, parent, parent_type, setcoeff!, tail, zeros

import ..Oscar: left_ideal, groebner_basis
import ..Oscar: add!, mul!, one, one!, sub!, submul!, zero, zero!, is_one, is_zero
import ..Oscar: divexact!, divrem, rem!, rem


include("PBWAlgebra.jl")
include("PBWAlgebraIdeal.jl")
include("Types.jl")
include("PBWAlgebraModule.jl")

export PBWAlg, PBWAlgebraElem
end

using .NCAlgebra

export PBWAlg, PBWAlgElem
