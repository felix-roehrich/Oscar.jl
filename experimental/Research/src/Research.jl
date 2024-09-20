module Research
using ..Oscar
import Oscar: GAPWrap

export canonical_basis

include("PathVector.jl")
include("CanonicalBasis.jl")

end

using .Research

export canonical_basis
