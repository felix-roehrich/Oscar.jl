module LSTheory

using ..Oscar
using AbstractAlgebra.PrettyPrinting

import Oscar: expressify, points

include("Poset.jl")
include("LSFan.jl")

export Poset
export LSLattice, LSLatticeElem
export LSMonoid
export ls_fan_of_monoids
export ls_lattice
export maximal_chains
export poset

end

using .LSTheory

export Poset
export LSLattice, LSLatticeElem
export LSMonoid
export ls_fan_of_monoids
export ls_lattice
export maximal_chains
export poset
