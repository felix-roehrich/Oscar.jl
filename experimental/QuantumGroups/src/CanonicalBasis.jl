@doc raw"""
    canonical_basis_elem(U::QuantumGroup, b::Vector{Int}) -> QuantumGroupElem
"""
function canonical_basis_elem(U::QuantumGroup, b::Vector{Int})
  return QuantumGroupElem(U, _canonical_basis_elem(U, b))
end

#function canonical_basis_elem(U::QuantumGroup, b::LusztigDatum)
#  d = deepcopy(b.datum)
#  for mv in braid_moves(weyl_group(U), U.w0, b._w0)
#    Oscar.LieAlgebras._move!(d, mv)
#  end
#  return QuantumGroupElem(U, _canonical_basis_elem(U, d))
#end

function _canonical_basis_elem(U::QuantumGroup, b::Vector{Int})
  bar = bar_involution(U)
  return get!(U.canonical_basis, b) do
    F = one(U)
    for i in 1:length(b)
      for _ in 1:b[i]
        F = mul!(F, gen(U, i))
      end
      mul!(F, inv(quantum_factorial(b[i], U.qi[i])))
    end

    G = F
    F = sub!(bar(F), F)
    while !iszero(F)
      elem = canonical_basis_elem(U, Singular.trailing_exponent_vector(F.elem))
      cf1 = numerator(div(trailing_coefficient(F), trailing_coefficient(elem)))

      # split the coefficient into positive and negative powers of q
      # use postive powers for the canoncial basis element and discard the negative powers
      # this corresponds to b and bar(b)
      cf2 = coefficient_ring(U)(
        parent(cf1.poly)([coeff(cf1, n) for n in 0:degree(cf1.poly)])
      )

      G = addmul!(G, elem, cf2)
      F = submul!(F, elem, coefficient_ring(U)(cf1))
    end

    return G.elem
  end
end
