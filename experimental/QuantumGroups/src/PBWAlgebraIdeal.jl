function pbw_algebra_ideal(gens::Vector{PBWAlgebraElem{T}}) where {T}
  return PBWAlgebraIdeal(gens, PBWAlgebraElem{T}[])
end

function groebner_basis(I::PBWAlgebraIdeal)
  if !isempty(I.gb)
    return I.gb
  end

  G = deepcopy(I.gens)
  r = zero(I.gens[1])
  ok = false
  while !ok
    ok = true
    for i in 1:length(G), j in (i + 1):length(G)
      s = s_polynomial(G[i], G[j])
      rem!(r, s, G)
      if !iszero(r)
        ok = false
        push!(G, divexact!(deepcopy(r), coeff(r, 1)))
      end
    end
  end

  #I.gb = G
  return G
end

function s_polynomial(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  A = parent(x)

  z = zeros(Int, ngens(A))
  monomial_lcm!(
    z, AbstractAlgebra.exponent_vector_ref(x.poly, 1),
    AbstractAlgebra.exponent_vector_ref(y.poly, 1),
  )

  px = zero(x.poly)
  _mul_m_p!(A, px, z - AbstractAlgebra.exponent_vector_ref(x.poly, 1), x.poly)

  py = zero(y.poly)
  _mul_m_p!(A, py, z - AbstractAlgebra.exponent_vector_ref(y.poly, 1), y.poly)

  c = deepcopy(coeff(px, 1))
  d = gcd(coeff(px, 1), coeff(py, 1))

  px = AbstractAlgebra.divexact!(mul!(px, coeff(py, 1)), d)
  py = AbstractAlgebra.divexact!(mul!(py, c), d)

  return PBWAlgebraElem(A, sub!(px, py))
end
