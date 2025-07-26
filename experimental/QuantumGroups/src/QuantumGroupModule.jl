function simple_module(U::QuantumGroup, wt::Vector{Int})
  return simple_module(U, WeightLatticeElem(root_system(U), wt))
end

function simple_module(U::QuantumGroup, wt::WeightLatticeElem)
  @req is_dominant(wt) "weight must be dominant"

  R = root_system(U)
  npos = number_of_positive_roots(R)
  I = pbw_algebra_ideal(
    [
      zero_pbw_gen(U, 2 * i) - coefficient_ring(U)(wt[i]) for i in 1:rank(R);
      zero_pbw_gen(U, 2 * i - 1) - q_integer(Int(wt[i]), U.qi) for i in 1:rank(R);
      positive_pbw_gen(U, i) for i in 1:npos
    ],
  )

  return QuantumGroupModule(U, I)
end
