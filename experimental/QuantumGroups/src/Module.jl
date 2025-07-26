function quantum_group(M::QuantumGroupModule)
  return M.U
end

function Base.show(io::IO, M::QuantumGroupModule)
  io = pretty(io)
  if is_terse(io)
    print(io, LowercaseOff(), "Quantum group module")
  else
    print(io, LowercaseOff(), "Simple module of highest weight  ")
  end
end

function parent(u::QuantumGroupModuleElem)
  return u.parent
end

function expressify(u::QuantumGroupModuleElem; context=nothing)
  return Expr(:call, :*, expressify(u.elem; context), "u0")
end

@enable_all_show_via_expressify QuantumGroupModuleElem

function Base.:*(x::QuantumGroupElem, u::QuantumGroupModuleElem)
  selem = Singular.reduce(
    x.elem.sdata * u.elem.elem.sdata, Oscar.singular_groebner_basis(parent(u).I)
  )
  return QuantumGroupModuleElem(
    parent(u), QuantumGroupElem(parent(x), PBWAlgElem(parent(x).alg, selem))
  )
end

function simple_module(U::QuantumGroup, hw::Vector{<:Integer})
  return simple_module(U, WeightLatticeElem(root_system(U), hw))
end

function simple_module(U::QuantumGroup, hw::WeightLatticeElem)
  R = root_system(U)
  np = number_of_positive_roots(R)
  I = left_ideal(
    [
      [gen(U.alg, U.cvx[i])^Int(dot(positive_root(R, i), hw) + 1) for i in 1:np];
      [gen(U.alg, np + 2 * i) - U.qi[i]^Int(hw[i]) for i in 1:rank(R)];
      [
        gen(U.alg, np + 2 * i - 1) - quantum_integer(Int(hw[i]), U.qi[i]) for i in 1:rank(R)
      ]
      [gen(U.alg, np + 2 * rank(R) + i) for i in 1:np]
    ],
  )
  return QuantumGroupModule(U, I)
end

function highest_weight_vector(M::QuantumGroupModule)
  return QuantumGroupModuleElem(M, one(quantum_group(M)))
end

function dot(x::QuantumGroupModuleElem, y::QuantumGroupModuleElem)
  S = _dot_antihomomorphism(quantum_group(parent(x)))
  r = S(x.elem) * y
  return leading_coefficient(r.elem)
end
