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
  selem = Singular.reduce(x.elem.sdata * u.elem.elem.sdata, parent(u).I.gb)
  return QuantumGroupModuleElem(
    parent(u), QuantumGroupElem(parent(x), PBWAlgElem(parent(x).alg, selem))
  )
end

function simple_module(U::QuantumGroup, hw::Vector{<:Integer})
  return simple_module(U, WeightLatticeElem(root_system(U), hw))
end

function simple_module(U::QuantumGroup, hw::WeightLatticeElem)
  I = left_ideal([
    chevalley_gen(U, i).elem^(Int(hw[i]) + 1) for i in 1:rank(root_system(U))
  ])
  Oscar.groebner_assure!(I)
  return QuantumGroupModule(U, I)
end

function highest_weight_vector(M::QuantumGroupModule)
  return QuantumGroupModuleElem(M, one(quantum_group(M)))
end
