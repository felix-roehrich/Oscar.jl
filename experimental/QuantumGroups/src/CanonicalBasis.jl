struct CanonicalBasisElem
  c::Vector{Int}
  f::Vector{Vector{Tuple{Int,Int}}}
  
  function CanonicalBasisElem()
    return new([], [])
  end
end

function (a::LSFanElem)(b::CanonicalBasisElem)
  z = 0
  for i in 1:length(b.f)
    _, c = tensor_coefficient(b.f[i], a)
    z += b.c[i]*c
  end
  return z
end

function canonical_basis(R::RootSystem, deg::Vector{Int}, pts::Vector{LSFanElem})
  gR = GAP.Globals.RootSystem(GAP.Obj("A"), rank(R))
  gU = GAP.Globals.QuantizedUEA(gR)
  gB = GAP.Globals.CanonicalBasis(gU)
  
  basis = CanonicalBasisElem[]
  
  elems = GAP.Globals.MonomialElements(gB, GAP.Obj(deg))
  strs = GAP.Globals.Strings(gB, GAP.Obj(deg))
  
  rdec = [1,2,1,3,2,1] # LongestWeylWord
  for i in 1:length(pts)
    s = adapted_string(pts[i], rdec)
    ss = [[rdec[j], s[j]] for j in 1:length(rdec) if s[j] != 0]
    
    j = 1
    while j < length(strs)
      if collect(Iterators.partition(strs[j], 2)) == ss
        break
      end
      j += 1
    end
    strs[i], strs[j] = strs[j], strs[i]
    elems[i], elems[j] = elems[j], elems[i]
  end
  
  #elems = elems[1:length(pts)]
  for el in elems
    b = CanonicalBasisElem()
    rep = GAP.Globals.ExtRepOfObj(el)
    for (f, c) in Iterators.partition(rep, 2)
      v = Tuple{Int,Int}[]
      p = i -> i == 1 ? 1 : i == 3 ? 2 : 3
      
      for (i, n) in Iterators.partition(f, 2)
        push!(v, (p(i), n))
      end
      push!(b.c, GAP.Globals.Value(c, 1)) # c is a laurent polynomial
      push!(b.f, v)
    end
    push!(basis, b)
  end
  
  return basis
end

function filter_sort_canonical_basis(can::Vector{CanonicalBasisElem}, pts::Vector{LSFanElem})  
  rdec = [1,2,1,3,2,1]
  for i in 1:length(pts)
    s = adapted_string(pts[i], rdec)
    j = findfirst(b -> b.f[1] == [(rdec[j], s[j]) for j in 1:length(rdec) if s[j] != 0], can)
    can[i], can[j] = can[j], can[i]
  end
  return can[1:length(pts)]
end
