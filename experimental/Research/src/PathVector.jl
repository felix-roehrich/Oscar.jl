struct PathVectorContext
  U::GAP.GapObj
  gens::Vector{GAP.GapObj}
  
  Mv::GAP.GapObj
  
  function PathVectorContext(wt::WeightLatticeElem)
    R = root_system(wt)
    GAP.Globals.LoadPackage(GAP.Obj("QuaGroup"))
    
    # for now only do type A
    U = GAP.Globals.QuantizedUEA(GAP.Globals.RootSystem(GAP.Obj("A"), rank(R)))
    Mv = GAP.Globals.HighestWeightModule(U, GAP.Obj(Int.(Oscar._vec(coefficients(wt)))))
    
    B = GAP.Globals.CanonicalBasis(U)
    gens = [ GAP.Globals.MonomialElements(B, GAP.Obj(Int.(i .== 1:rank(R))))[1] for i in 1:rank(R)]
  
    return new(U, gens, Mv)
  end
end

function PathVectorContext(R::RootSystem, wt::Vector{Int})
  return PathVectorContext(WeightLatticeElem(R, wt))
end

struct PathVector
  ctx::PathVectorContext
  path::LSPathModelElem
end

function (p::PathVector)(i::Vector{Int}, n::Vector{Int})
  iter = PathVectorSummandIterator(p, i, n)

  c = 0
  s = zero(iter.field)
  for m in iter
    c += 1
    s += coefficient(iter, m)
  end

  return c, s
end

struct PathVectorSummandIterator
  l::Int
  i::Vector{Int}
  n::Vector{Int}
  p::PathVector
  
  na::Vector{Vector{Int}} # number of alpha
  
  # field
  field::AbsSimpleNumField
  
  # caching
  act::Dict{Vector{Int}, GAP.GapObj}
  bound::Dict{Tuple{Int, Vector{Int}}, Int}
  coeff::Dict{Vector{Int}, AbsSimpleNumFieldElem}
end

function PathVectorSummandIterator(p::PathVector, i::Vector{Int}, n::Vector{Int})
  l = lcm(2, map(s -> Int(denominator(s.t)), p.path.s)...)
  n = [l*n[i] for i in length(n):-1:1]
  
  na = sizehint!(Vector{Int}[], l)
  for i in 1:length(p.path.s)
    num = Int(l * p.path.s[i].t)
    
    nai = zeros(Int, rank(root_system(parent(p.path).wt)))
    wt = deepcopy(parent(p.path).wt)
    for s in Iterators.reverse(word(p.path.s[i].w))
      si = Int(s)
      nai[si] += Int(wt[si])
      reflect!(wt, si)
    end
    
    while num > 0
      push!(na, nai)
      num -= 1
    end
  end
  
  F, _ = cyclotomic_field(2*l)
  return PathVectorSummandIterator(l, reverse(i), n, p, na, F, Dict([0] => GAPWrap.Basis(p.ctx.Mv)[1]), Dict(), Dict())
end

function _action(iter::PathVectorSummandIterator, n::Vector{Int})
  v = get(iter.act, n, nothing)
  if !isnothing(v)
    return v
  end
  
  l = length(n)
  n2 = copy(n)
  while isnothing(v)
    if n2[l] > 0
      n2[l] -= 1
    else
      l -= 1
    end
    v = get(iter.act, view(n2, 1:l), nothing)
  end
  
  while l <= length(n)
    while n2[l] < n[l]
      n2[l] += 1
      v = iter.p.ctx.gens[iter.i[l]]^v / GAP.Globals.GaussNumber(n2[l], GAP.Globals._q)
      iter.act[n2[1:l]] = v
    end
    l += 1
  end
  
  return v
end

function Base.iterate(iter::PathVectorSummandIterator)
  m = zeros(Int, length(iter.n), iter.l)
  return initialize(iter, m, 1)
end

function Base.iterate(iter::PathVectorSummandIterator, m::Matrix{Int})
  for i in length(iter.i)-2:-1:1
    nr = next(m[i, :], bound(iter, m, i))
    if isnothing(nr)
      continue
    end
    
    m[i, :] = nr
    return initialize(iter, m, i+1)
  end

  return nothing
end

function initialize(iter::PathVectorSummandIterator, m::Matrix{Int}, row::Int)
  i = row
  while i <= length(iter.i)
    m[i, :] = bound(iter, m, i)
    
    s = 0
    for j in 1:iter.l
      s += m[i, j]
      if s >= iter.n[i]
        m[i, j] -= s - iter.n[i]
        m[i, j+1:iter.l] .= 0
        break
      end
    end
    
    if s < iter.n[i]
      i -= 1
      np = Vector{Int}()
      while i > 0
        np = next(m[i, :], bound(iter, m, i))
        if !isnothing(np)
          break
        end
        i -= 1
      end
      if i == 0
        return nothing
      end
      m[i, :] = np
    end
    
    i += 1
  end
  
  return m, m
end

function bound(iter::PathVectorSummandIterator, m::Matrix{Int}, row::Int)
  b = zeros(Int, iter.l)
  c = zeros(Int, row)
  for j in 1:ncols(m)
    copy!(c, @view m[1:row, j])
    b[j] = get!(iter.bound, (j, m[1:row-1, j])) do
      f = iter.na[j][iter.i[row]] - sum(m[i, j] for i in 1:row-1 if iter.i[i] == iter.i[row]; init=0)
      c[row] = 1
      while f > 0
        v = _action(iter, c)
        if GAPWrap.IsZero(v)
          break
        end
        f -= 1
        c[row] += 1
      end
      return c[row]-1
    end
  end
  
  return b
end

function next(p::Vector{Int}, b::Vector{Int})
  np = copy(p)
  len = length(np)

  r = np[len] + 1
  np[len] = 0
  for i in (len - 1):-1:1
    if np[i] == 0
      continue
    end
    
    # overflow over the bound
      of = np[i + 1] + r - b[i + 1]
      if of <= 0
        np[i] -= 1
        np[i + 1] += r
        return np
      elseif of > sum(b[j] for j in (i + 2):len; init=0)
        r += np[i]
        np[i] = 0
        continue
      end

      np[i] -= 1
      np[i + 1] += r
      for j in (i + 2):len
        np[j - 1] -= of
        np[j] += of
        of -= b[j]
        if of <= 0
          return np
        end
      end
  end

  return nothing
end

function coefficient(iter::PathVectorSummandIterator, m::Matrix{Int})
  P = parent(iter.p.path)
  cm = cartan_matrix(root_system(P))
  
  A = iter.field
  q = gen(A)
  coeff = one(A)
  wt = [deepcopy(P.wt) for _ in 1:ncols(m)]
  for i in 1:nrows(m)
    s = iter.n[i]
    for j in ncols(m):-1:2
      s -= m[i, j]
      addmul!(wt[j].vec, cm[:, iter.i[i]:iter.i[i]], -m[i, j])
      coeff *= q^-(s*(m[i, j]+Int(wt[j].vec[iter.i[i]])))
    end
  end
  
  for j in 1:ncols(m)
    if haskey(iter.coeff, m[:, j])
      coeff *= iter.coeff[m[:, j]]
      if iszero(coeff)
        break
      end
      continue
    end

    v = _action(iter, m[:, j])
    if GAPWrap.IsZero(v)
      error("should not happen")
    end
    
    rep = GAP.Globals.ExtRepOfObj(GAP.Globals.ExtRepOfObj(v))
    coeffs = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[2])
    coeff *= get!(iter.coeff, m[:, j], sum(l -> coeffs[1][l] * q^(l - 1 + coeffs[2]), 1:length(coeffs[1])))
  end

  return q
end

function canonical_basis(R::RootSystem, deg::Vector{Int}, pts::Vector{LSPathModelElem})
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
