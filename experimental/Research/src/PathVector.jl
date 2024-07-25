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

struct PathVectorSummandIterator
  l::Int
  i::Vector{Int}
  n::Vector{Int}
  p::PathVector
  
  na::Vector{Vector{Int}} # number of alpha
  
  # field
  field::AbsSimpleNumField
  
  # caching
  bound::Dict{Vector{Int}, Int}
  coeff::Dict{Vector{Int}, AbsSimpleNumFieldElem}
end

function PathVectorSummandIterator(p::PathVector, i::Vector{Int}, n::Vector{Int})
  l = lcm(2, map(d -> Int(denominator(d)), p.path.d)...)
  n = [l*n[i] for i in length(n):-1:1]
  
  na = sizehint!(Vector{Int}[], l)
  for i in 1:length(p.path.d)
    num = Int(l * p.path.d[i])
    
    nai = zeros(Int, rank(root_system(parent(p.path).wt)))
    wt = deepcopy(parent(p.path).wt)
    for s in Iterators.reverse(word(p.path.s[i]))
      si = Int(s)
      nai[si] += Int(wt[si])
      reflect!(wt, si)
    end
    
    while num > 0
      push!(na, nai)
      num -= 1
    end
  end
  
  F, _ = cyclotomic_field(2*iter.l)
  return PathVectorSummandIterator(l, reverse(i), n, p, na, F, Dict(), Dict())
end

function Base.iterate(iter::PathVectorSummandIterator)
  m = zeros(Int, length(iter.n), iter.l)
  return initialize(iter, m, 1)
end

function Base.iterate(iter::PathVectorSummandIterator, m::Matrix{Int})
  for i in length(iter.i)-2:-1:1
    nr = next(m[i, :], bound(iter, m, i))
    if isnothing(nr)
      i -= 1
      continue
    end
    
    m[i, :] = nr
    return initialize(m)
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
  for j in 1:ncols(m)
    fa = copy(iter.na[j])
    v = GAPWrap.Basis(iter.p.ctx.Mv)[1]
    for i in 1:row-1
      fa[iter.i[i]] -= m[i, j]
      for _ in 1:m[i, j]
        v = iter.p.ctx.gens[iter.i[i]]^v
      end
    end
    
    while !GAPWrap.IsZero(v) && fa[iter.i[row]] > 0
      v = iter.p.ctx.gens[iter.i[row]]^v
      fa[iter.i[row]] -= 1
      b[j] += 1
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
  cm = cartan_matrix(root_system(iter.p.wt))
  
  A = iter.field
  q = gen(A)
  coeff = one(A)
  for i in 1:nrows(m)
    s = iter.n[i]
    for j in 1:ncols(m)
      s -= m[i, j]
      q *= q^-(s*(m[i, j]+w[j].vec[iter.rdec[i]]))
    end
  end
  
  for j in 1:ncols(mat)
    if haskey(iter.coeff, mat[:, j])
      q *= iter.coeff[mat[:, j]]
      if q == 0
        break
      end
      continue
    end

    v = GAPWrap.Basis(iter.p.ctx.Mv)[1]
    for i in 1:nrows(mat)
      for _ in 1:mat[i, j]
        v = iter.p.ctx.gens[iter.i[i]]^v
      end
      v /= GAP.Globals.GaussianFactorial(mat[i, j], GAP.Globals._q)
    end
    if GAPWrap.IsZero(v)
      q = 0
      iter.coeff[mat[:, i]] = zero(A)
      break
    end
    
    rep = GAP.Globals.ExtRepOfObj(GAP.Globals.ExtRepOfObj(wv))
    coeffs = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[2])
    coeff = sum(l -> coeffs[1][l] * v^(l - 1 + coeffs[2]), 1:length(coeffs[1]))
    iter.cc_cache[mat[:, i]] = coeff
    q *= coeff
  end

  return q
end
