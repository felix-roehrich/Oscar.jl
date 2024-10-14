##### Kashiwara Crystal #####

abstract type AbstractCrystal <: AbstractAlgebra.Set end

abstract type AbstractCrystalElem <: AbstractAlgebra.SetElem end

# AbstractCrystalElem required methods

function Base.eps(b::AbstractCrystalElem, i::Int)
  error("not implemented")
end

function phi(b::AbstractCrystalElem, i::Int)
  error("not implemented")
end

@doc raw"""
    Apply the root operator $\tilde e_i$ to `b` in place and return the result.
"""
function ealpha!(b::AbstractCrystalElem, i::Int)
  error("not implemented")
end

@doc raw"""
    Apply the root operator $\tilde f_i$ to `b` in place and return the result.
"""
function falpha!(b::AbstractCrystalElem, i::Int)
  error("not implemented")
end

@doc raw"""
    Return the weight of `b`.
"""
function weight(b::AbstractCrystalElem)
  error("not implemented")
end

# AbstractCrystal optional methods

@doc raw"""
    is_upper_normal(B::AbstractCrystal) -> bool

Return `true` if `B` is upper normal, and `false` if `B` is not upper normal or it is unknown.
In the case of `false` consult the specific implemenation for details.
"""
function is_upper_normal(B::AbstractCrystal)
  return true
end

@doc raw"""
    is_lower_normal(B::AbstractCrystal) -> bool

Return `true` if `B` is lower normal, and `false` if `B` is not lower normal or it is unknown.
In the case of `false` consult the specific implemenation for details.
"""
function is_lower_normal(B::AbstractCrystal)
  return true
end

@doc raw"""
    is_normal(B::AbstractCrystal) -> bool

Return `true` if `B` is normal, and `false` if `B` is not normal or it is unknown,
whether `B` is upper normal or lower normal.
In the case of `false` consult the specific implemenation for details.
"""
function is_normal(B::AbstractCrystal)
  return is_upper_normal(B) && is_lower_normal(B)
end

# AbstractCrystalElem provided methods

# ealpha

@doc raw"""
    ealpha(b::AbstractCrystalElem, i::Int) -> AbstractCrystalElem

Return the result of applying the root operator $\tilde e_i$ to `b`.
"""
function ealpha(b::AbstractCrystalElem, i::Int)
  return ealpha!(deepcopy(b), i)
end

@doc raw"""
    ealpha(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem

Return the result of applying the root operator $\tilde e_i$ to `b` `n`-times.
"""
function ealpha(b::AbstractCrystalElem, i::Int, n::Int)
  return ealpha!(deepcopy(b), i, n)
end

@doc raw"""
    ealpha(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int}) -> AbstractCrystalElem

Return the result of applying the root operator $\tilde e_i$ to `b`.
"""
function ealpha(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  return ealpha!(deepcopy(b), i, n)
end

@doc raw"""
    ealpha!(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem
"""
function ealpha!(b::AbstractCrystalElem, i::Int, n::Int)
  for _ in 1:n
    ealpha!(b, i)
  end
  return b
end

@doc raw"""
    ealpha!(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int}) -> AbstractCrystalElem
"""
function ealpha!(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  @req length(i) == length(n) "i and n must have same length"
  
  for l in length(i):-1:1
    ealpha!(b, i[l], n[l])
  end
  return b
end

# falpha

@doc raw"""
    falpha(b::AbstractCrystalElem, i::Int) -> AbstractCrystalElem

Return the result of applying the root operator $\tilde f_i$ to `b`.
"""
function falpha(b::AbstractCrystalElem, i::Int)
  return falpha!(deepcopy(b), i)
end

@doc raw"""
    falpha(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem

Return the result of applying the root operator $\tilde f_i$ to `b` `n`-times.
"""
function falpha(b::AbstractCrystalElem, i::Int, n::Int)
  return falpha!(deepcopy(b), i, n)
end

@doc raw"""
    falpha(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int}) -> AbstractCrystalElem

Return the result of applying the root operator $\tilde f_i$ to `b`.
"""
function falpha(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  return falpha!(deepcopy(b), i, n)
end

function falpha!(b::AbstractCrystalElem, i::Int, n::Int)
  for _ in 1:n
    falpha!(b, i)
  end
  return b
end

function falpha!(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  @req length(i) == length(n) "i and n must have same length"
  
  for l in length(i):-1:1
    falpha!(b, i[l], n[l])
  end
  return b
end

# other

function adapted_string(p::AbstractCrystalElem, rdec::Vector{Int})
  s = zero(rdec)
  b = deepcopy(p)
  for i in 1:length(rdec)
    s[i] = Base.eps(b, rdec[i]) # TODO: this will not work for non normal crystals
    ealpha!(b, rdec[i], s[i])
  end
  return s
end

# LS Path Model 

struct LSPathModel <: AbstractCrystal
  wt::WeightLatticeElem
  ext::Dict{Vector{UInt8}, WeightLatticeElem} # map from Weyl group elements to extremal weights
  
  function LSPathModel(wt::WeightLatticeElem)
    @req is_dominant(wt) "weight must be dominant"
    return new(wt, Dict(UInt8[] => wt))
  end
end

# Young Tableu Path Model

struct YTPathModel <: AbstractCrystal
  wt::WeightLatticeElem
end

struct YTPathModelElem <: AbstractCrystalElem
  parent::YTPathModel
  t::YoungTableau
end

function Base.show(io::IO, mime::MIME"text/plain", P::LSPathModel)
  #@show_name(io, P)
  #@show_special(io, mime, P)
  #io = pretty(io)
  print(io, "LS path model for ")
  show(io, mime, P.wt)
end

function Base.show(io::IO, P::LSPathModel)
  #@show_name(io, R)
  #@show_special(io, R)
  if is_terse(io)
    print(io, "LS path model")
  else
    print(io, "LS path model for $(P.wt)")
  end
end

@doc raw"""
    ls_path_model(wt::WeightLatticeElem) -> LSPathModel

Return the LS path model for the dominant weight `wt`.
"""
function ls_path_model(wt::WeightLatticeElem)
  return LSPathModel(wt)
end

@doc raw"""
    ls_path_model(R::RootSystem, wt::Vector{<:IntegerUnion}) -> LSPathModel

Return the LS path model for the root system `R` and dominant weight `wt`.
"""
function ls_path_model(R::RootSystem, wt::Vector{<:IntegerUnion})
  return ls_path_model(WeightLatticeElem(R, wt))
end

function _extremal_weight(P::LSPathModel, w::WeylGroupElem)
  wt = get(P.ext, word(w), nothing)
  if isnothing(wt)
    # we need to make a copy, because w may be modified
    return get!(P.ext, deepcopy(word(w)), w*P.wt)
  end
  
  return wt
end

function (P::LSPathModel)(wv::Vector{WeylGroupElem}, tv::Vector{<:RationalUnion})
  return LSPathModelElem(P, [LSPathSegment(t, w) for (t, w) in Iterators.zip(tv, wv)])
end

function (P::LSPathModel)(p::YTPathModelElem)
  @req P.wt == parent(p).wt "dominant weights must match"
  
  pp = dominant_path(P)
  i = Int.(word(longest_element(weyl_group(P.wt))))
  falpha!(pp, i, adapted_string(p, i))
  return pp
end

function root_system(P::LSPathModel)
  return root_system(P.wt)
end

# LSPathModelElem

struct LSPathSegment
  t::QQFieldElem
  w::WeylGroupElem
end

function Base.:(==)(s1::LSPathSegment, s2::LSPathSegment)
  return s1.t == s2.t && s1.w == s2.w
end

struct LSPathModelElem <: AbstractCrystalElem
  parent::LSPathModel
  s::Vector{LSPathSegment}
end

function parent(p::LSPathModelElem)
  return p.parent
end

function Base.:(==)(p::LSPathModelElem, q::LSPathModelElem)
  return parent(p) === parent(q) && p.s == q.s
end

function Base.:(<)(p::LSPathModelElem, q::LSPathModelElem)
  for i in 1:min(length(p.s), length(q.s))
    if p.s[i].w < q.s[i].w
      return true
    elseif p.s[i].w == q.s[i].w
      if p.s[i].t < q.s[i].t
        return true
      elseif p.s[i].t > q.s[i].t
        return false
      end
    else
      return false
    end
  end
  return false
end

function Base.deepcopy_internal(p::LSPathModelElem, dict::IdDict)
  if haskey(dict, p)
    return dict[p]
  end

  p2 = LSPathModelElem(parent(p), deepcopy_internal(p.s, dict))
  dict[p] = p2
  return p2
end

function Base.hash(p::LSPathModelElem, h::UInt)
  b = 0xd379e91be1f36479 % UInt
  h = hash(parent(p), h)
  h = hash(p.d, h)
  h = hash(p.s, h)

  return xor(h, b)
end

function expressify(p::LSPathModelElem, s=:e; context=nothing)
  sum = Expr(:call, :+)
  for seg in p.s
    push!(
      sum.args,
      Expr(:call, :*, expressify(seg.t; context=context), "$s($(seg.w))"),
    )
  end
  return sum
end
@enable_all_show_via_expressify LSPathModelElem

@doc raw"""
    max(p::LSPathModelElem) -> WeylGroupElem

  Return the maximal element in the support of `p`.
"""
function Base.max(p::LSPathModelElem)
  return p.s[1].w
end

@doc raw"""
    ls_sequence(p::LSPathModelElem) -> Tuple{Vector{WeylGroupElem}, Vector{QQFieldElem}}
    
Return the LS sequence for `p`.
"""
function ls_sequence(p::LSPathModelElem)
  return map(s -> s.w, p.s), [zero(QQ); accumulate((t, s) -> t + s.t, p.s; init=zero(QQ))]
end

@doc raw"""
    dominant_path(P::LSPathModel) -> LSPathModelElem
    
Return the dominant path for `P`.
"""
function dominant_path(P::LSPathModel)
  return LSPathModelElem(P, [LSPathSegment(one(QQ), one(weyl_group(root_system(P))))])
end

@doc raw"""
    halpha(p::LSPathModelElem, i::Int) -> Vector{QQFieldElem}

Return the rational turning points of the function $h_\alpha$ for `p` where $\alpha$ is the `i`th simple root.
"""
function halpha(p::LSPathModelElem, i::Int)
  h = sizehint!([zero(QQ)], length(p.s)+1)
  for k in 1:length(p.s)
    push!(h, h[k] + p.s[k].t * _extremal_weight(parent(p), p.s[k].w)[i])
  end
  return h
end

# implement AbstractCrystalElem methods

function ealpha!(p::LSPathModelElem, i::Int)
  h = halpha(p, i)
  
  j = argmin(h) # first time the global min is assumed
  k = findprev(>=(h[j]+1), h, j-1)
  
  if isnothing(k)
    empty!(p.s)
    return p
  end
  
  of = (h[j]+1-h[k])//_extremal_weight(parent(p), p.s[k].w)[i]
  # of > 0, so we need to split this segment
  if !iszero(of)
    insert!(p.s, k, LSPathSegment(of, deepcopy(p.s[k].w)))
    sub!(p.s[k+1].t, p.s[k+1].t, of)
    k += 1
    j += 1
  end
  if j <= length(p.s) && length(p.s[j-1].w) == length(p.s[j].w)+1
    add!(p.s[j].t, p.s[j].t, p.s[j-1].t)
    deleteat!(p.s, j-1)
    j -= 1
  end
  
  for l in k:j-1
    lmul!(p.s[l].w, i)
  end

  return p
end

function falpha!(p::LSPathModelElem, i::Int)
  h = halpha(p, i)
  
  # find the last time the global minimum is assumed
  j = 1
  for l in 2:length(h)
    if h[l] <= h[j]
      j = l
    end
  end
  
  k = findnext(>=(h[j]+1), h, j+1)
  if isnothing(k)
    empty!(p.s)
    return p
  end
  
  of = (h[k]-h[j]-1)//_extremal_weight(parent(p), p.s[k-1].w)[i]
  # of > 0, we need to cut the segment
  if !iszero(of) 
    insert!(p.s, k, LSPathSegment(of, deepcopy(p.s[k-1].w)))
    sub!(p.s[k-1].t, p.s[k-1].t, of)
  end
  if j > 1 && length(p.s[j-1].w) == length(p.s[j].w)+1
    add!(p.s[j-1].t, p.s[j-1].t, p.s[j].t)
    deleteat!(p.s, j)
    k -= 1
  end
  
  for l in j:k-1
    lmul!(p.s[l].w, i)
  end

  return p
end

function Base.eps(p::LSPathModelElem, i::Int)
  return -Int(minimum(halpha(p, i)))
end

function phi(p::LSPathModelElem, i::Int)
  h = halpha(p, i)
  return Int(h[end]) - Int(minimum(h))
end

function weight(p::LSPathModelElem)
  v = sum(s -> mul!(QQ.(coefficients(_extremal_weight(parent(p), s.w))), s.t), p.s)
  return WeightLatticeElem(root_system(parent(p)), ZZ.(v))
end

# LSPathModel specific methods

function iszero(p::LSPathModelElem)
  return isempty(p.d)
end

function (P::LSPathModel)(v::Vector{<:IntegerUnion})
  return P(WeightLatticeElem(root_system(P), v))
end

function (P::LSPathModel)(w::WeightLatticeElem)
  nf = inv(QQ.(cartan_matrix(root_system(P)))) * coefficients(P.wt - w)
  if !all(is_integral, nf)
    return LSPathModelElem[]
  end
  nf = [Int(nf[i]) for i in 1:length(nf)]

  W = weyl_group(root_system(P))
  w0 = longest_element(W)
  path = zeros(Int, length(w0))
  a = dominant_path(P)

  points = LSPathModelElem[]
  num = zeros(Int, length(nf)) # number of times acted with fi
  
  G = gens(W)
  
  i = length(w0)
  ret = false
  while true
    s = Int(w0[i])
    
    mi = min(phi(a, s), nf[s]-num[s])
    falpha!(a, s, mi)
    num[s] += mi
    path[i] += mi
    
    ok = num == nf
    if i == 1 || ok
      if ok && !(a in points)
        push!(points, deepcopy(a))
      end
      
      ealpha!(a, s, path[i])
      num[s] -= path[i]
      path[i] = 0
      
      while path[i] == 0
        i += 1
        if i > length(w0)
          @goto done
        end
      end
      
      s = Int(w0[i])
      ealpha!(a, s)
      num[s] -= 1
      path[i] -= 1
    end
    i -= 1
  end
  @label done

  return points
end

# Young Tableau path model implemenation

function (P::YTPathModel)(p::LSPathModelElem)
  @req P.wt == parent(p).wt "dominant weights must match"
  
  pp = dominant_path(P)
  i = Int.(word(max(p)))
  falpha!(pp, i, adapted_string(p, i))
  return pp
end

function yt_path_model(wt::WeightLatticeElem)
  return YTPathModel(wt)
end

function yt_path_model(R::RootSystem, wt::Vector{Int})
  return YTPathModel(WeightLatticeElem(R, wt))
end

function Base.show(io::IO, mime::MIME"text/plain", p::YTPathModelElem)
  show(io, mime, p.t)
end

function dominant_path(P::YTPathModel)
  rk= rank(root_system(P.wt))
  shape = zeros(Int, rk)
  shape[rk] = Int(P.wt[rk])
  for i in rk-1:-1:1
    shape[i] = shape[i+1] + Int(P.wt[i])
  end
  
  return YTPathModelElem(P, young_tableau([fill(i, shape[i]) for i in 1:length(shape)]))
end

function ealpha!(p::YTPathModelElem, i::Int)
  t = p.t.t
  
  i = (1, length(t[r]))
  m = 0
  h = 0
  for r in 1:length(t)
    for c in length(t[r]):-1:1
      if t[r][c] == i
        h += 1
      elseif t[r][c] == i+1
        h -= 1
        if h <= m
          m = h
          i = (r, c)
        end
      end
    end
  end
  
  
end

function falpha!(p::YTPathModelElem, i::Int)
  t = p.t.t
  
  # find last time the global minimum is assumed
  k, l = 1, length(t[1])+1
  m = 0
  h = 0
  for r in 1:length(t)
    for c in length(t[r]):-1:1
      if t[r][c] == i
        h += 1
      elseif t[r][c] == i+1
        h -= 1
      end
      
      if h <= m
        m = h
        k, l = r, c
      end
    end
  end
  
  if m < h
    if l == 1
      k += 1
      l = length(t[k])+1
    end
    t[k][l-1] += 1
  end
  
  return p
end

function global_eps(p::LSPathModelElem, i::Integer)
  n = QQFieldElem(0)
  for s in p.s
    b, _, _ = explain_lmul(s.w, i)
    if !b
      n -= s.t*_extremal_weight(parent(p), s.w)[i]
    else
      n = floor(n)
    end
  end
  return floor(Int, n)
end

function total_eps(p::LSPathModelElem, i::Integer)
  n = QQFieldElem(0)
  for s in p.s
    b, _, _ = explain_lmul(s.w, i)
    if !b
      n -= s.t*_extremal_weight(parent(p), s.w)[i]
    end
  end
  return n
end

function is_balanced(p::LSPathModelElem)
  for w in reduced_expressions(longest_element(weyl_group(root_system(parent(p)))))
    t = deepcopy(p)
    for i in w
      ii = Int(i)
      n = eps(t, ii)
      if global_eps(t, ii) != n
        return false
      end
      ealpha!(t, ii, n)
    end
  end
  return true
end

function is_w0_compatible(p::LSPathModelElem)
  rk = rank(root_system(parent(p)))
  q = deepcopy(p)
  for i in length(q.s):-1:1
    while length(q.s[i].w) > 0
      s = findfirst(!first(explain_lmul(q.s[i].w, s)) for s in 1:rk)
      for j in 1:i
        l = length(q.s[j].w)
        if length(lmul!(q.s[j].w, s)) > l
          return false
        end
      end
    end
  end
  
  return true
end

function rev_triangle_greater(p::LSPathModelElem, q::LSPathModelElem)
  for i in 0:min(length(p.s), length(q.s))-1
    if p.s[end-i].w < q.s[end-i].w
      return true
    elseif p.s[end-i].w == q.s[end-i].w
      if p.s[end-i].t > q.s[end-i].t
        return true
      elseif p.s[end-i].t < q.s[end-i].t
        return false
      end
    else
      return false
    end
  end
  return false
end

function initial_string(p::LSPathModelElem, w::Vector{Int})
  str = zeros(Int, length(w))
  n = zero(QQ)
  
  t = deepcopy(p)
  for i in 1:length(w)
    for s in t.s
      b, j, _ = explain_lmul(s.w, w[i])
      if !b
        deleteat!(word(s.w), j)
        n += s.t*_extremal_weight(parent(p), s.w)[w[i]]
      else
        break
      end
    end
    
    str[i] = Int(n)
    zero!(n)
  end

  return str
end

function test_idea(t::Symbol, r::Int, wt::Vector{Int})
  GAP.Packages.load("QuaGroup")
  
  R = root_system(t, r)
  W = weyl_group(R)
  P = ls_path_model(R, wt)
  L = lie_algebra(QQ, t, r)

  char = character(L, wt)
  w0 = longest_element(weyl_group(R))
  for k in keys(char)
    pts = P(k)
    mx = filter(p -> !any(>(p), pts), pts)
    for p in mx
      dem = filter(q -> q != p && max(q) == max(p), pts)
      for q in dem
        if p.s[1].t != q.s[1].t && all(eps(q, i) >= -p.s[1].t*_extremal_weight(P, p.s[1].w)[i] for i in 1:r)
          return p, q
        end
      end
    end
  end
end

#=
if !is_w0_compatible(p) # || !isone(p.s[end].w)
        continue
      end
      
      rd = Int[]
      q = deepcopy(p)
      for i in length(q.s):-1:1
        while !isone(q.s[i].w)
          push!(rd, q.s[i].w[1])
          for k in 1:i
            lmul!(q.s[k].w, q.s[i].w[1])
          end
        end
      end
      append!(rd, word(inv(W(rd; normalize=false))*w0))
      
      s = adapted_string(p, rd)
      deg = [sum(s[j] for j in 1:length(s) if rd[j] == i) for i in 1:r]
      
      gR = GAP.Globals.RootSystem(GAP.Obj(t), r)
      GAP.Globals.SetLongestWeylWord(gR, GAP.Obj(rd))
      gU = GAP.Globals.QuantizedUEA(gR)
      gB = GAP.Globals.CanonicalBasis(gU)  
      strs = Vector{Int}.(GAP.Globals.Strings(gB, GAP.Obj(deg)))
      
      gs = Int[]
      for j in 1:length(rd)
        if s[j] == 0
          continue
        end
        push!(gs, Int(rd[j]), s[j])
      end
      ind = findfirst(==(gs), strs)
      if length(GAP.Globals.ExtRepOfObj(GAP.Globals.MonomialElements(gB, GAP.Obj(deg))[ind])) > 2
        return p
      end
      =#

struct TensorProductCrystal <: AbstractCrystal
  B1::AbstractCrystal
  B2::AbstractCrystal
end

struct TensorProductCrystalElem <: AbstractCrystalElem
  parent::TensorProductCrystal

  b1::AbstractCrystalElem
  b2::AbstractCrystalElem
end

function expressify(b::TensorProductCrystalElem; context=nothing)
  t = Expr(:call, :*)
  push!(t.args, expressify(b.b1; context))
  push!(t.args, expressify(b.b2; context))
  return t
end
@enable_all_show_via_expressify TensorProductCrystalElem

function Base.eps(b::TensorProductCrystalElem, i::Int)
  return max(eps(b.b1, i), eps(b.b2, i) - weight(b.b1)[i])
end

function phi(b::TensorProductCrystalElem, i::Int)
  return max(phi(b.b2, i), eps(b.b1, i) + weight(b.b2)[i])
end

function falpha!(b::TensorProductCrystalElem, i::Int)
  if phi(b.b1, i) > eps(b.b2, i)
    falpha!(b.b1, i)
  else
    falpha!(b.b2, i)
  end
  return b
end

function ealpha!(b::TensorProductCrystalElem, i::Int)
  if phi(b.b1, i) >= eps(b.b2, i)
    ealpha!(b.b1, i)
  else
    ealpha!(b.b2, i)
  end
  return b
end


