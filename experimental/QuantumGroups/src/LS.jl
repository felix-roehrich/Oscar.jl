###############################################################################
#
#   Lakshmibai-Seshadri things
#
###############################################################################

L, v = laurent_polynomial_ring(ZZ, "v")

struct LSFan
  gap_U::GAP.GapObj
  gap_V::GAP.GapObj
  highest_weight::WeightLatticeElem
  d::Vector{ZZRingElem}

  dual_root_system::RootSystem
  dual_heighst_weight::WeightLatticeElem
end

function ls_fan(w::WeightLatticeElem)
  GAP.Globals.LoadPackage(GAP.Obj("QuaGroup"))

  cm = cartan_matrix(root_system(w))
  d = cartan_symmetrizer(cm)

  dual_rs = root_system(transpose(cm))
  dual_hw = WeightLatticeElem(dual_rs, d .* vec(w))

  # gap uses the transposed version for Cartan matrices
  type = cartan_type(cartan_matrix(dual_rs))
  U = GAP.Globals.QuantizedUEA(
    GAP.Globals.RootSystem(GAP.Obj(map(x -> GAP.Obj(x), Iterators.flatten(type))))
  )
  V = GAP.Globals.HighestWeightModule(U, GAP.Obj(Int.(vec(dual_hw))))

  return LSFan(U, V, w, d, dual_rs, dual_hw)
end

function ls_fan(R::RootSystem, v::Vector{<:IntegerUnion})
  return ls_fan(WeightLatticeElem(R, matrix(ZZ, rank(R), 1, v)))
end

function (LS::LSFan)()
  return LSFanElem(LS, (1, one(weyl_group(LS))))
end

# !!! missing check that the constructed point belongs to LS
function (LS::LSFan)(v::Tuple{RationalUnion,WeylGroupElem}...; check=true)
  return LSFanElem(LS, v...; check=check)
end

function (LS::LSFan)(v::Tuple{RationalUnion, Vector{<:Integer}}...; check=true)
  return LSFanElem(LS, map(t -> (t[1], WeylGroupElem(root_system(LS), t[2])), v)...; check=check)
end

function bonds(hw::WeightLatticeElem, x::WeylGroupElem)
  b = zeros(Int, length(x))
  w = deepcopy(hw)
  for i in length(x):-1:1
    s = Int(word(x)[i])
    b[i] += w[s]
    reflect!(w, s)
  end
  return b
end

function highest_weight(LS::LSFan)
  return LS.highest_weight
end

function points_of_weight(LS::LSFan, v::Vector{Int})
  return points_of_weight(LS, WeightLatticeElem(root_system(LS), v))
end

function points_of_weight(LS::LSFan, w::WeightLatticeElem)
  nf = inv(change_base_ring(QQ, cartan_matrix(root_system(LS)))) * coefficients(highest_weight(LS) - w)
  if !all(is_integral, nf)
    return LSFanElem[]
  end
  nf = [Int(nf[i]) for i in 1:length(nf)]

  W = weyl_group(root_system(LS))
  w0 = longest_element(W)
  path = zeros(Int, length(w0))
  a = LS()

  points = LSFanElem[]
  num = zeros(Int, length(nf)) # number of times acted with fi

  for rdec in reduced_expressions(w0; up_to_commutation=true)
    i = length(w0)
    while true
      s = Int(rdec[i])
      while num[s] < nf[s] && !isnothing(fi!(a, s))
        num[s] += 1
        path[i] += 1
      end
      
      if num == nf
        if !(a in points)
          push!(points, deepcopy(a))
        end
        ei!(a, s, path[i])
        num[s] -= path[i]
        path[i] = 0
      elseif i == 1
        ei!(a, s, path[i])
        num[s] -= path[i]
        path[i] = 0
      end

      while path[i] == 0
        i += 1
        if i > length(w0)
          @goto done
        end

        s = Int(rdec[i])
        ei!(a, s)
        num[s] -= 1
        path[i] -= 1
      end
      i -= 1
    end

    @label done
  end

  return points
end

function apply_sequence(LS::LSFan, seq::Vector{Int})
  W = weyl_group(root_system(LS))
  w0 = longest_element(W)

  a = LS()
  for i in length(w0):-1:1
    if isnothing(fi!(a, Int(w0[i]), seq[i]))
      return nothing
    end
  end
  return a
end

#=
function points_of_weight(LS::LSFan, w::WeightLatticeElem)
  nf = inv_cartan_matrix(root_system(LS)) * coefficients(highest_weight(LS) - w)
  if !all(is_integral, nf)
    return LSFanElem[]
  end
  nf = [Int(nf[i]) for i in 1:length(nf)]

  W = weyl_group(root_system(LS))
  w0 = longest_element(W)
  path = Vector{LSFanElem}(undef, length(w0) + 1)
  path[end] = LS((1, one(W)))

  points = LSFanElem[]
  num = zeros(Int, length(nf)) # number of times acted with fi
  nfi = zeros(Int, length(w0)) # number of times acted at i
  for rdec in reduced_expressions(w0; up_to_commutation=true)
    i = length(rdec)
    while i < length(path)
      while i > 0
        s = Int(rdec[i])
        if num[s] == nf[s]
          path[i] = path[i + 1]
          i -= 1
          continue
        end

        path[i] = deepcopy(path[i + 1])
        while num[s] < nf[s] && fi!(path[i], s) != nothing
          num[s] += 1
          nfi[i] += 1
        end
        if nfi[i] > 0
          i -= 1
          continue
        end
        break
      end

      i += 1
      if num == nf && !(path[1] in points)
        push!(points, deepcopy(path[1]))
      end
      while i < 2
        num[Int(rdec[i])] -= nfi[i]
        nfi[i] = 0
        i += 1
      end
      while i < length(path)
        if nfi[i] > 0
          num[Int(rdec[i])] -= 1
          nfi[i] -= 1
          if nfi[i] > 0
            ei!(path[i], Int(rdec[i]))
            i -= 1
            break
          end
        end
        i += 1
      end
    end
    zero!(num)
    zero!(nfi)
  end

  return points
end
=#

function root_system(LS::LSFan)
  return root_system(highest_weight(LS))
end

function weyl_group(LS::LSFan)
  return weyl_group(root_system(LS))
end

###############################################################################
# LSFanElem

# internal type
mutable struct lsPair
  a::QQFieldElem
  w::WeylGroupElem
  weight::WeightLatticeElem
end

function Base.hash(p::lsPair, h::UInt)
  b = 0x6986f85c9414efe1 % UInt
  h = hash(p.a, h)
  h = hash(p.w, h)

  return xor(h, b)
end

function Base.:(==)(a::lsPair, b::lsPair)
  return a.a == b.a && a.w == b.w
end

struct LSFanElem
  parent::LSFan
  vec::Vector{lsPair}
end

function LSFanElem(LS::LSFan, v::Tuple{RationalUnion,WeylGroupElem}...; check::Bool=true)
  @req !check || sum(a for (a, _) in v) == 1 "the sum of all coefficients needs to be 1"
  @req !check || allunique(length(w) for (_, w) in v) "the support must be thin"

  return LSFanElem(LS, [lsPair(a, w, w * highest_weight(LS)) for (a, w) in v])
end

function Base.deepcopy_internal(a::LSFanElem, dict::IdDict)
  if haskey(dict, a)
    return dict[a]
  end

  a2 = LSFanElem(parent(a), deepcopy_internal(a.vec, dict))
  dict[a] = a2
  return a2
end

function Base.getindex(a::LSFanElem, x::WeylGroupElem)
  i = findfirst(p -> p.w == x, a.vec)
  if isnothing(i)
    return 0
  end

  return a.vec[i].a
end

function Base.hash(a::LSFanElem, h::UInt)
  b = 0xd379e91be1f36479 % UInt
  h = hash(parent(a), h)
  h = hash(a.vec, h)

  return xor(h, b)
end

function Base.iszero(a::LSFanElem)
  return isempty(a.vec)
end

function Base.:(<)(a::LSFanElem, b::LSFanElem)
  for i in 1:min(length(a.vec), length(b.vec))
    if a.vec[i].w < b.vec[i].w
      return true
    elseif a.vec[i].w == b.vec[i].w
      if a.vec[i].a < b.vec[i].a
        return true
      elseif a.vec[i].a > b.vec[i].a
        return false
      end
    else
      return false
    end
  end

  return false
end

function Base.:(==)(a::LSFanElem, b::LSFanElem)
  return parent(a) === parent(b) && a.vec == b.vec
end

function Base.copy(a::LSFanElem)
  return LSFanElem(a.parent, deepcopy(a.vec))
end

@doc raw"""

Returns the maximal Weyl group element in the support of `a`.
"""
function Base.max(a::LSFanElem)
  return a.vec[1].w
end

function Base.parent(a::LSFanElem)
  return a.parent
end

function Base.vec(a::LSFanElem)
  return a.vec
end

function expressify(a::LSFanElem, s=:e; context=nothing)
  sum = Expr(:call, :+)
  for i in 1:length(a.vec)
    push!(
      sum.args,
      Expr(:call, :*, expressify(a.vec[i].a; context=context), "$s($(a.vec[i].w))"),
    )
  end
  return sum
end
@enable_all_show_via_expressify LSFanElem

function weight(a::LSFanElem)
  vec = sum(x -> x.a * change_base_ring(QQ, x.weight.vec), a.vec)
  return WeightLatticeElem(root_system(parent(a)), change_base_ring(ZZ, vec))
end

function height(a::LSFanElem, i::Int)
  h = sizehint!([zero(QQ)], length(a.vec))
  for r in 1:length(a.vec)
    push!(h, h[r] + a.vec[r].a*a.vec[r].weight[i])
  end
  return h
end

function eps(a::LSFanElem, i::Int)
  return -Int(minimum(height(a, i)))
end

function fi!(a::LSFanElem, i::Int)
  h = height(a, i)

  mi = minimum(h)
  j = findlast(==(mi), h)
  k = findnext(>=(mi+1), h, j+1)
  if isnothing(k)
    return nothing
  end

  diff = 1//a.vec[k-1].weight[i]*(h[j]+1-h[k-1])
  if a.vec[k-1].a > diff
    k -= 1
    wt = reflect(a.vec[k].weight, i)
    if k == 1 || a.vec[k-1].weight != wt
      insert!(a.vec, k, lsPair(diff, lmul(a.vec[k].w, i), wt))
      a.vec[k+1].a -= diff
    else
      a.vec[k-1].a += diff
      a.vec[k].a -= diff
    end
  end
  if 1 < j < k && length(a.vec[j-1].w) == length(a.vec[j].w)+1
    a.vec[j-1].a += a.vec[j].a
    deleteat!(a.vec, j)
    k -= 1
  end

  for r in j:k-1
    reflect!(a.vec[r].weight, i)
    lmul!(a.vec[r].w, i)
  end

  return a
end

function fi!(a::LSFanElem, i::Int, n::Int)
  for _ in 1:n
    if isnothing(fi!(a, i))
      return nothing
    end
  end
  return a
end

function fi!(a::LSFanElem, i::Vector{<:Integer}, n::Vector{<:Integer})
  @req length(i) === length(n) "i and n must have the same length"

  for l in length(i):-1:1
    fi!(a, i[l], n[l])
  end
  return a
end

function ei!(a::LSFanElem, i::Int)
  h = height(a, i)
  mi = minimum(h)
  j = findfirst(==(mi), h)

  k = findprev(>=(mi+1), h, j-1)
  if isnothing(k)
    return nothing
  end

  diff = 1//a.vec[k].weight[i]*(h[j]+1-h[k])
  if diff > 0 && a.vec[k].a > diff
    wt = reflect(a.vec[k].weight, i)
    if k == length(a.vec) || a.vec[k+1].weight != wt
      insert!(a.vec, k+1, lsPair(0, lmul(a.vec[k].w, i), wt))
    end
    a.vec[k+1].a += a.vec[k].a - diff
    a.vec[k].a = diff
    k += 2
    j += 1
  end
  if k < j <= length(a.vec) && length(a.vec[j-1].w) == length(a.vec[j].w)+1
    a.vec[j].a += a.vec[j-1].a
    deleteat!(a.vec, j-1)
    j -= 1
  end

  for r in k:j-1
    reflect!(a.vec[r].weight, i)
    lmul!(a.vec[r].w, i)
  end

  return a
end

function ei!(a::LSFanElem, i::Int, n::Int)
  for _ in 1:n
    if isnothing(ei!(a, i))
      return nothing
    end
  end
  return a
end

function Base.adjoint(a::LSFanElem)
  b = deepcopy(a)
  s = Int(max(a)[1])
  for i in 1:length(b.vec)
    if b.vec[i].weight[s] >= 0
      if b.vec[i].weight == b.vec[i - 1].weight
        b.vec[i - 1].a += b.vec[i].a
        deleteat!(b.vec, i)
      end
      break
    end

    reflect!(b.vec[i].weight, s)
    lmul!(b.vec[i].w, s)
  end

  return b
end

function sequence(a::LSFanElem, rdec::Vector{<:Integer}=word(max(a)))
  #@req max(a) == WeylGroupElem(weyl_group(root_system(parent(a))) , rdec) "not a reduced decomposition for the maximal element in supp(a)"
  b = deepcopy(a)

  seq = Tuple{Int,Int}[]
  for s in rdec
    w = weight(b)
    for i in 1:length(b.vec)
      if b.vec[i].weight[Int(s)] == 0
        continue
      elseif b.vec[i].weight[Int(s)] > 0
        if i > 1 && b.vec[i].weight == b.vec[i - 1].weight
          b.vec[i - 1].a += b.vec[i].a
          deleteat!(b.vec, i)
        end
        break
      end

      reflect!(b.vec[i].weight, Int(s))
      lmul!(b.vec[i].w, Int(s))
    end

    push!(seq, (Int(s), (weight(b)[Int(s)] - w[Int(s)]) / 2))
  end

  return seq
end

# w highest weight
function sequence(w::WeightLatticeElem, word::Vector{UInt8})
  v = zeros(Int, length(w.vec))
  w2 = deepcopy(w)
  for s in Iterators.reverse(word)
    v[Int(s)] += w2[Int(s)]
    reflect!(w2, Int(s))
  end
  return v
end

###############################################################################
# SequenceIterator

struct SequenceIterator
  LS::LSFan
  seq::Vector{Tuple{Int,Int}}
end

function SequenceIterator(a::LSFanElem)
  LS = parent(a)
  seq = sequence(a)

  exp = zeros(Int, rank(root_system(LS)))
  seen = zeros(Int, rank(root_system(LS)))
  for i in length(seq):-1:1
    exp[seq[i][1]] += seq[i][2]
    seen[seq[i][1]] = i
  end

  return SequenceIterator(
    LS, [(seq[i][1], seen[seq[i][1]] == i ? exp[seq[i][1]] : 0) for i in 1:length(seq)]
  )
end

function Base.IteratorSize(::Type{SequenceIterator})
  return Base.SizeUnknown()
end

function Base.iterate(iter::SequenceIterator)
  return iter.seq, iter.seq
end

function Base.iterate(iter::SequenceIterator, seq::Vector{Tuple{Int,Int}})
  new_seq = copy(seq)
  exp = zeros(Int, rank(root_system(iter.LS)))

  for i in length(seq):-1:1
    s, n = seq[i]
    if n == 0
      continue
    end

    j = findnext(p -> p[1] == s, seq, i + 1)
    if isnothing(j)
      new_seq[i] = (seq[i][1], 0)
      exp[seq[i][1]] += seq[i][2]
      continue
    end

    new_seq[i] = (s, n - 1)
    new_seq[j] = (s, 1)
    while any(exp .> 0)
      i += 1
      new_seq[i] = (seq[i][1], exp[seq[i][1]] + new_seq[i][2])
      exp[seq[i][1]] = 0
    end
    return new_seq, new_seq
  end
end

###############################################################################
# TensorIterator

struct TensorIterator
  LS::LSFan
  a::LSFanElem
  b::LSFanElem

  bound_cache::Dict{Tuple{Int,Vector{Int}},Int}
  cc_cache::Dict{Vector{Int},nf_elem}

  rdec::Vector{Int} #
  exp_a::Vector{Int} # row data
  exp_b::Vector{Vector{Int}} # column data
end

function TensorIterator(a::LSFanElem, b::LSFanElem; rdec::Vector{<:Integer}=word(max(a)))
  @req parent(a) === parent(b) "$a and $b must belong to the same LS fan of monoids"
  @req weight(a) == weight(b) "$a and $b must have the same weight"

  LS = parent(a)
  #l = lcm(
  #  2, map(p -> Int(denominator(p.a)), a.vec)..., map(p -> Int(denominator(p.a)), b.vec)... # TODO: l is not chosen minimal
  #) # == lbar in Littelmann paper
  l = lcm(2, map(p -> Int(denominator(p.a)), b.vec)...)

  seq = sequence(a, rdec)
  rdec = [i for (i, _) in Iterators.reverse(seq)]
  exp_a = [n * l * LS.d[i] for (i, n) in Iterators.reverse(seq)]
  exp_b = Vector{Int}[]
  sizehint!(exp_b, l)

  for p in vec(b)
    n = Int(l * p.a)
    s = sequence(LS.dual_heighst_weight, word(p.w))
    while n > 0
      push!(exp_b, s)
      n -= 1
    end
  end

  return TensorIterator(
    LS,
    a,
    b,
    Dict{Tuple{Int,Vector{Int}},Int}(),
    Dict{Vector{Int},nf_elem}(),
    rdec,
    exp_a,
    exp_b,
  )
end

function TensorIterator(seq::Vector{Tuple{Int, Int}}, b::LSFanElem)
  LS = parent(b)
  l = lcm(2, map(p -> Int(denominator(p.a)), b.vec)...) # TODO: l is not chosen minimal # == lbar in Littelmann paper

  rdec = [i for (i, _) in Iterators.reverse(seq)]
  exp_a = [n * l * LS.d[i] for (i, n) in Iterators.reverse(seq)]
  exp_b = Vector{Int}[]
  sizehint!(exp_b, l)

  for p in vec(b)
    n = Int(l * p.a)
    s = sequence(LS.dual_heighst_weight, word(p.w))
    while n > 0
      push!(exp_b, s)
      n -= 1
    end
  end

  return TensorIterator(
  LS,
  b,
  b,
  Dict{Tuple{Int,Vector{Int}},Int}(),
  Dict{Vector{Int},nf_elem}(),
  rdec,
  exp_a,
  exp_b,
  )
end

# deduplicate with iterate(iter, mat) ?
function Base.iterate(iter::TensorIterator)
  l = length(iter.exp_b)
  mat = zeros(Int, length(iter.rdec), l)
  
  i = 1
  while i <= length(iter.rdec)
    mat[i, :] = compute_bound(iter, mat, i)
    
    s = 0
    for j in 1:l
      s += mat[i, j]
      if s >= iter.exp_a[i]
        mat[i, j] -= s - iter.exp_a[i]
        for k in j+1:l
          mat[i, k] = 0
        end
        break
      end
    end
    
    if s < iter.exp_a[i]
      i -= 1
      np = Vector{Int}()
      while i > 0
        np = next_partition(mat[i, :], compute_bound(iter, mat, i))
        if !isnothing(np)
          break
        end
        i -= 1
      end
      if i == 0
        return nothing
      end
      
      mat[i, :] = np
    end
    
    i += 1
  end

  return mat, mat
end

function Base.iterate(iter::TensorIterator, t::Matrix{Int})
  nt = copy(t)
  for i in (nrows(t) - 2):-1:1
    while true
      nr = next_partition(nt[i, :], compute_bound(iter, nt, i))
      if isnothing(nr)
        break
      end

      nt[i, :] = nr
      r = 0
      for j in (i + 1):nrows(t)
        b = compute_bound(iter, nt, j)
        r = iter.exp_a[j]
        for k in 1:ncols(t)
          m = min(r, b[k])
          nt[j, k] = m
          r -= m
        end
        if r > 0
          break
        end
      end
      if r == 0
        return nt, nt
      end
    end
  end

  return nothing
end


function tensor_coefficient(
  a::LSFanElem, b::LSFanElem; rdec::Vector{<:Integer}=word(max(a))
)
  iter = TensorIterator(a, b; rdec=rdec)
  return _tensor_coefficient(iter)
end

function tensor_coefficient(seq::Vector{Tuple{Int, Int}}, b::LSFanElem)
  iter = TensorIterator(seq, b)
  return _tensor_coefficient(iter)
end

function tensor_coefficient(w::Vector{Int}, s::Vector{Int}, b::LSFanElem)
  iter = TensorIterator([(w[i], s[i]) for i in 1:length(s)], b)
  return _tensor_coefficient(iter)
end

function _tensor_coefficient(iter::TensorIterator)
  coeff = 0
  num = 0

  rk = rank(root_system(iter.LS))
  pbw = [
  GAP.Globals.MonomialElements(
  GAP.Globals.CanonicalBasis(iter.LS.gap_U), GAP.Obj(Int.(i .== 1:rk))
  )[1] for i in 1:rk
  ]
  l = length(iter.exp_b)
  A, q = cyclotomic_field(2*l)

  d = 1
  for p in vec(iter.b)
    n = Int(l * p.a)
    v = GAP.Globals.Basis(iter.LS.gap_V)[1]
    for (i, n) in Iterators.reverse(Iterators.zip(word(p.w), bonds(iter.LS.dual_heighst_weight, p.w)))
      v = pow(iter.LS.gap_U, pbw[Int(i)], n)^v
    end

    rep = GAP.Globals.ExtRepOfObj(GAP.Globals.ExtRepOfObj(v))
    coeffs = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[2])
    d *= sum(k -> coeffs[1][k]*q^(k-1+coeffs[2]), 1:length(coeffs[1]))^n
  end

  for mat in iter
    c = compute_coefficient(iter, mat)
    num += 1
    coeff += c
  end

  return num, A(coeff) / d
end

function pow(U::GAP.GapObj, f::GAP.GapObj, n::Int)
  if n == 0
    return GAP.Globals.One(U)
  end

  fam = GAP.Globals.ElementsFamily(GAP.Globals.FamilyObj(U))
  rep = GAP.Globals.List(GAP.Globals.ExtRepOfObj(f), GAP.Globals.ShallowCopy)
  rep[1][2] = n
  return GAP.Globals.ObjByExtRep(fam, rep)
end

function compute_coefficient(iter::TensorIterator, mat::Matrix{Int})
  l = length(iter.exp_b)
  cm = cartan_matrix(iter.LS.dual_root_system)
  d = cartan_symmetrizer(cm)

  rk = rank(root_system(iter.LS))
  pbw = [
    GAP.Globals.MonomialElements(
      GAP.Globals.CanonicalBasis(iter.LS.gap_U), GAP.Obj(Int.(i .== 1:rk))
    )[1] for i in 1:rk
  ]
  A, v = cyclotomic_field(2*l)

  q = 1
  w = [deepcopy(iter.LS.dual_heighst_weight) for _ in 1:ncols(mat)]
  for i in 1:nrows(mat)
    s = iter.exp_a[i]
    for j in ncols(mat):-1:2
      s -= mat[i, j]
      addmul!(w[j].vec, cm[:, iter.rdec[i]:iter.rdec[i]], -mat[i, j])
      q *= (v^d[iter.rdec[i]])^-(s * (mat[i, j] + w[j].vec[iter.rdec[i]]))
    end
    addmul!(w[1].vec, cm[:, iter.rdec[i]:iter.rdec[i]], -mat[i, 1])
  end

  for i in 1:ncols(mat)
    if haskey(iter.cc_cache, mat[:, i])
      q *= iter.cc_cache[mat[:, i]]
      if q == 0
        break
      end

      continue
    end

    wv = GAP.Globals.Basis(iter.LS.gap_V)[1]
    for l in 1:nrows(mat)
      wv = pow(iter.LS.gap_U, pbw[iter.rdec[l]], mat[l, i])^wv
    end
    rep = GAP.Globals.ExtRepOfObj(GAP.Globals.ExtRepOfObj(wv))
    if iszero(rep)
      q = 0
      iter.cc_cache[mat[:, i]] = zero(A)
      break
    end

    coeffs = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[2])
    coeff = sum(l -> coeffs[1][l] * v^(l - 1 + coeffs[2]), 1:length(coeffs[1]))
    iter.cc_cache[mat[:, i]] = coeff
    q *= coeff
  end

  return q
end

function compute_bound(iter::TensorIterator, mat::Matrix{Int}, row::Int)
  b = zeros(Int, ncols(mat))
  for j in 1:ncols(mat)
    if haskey(iter.bound_cache, (j, mat[1:(row - 1), j]))
      b[j] = iter.bound_cache[(j, mat[1:(row - 1), j])]
      continue
    end

    s = copy(iter.exp_b[j])
    v = GAP.Globals.Basis(iter.LS.gap_V)[1] # highest weight vector
    for i in 1:(row - 1)
      s[iter.rdec[i]] -= mat[i, j]
      for l in 1:mat[i, j]
        v = GAP.Globals.Falpha(iter.LS.gap_V, v, iter.rdec[i])
      end
    end

    k = 0
    while k < s[iter.rdec[row]]
      k += 1
      temp = GAP.Globals.Falpha(iter.LS.gap_V, v, iter.rdec[row])
      if iszero(temp)
        k -= 1
        break
      end
      v = temp
    end

    b[j] = k
    iter.bound_cache[(j, mat[1:(row - 1), j])] = b[j]
  end

  return b
end

# p partition, b bound
function next_partition(p::Vector{Int}, b::Vector{Int})
  np = copy(p)
  len = length(np)

  r = np[len] + 1
  np[len] = 0
  for i in (len - 1):-1:1
    if np[i] != 0
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
  end

  return nothing
end

###############################################################################
# ???

function leq_points_of_weight(a::LSFanElem)
  points = LSFanElem[]
  for seq in SequenceIterator(a)
    b = LSFanElem(parent(a), (1, one(parent(max(a)))); check=false)
    for (i, n) in Iterators.reverse(seq)
      for _ in 1:n
        fi!(b, i)
        if iszero(b)
          @goto next
        end
      end
    end
    if b <= a && !(b in points)
      push!(points, b)
    end

    @label next
  end

  return points
end

function sequence(a::LSFanElem, b::LSFanElem; rdec=word(max(a)))
  b2 = deepcopy(b)
  for (i, n) in sequence(a, rdec)
    if isnothing(ei!(b2, i, n))
      return false
    end
  end

  return true
end

function adapted_string(a::LSFanElem, rdec::Vector{Int})
  s = zero(rdec)

  b = deepcopy(a)
  for i in 1:length(rdec)
    s[i] = eps(b, rdec[i])
    ei!(b, rdec[i], s[i])
  end
  return s
end

function global_eps(a::LSFanElem, i::Integer)
  n = QQFieldElem(0)
  for j in 1:length(a.vec)+1
    if j <= length(a.vec) && a.vec[j].weight[i] < 0
      n -= a.vec[j].a*a.vec[j].weight[i]
    end
  end
  return floor(Int, n)
end

function global_string(a::LSFanElem, rdec::Vector{Int})
  s = zero(rdec)
  b = deepcopy(a)

  n = QQFieldElem(0)
  for i in 1:length(rdec)
    for j in 1:length(b.vec)+1
      if j <= length(b.vec) && b.vec[j].weight[rdec[i]] < 0
        n -= b.vec[j].a*b.vec[j].weight[rdec[i]]
      else
        nf = floor(n)
        s[i] += Int(nf)
        zero!(n)

        k = j - 1
        while k >= 1 && b.vec[k].weight[rdec[i]] < 0
          a = b.vec[k].a*b.vec[k].weight[rdec[i]] # a is negative
          if nf + a > 0
            nf += a
            lmul!(b.vec[k].w, rdec[i])
            reflect!(b.vec[k].weight, rdec[i])
          elseif iszero(nf + a)
            lmul!(b.vec[k].w, rdec[i])
            reflect!(b.vec[k].weight, rdec[i])
            break
          else
            insert!(b.vec, k+1, lsPair(0, lmul(b.vec[k].w, rdec[i]), reflect(b.vec[k].weight, rdec[i])))
            b.vec[k+1].a = 1//b.vec[k+1].weight[rdec[i]]*nf
            b.vec[k].a -= b.vec[k+1].a
            break
          end
          k -= 1
        end
      end
    end
  end

  return s
end

# ============================================================

# number of partitions of n into m nonnegative numbers leq to bound
function partition(n::Int, m::Int, bound::Int = n)
  d = min(n, bound)

  if (m-1)*bound > n-d
    return sum(partition(n-i, m-1, min(n-i, bound)) for i in 0:d)
  elseif (m-1)*bound == n-d
    return 1
  else
    return 0
  end
end

###############################################################################
# Foldings

function down_foldings(a::LSFanElem, dir::Int)
  fds = []

  if length(a.vec) < 2
    return fds
  end

  for i in 1:length(a.vec)
    # find local minimum
    if (i > 1 && a.vec[i-1].w[dir] <= a.vec[i].w[dir]) && (i < length(a.vec) && a.vec[i].w[dir] >= a.vec[i+1].w[dir])
      continue
    end

    for j in i+1:length(a.vec)
      # find local minimum
      if (j > 1 && a.vec[j-1].w[dir] <= a.vec[j].w[dir]) && (j < length(a.vec) && a.vec[j].w[dir] >= a.vec[j+1].w[dir])
        continue
      end
    end
  end
end

function is_balanced(p::LSFanElem)
  for w in reduced_expressions(longest_element(weyl_group(root_system(parent(p)))))
    t = deepcopy(p)
    for i in w
      ii = Int(i)
      n = eps(t, ii)
      if global_eps(t, ii) != n
        return false
      end
      ei!(t, ii, n)
    end
  end
  return true
end

function test_idea(pts::Vector{LSFanElem})
  for p in pts
    for w in reduced_expressions(longest_element(weyl_group(parent(p))))
      iw = Int.(w)
      as = adapted_string(p, iw)
      if global_string(p, iw) == as
        _, c = tensor_coefficient(iw, as, p)
        if c != 1
          return p
        end
        break
      end
    end
  end
end
