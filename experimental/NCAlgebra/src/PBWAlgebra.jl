const pbwAlg_multGrow = 5

@enum NCType begin
  NCTypeLie
  NCTypeSkew
  NCTypeComm
end

struct PBWAlg{T<:FieldElem} <: NCRing
  R::MPolyRing{T}
  rels::Vector{MPolyRingElem{T}}
  mult::Vector{Matrix{MPolyRingElem{T}}}

  function PBWAlg(R::MPolyRing{T}, rels::Vector{MPolyRingElem{T}}) where {T<:FieldElem}
    mult = Vector{Matrix{MPolyRingElem{T}}}(undef, length(rels))
    for i in 1:length(rels)
      if length(rels[i]) == 1 # quasi-commuative case
        mult[i] = Matrix{MPolyRingElem{T}}(undef, 1, 1)
      else
        mult[i] = Matrix{MPolyRingElem{T}}(undef, pbwAlg_multGrow, pbwAlg_multGrow)
      end
      mult[i][1, 1] = rels[i]
    end

    return new{T}(R, rels, mult)
  end
end

function pbw_algebra(R::MPolyRing{T}, rels::Vector{MPolyRingElem{T}}) where {T<:FieldElem}
  return PBWAlg(R, rels)
end

function pbw_algebra(R::MPolyRing{T}, rels::Matrix{MPolyRingElem{T}}) where {T<:FieldElem}
  N = ngens(R)
  rels2 = [rels[i, j] for i in 1:N, j in (i + 1):N]
  return PBWAlg(R, rels2)
end

function _linear_index(A::PBWAlg, i::Int, j::Int)
  return div((i - 1) * (2 * ngens(A) - i - 2), 2) + j - 1
end

mutable struct PBWAlgebraElem{T<:FieldElem}
  parent::PBWAlg{T}
  poly::MPolyRingElem{T}
end

function elem_type(::Type{PBWAlg{T}}) where {T}
  return PBWAlgebraElem{T}
end

function parent_type(::Type{PBWAlgebraElem{T}}) where {T}
  return PBWAlgebraElem{T}
end

function is_admissible_ordering(
  vars::Vector{MPolyRingElem{T}}, rels::Vector{MPolyRingElem{T}}
) where {T}
  N = length(vars)
  m = Vector{Int}(undef, N)
  lead = zeros(Int, N)
  for i in 1:N, j in (i + 1):N
    k = _linear_index(vars, i, j)
    lead[i], lead[j] = 1, 1
    unpack_monomial!(m, rels[k], 1)
    if monomial_cmp(parent(rels[k]), m, lead) != 0
      return false, (i, j)
    end
    lead[i], lead[j] = 0, 0
  end

  return true, (0, 0)
end

###############################################################################
#
#   PBWAlgebra
#
###############################################################################

function Base.show(io::IO, A::PBWAlg)
  io = pretty(io)
  print(io, LowercaseOff(), "PBW algebra over ", Lowercase(), coefficient_ring(A))
end

function coefficient_ring(A::PBWAlg{T}) where {T}
  return coefficient_ring(A.R)
end

function gen(A::PBWAlg, i::Int)
  return PBWAlgebraElem(A, gen(A.R, i))
end

function gens(A::PBWAlg)
  return [gen(A, i) for i in 1:ngens(A)]
end

function ngens(A::PBWAlg)
  return ngens(A.R)
end

function zero(A::PBWAlg)
  return PBWAlgebraElem(A, zero(A.R))
end

###############################################################################
#
#   PBWAlgebraElem
#
###############################################################################

function parent(x::PBWAlgebraElem)
  return x.parent
end

function Base.show(io::IO, x::PBWAlgebraElem)
  #AbstractAlgebra.combine_like_terms!(x.poly)
  show(io, x.poly)
end

function Base.deepcopy_internal(x::PBWAlgebraElem, dict::IdDict)
  return get!(dict, x, PBWAlgebraElem(parent(x), Base.deepcopy_internal(x.poly, dict)))
end

function Base.hash(x::PBWAlgebraElem, h::UInt)
  b = 0xa080ada44dcea378 % UInt
  h = hash(parent(x), h)
  h = hash(x.poly, h)

  return xor(h, b)
end

###############################################################################

function add!(z::PBWAlgebraElem, x::PBWAlgebraElem, y::PBWAlgebraElem)
  z.poly = add!(z.poly, x.poly, y.poly)
  return z
end

function mul!(
  z::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}
) where {T<:FieldElem}
  _mul_p_p!(parent(z), z.poly, x.poly, y.poly)
  return z
end

function mul!(z::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}, y::T) where {T<:FieldElem}
  z.poly = mul!(z.poly, x.poly, y)
  return z
end

function neg!(z::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}) where {T}
  z.poly = neg!(z.poly, x.poly)
end

function pow!(z::PBWAlgebraElem, x::PBWAlgebraElem, n::Int)
  if is_zero(x)
    return zero!(z)
  end

  if length(x) == 1
    exp = AbstractAlgebra.exponent_vector_ref(x.poly, 1)
    c = 0
    for i in 1:length(exp)
      if exp[i] != 0
        c += 1
        if c > 1
          break
        end
      end
    end
    if c <= 1
      z.poly = pow!(z.poly, x.poly, n)
      return z
    end
  end

  z = one!(z)
  for _ in 1:n
    z = z * x
  end
  return z
end

function sub!(z::PBWAlgebraElem, x::PBWAlgebraElem, y::PBWAlgebraElem)
  z.poly = sub!(z.poly, x.poly, y.poly)
  return z
end

function submul!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}, a::T) where {T}
  x.poly = submul!(x.poly, y.poly, a)
  return x
end

function swap!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T<:FieldElem}
  x.poly, y.poly = y.poly, x.poly
  return x
end

###############################################################################

function Base.:(==)(x::PBWAlgebraElem, y::PBWAlgebraElem)
  @req parent(x) == parent(y) "parent mismatch"
  return x.poly == y.poly
end

function Base.:+(x::PBWAlgebraElem, y::PBWAlgebraElem)
  @req parent(x) == parent(y) "parent mismatch"
  return add!(zero(x), x, y)
end

function Base.:*(x::PBWAlgebraElem, y::PBWAlgebraElem)
  @req parent(x) == parent(y) "parent mismatch"
  return mul!(zero(x), x, y)
end

function Base.:*(x::PBWAlgebraElem{T}, y::T) where {T<:FieldElem}
  return mul!(zero(x), x, y)
end

function Base.:*(x::T, y::PBWAlgebraElem{T}) where {T<:FieldElem}
  return mul!(zero(y), y, x)
end

function Base.:^(x::PBWAlgebraElem, n::Int)
  if n < 0
    throw(DomainError(n, "exponent must be >= 0"))
  elseif n == 0
    return one(x)
  end

  return pow!(zero(x), x, n)
end

###############################################################################

function is_one(x::PBWAlgebraElem)
  return is_one(x.poly)
end

function one(x::PBWAlgebraElem)
  return PBWAlgebraElem(parent(x), one(x.poly))
end

function one!(x::PBWAlgebraElem)
  x.poly = one!(x.poly)
  return x
end

function is_zero(x::PBWAlgebraElem)
  return is_zero(x.poly)
end

function zero(x::PBWAlgebraElem)
  return PBWAlgebraElem(parent(x), zero(x.poly))
end

function zero!(x::PBWAlgebraElem)
  x.poly = zero!(x.poly)
  return x
end

function length(x::PBWAlgebraElem)
  return length(x.poly)
end

function coeff(x::PBWAlgebraElem, i::Int)
  return coeff(x.poly, i)
end

function setcoeff!(x::PBWAlgebraElem{T}, i::Int, a::T) where {T}
  return setcoeff!(x.poly, i, a)
end

function leading_exponent_vector(x::PBWAlgebraElem)
  return leading_exponent_vector(x.poly)
end

###############################################################################
#
#   multiplication implementation
#
###############################################################################

# multiply polynomials x and y and store the result in z
function _mul_p_p!(
  A::PBWAlg{T}, z::MPolyRingElem{T}, x::MPolyRingElem{T}, y::MPolyRingElem{T}
) where {T<:FieldElem}
  if is_one(x)
    z = zero!(z)
    add!(z, y)
    return z
  elseif is_one(y)
    z = zero!(z)
    add!(z, x)
    return z
  end

  z = zero!(z)
  res = zero(z)
  for i in 1:length(x)
    for j in 1:length(y)
      _mul_m_m!(
        A,
        res,
        AbstractAlgebra.exponent_vector_ref(x, i),
        AbstractAlgebra.exponent_vector_ref(y, j),
      )
      mul!(res, coeff(x, i))
      mul!(res, coeff(y, j))
      add!(z, res)
    end
  end
end

function _mul_p_m!(
  A::PBWAlg{T}, z::MPoly{T}, x::MPoly{T}, y::AbstractVector{Int}
) where {T<:FieldElem}
  z = zero!(z)
  res = zero(z)
  for i in 1:length(x)
    _mul_m_m!(A, res, AbstractAlgebra.exponent_vector_ref(x, i), y)
    mul!(res, coeff(x, i))
    add!(z, res)
  end
end

function _mul_m_p!(A::PBWAlg, z::MPoly, x::AbstractVector{Int}, y::MPoly)
  z = zero!(z)
  res = zero(z)
  for j in 1:length(y)
    _mul_m_m!(A, res, x, AbstractAlgebra.exponent_vector_ref(y, j))
    mul!(res, coeff(y, j))
    add!(z, res)
  end
end

# multiply monomials x and z and store result in z
function _mul_m_m!(
  A::PBWAlg, z::MPolyRingElem, x::AbstractVector{Int}, y::AbstractVector{Int}
)
  xl = findlast(!iszero, x)
  if isnothing(xl)
    z = one!(z)
    AbstractAlgebra.add_exponent_vector!(z, 1, y)
    return z
  end

  yf = findfirst(!iszero, y)
  if isnothing(yf)
    z = one!(z)
    AbstractAlgebra.add_exponent_vector!(z, 1, x)
    return z
  end

  # monomials are ordered, no need for exchange relations
  if xl <= yf
    one!(z)
    AbstractAlgebra.add_exponent_vector!(z, 1, x)
    AbstractAlgebra.add_exponent_vector!(z, 1, y)
    return z
  end

  tmp = zero(z)
  mon = zeros(Int, ngens(A))

  # apply exchange relation
  _mul_gens(A, z, xl, x[xl], yf, y[yf])
  xl = findprev(!iszero, x, xl - 1)
  yf = findnext(!iszero, y, yf + 1)

  if !isnothing(yf)
    copyto!(mon, yf, y, yf, length(y) - yf + 1)
    _mul_p_m!(A, tmp, z, mon)
  else
    swap!(tmp, z)
  end
  if !isnothing(xl)
    zero!(mon)
    copyto!(mon, 1, x, 1, xl)
    _mul_m_p!(A, z, mon, tmp)
  else
    swap!(z, tmp)
  end
  return z
end

# i > j
function _mul_gens(
  A::PBWAlg{T}, z::MPolyRingElem{T}, i::Int, n::Int, j::Int, m::Int
) where {T<:FieldElem}
  ind = _linear_index(A, j, i)

  # quasi-commutative case
  if length(A.rels[ind]) == 1
    one!(z)
    AbstractAlgebra.add_exponent!(z, 1, i, n)
    AbstractAlgebra.add_exponent!(z, 1, j, m)
    cf = coeff(A.mult[ind][1, 1], 1)
    if !is_one(cf)
      pow!(coeff(z, 1), cf, n * m)
    end
    return z
  end

  # current and required multiplication table size
  curSize = size(A.mult[ind], 1)
  reqSize = max(n, m)

  z = zero!(z)
  if (curSize >= reqSize)
    if isassigned(A.mult[ind], n, m)
      return add!(z, A.mult[ind][n, m])
    end
  else
    newSize = reqSize + pbwAlg_multGrow
    mult = Matrix{MPolyRingElem{T}}(undef, newSize, newSize)
    for k in 1:curSize
      copyto!(mult, (k - 1) * newSize + 1, A.mult[ind], (k - 1) * curSize + 1, curSize)
    end
    A.mult[ind] = mult
  end

  mon = zeros(Int, ngens(A))
  mon[i] = 1
  for k in 2:n
    if !isassigned(A.mult[ind], k, 1)
      A.mult[ind][k, 1] = zero(A.R)
      _mul_m_p!(A, A.mult[ind][k, 1], mon, A.mult[ind][k - 1, 1])
    end
  end

  mon[i], mon[j] = 0, 1
  for k in 2:m
    if !isassigned(A.mult[ind], n, k)
      A.mult[ind][n, k] = zero(A.R)
      _mul_p_m!(A, A.mult[ind][n, k], A.mult[ind][n, k - 1], mon)
    end
  end

  return add!(z, A.mult[ind][n, m])
end

###############################################################################
#
#   division
#
###############################################################################

function monomial_lcm!(z::Vector{Int}, x::AbstractVector{Int}, y::AbstractVector{Int})
  for i in 1:length(z)
    z[i] = max(x[i], y[i])
  end
end

function monomial_divides(
  z::AbstractVector{Int}, x::AbstractVector{Int}, y::AbstractVector{Int}
)
  for i in 1:length(x)
    z[i] = x[i] - y[i]
    if z[i] < 0
      return false
    end
  end
  return true
end

function divexact!(x::PBWAlgebraElem{T}, a::T) where {T}
  for i in 1:length(x)
    AbstractAlgebra.divexact!(coeff(x, i), a)
  end
  return x
end

function rem(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  return rem!(zero(x), x, [y])
end

function rem!(
  r::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}, yi::Vector{PBWAlgebraElem{T}}
) where {T}
  A = parent(x)
  m = one(x)
  q = zero(x)

  r = set!(r, x)
  em = AbstractAlgebra.exponent_vector_ref(m.poly, 1)

  i = 1
  while i <= length(r) && !iszero(r)
    j = 1
    while j <= length(yi) &&
      !monomial_divides(
        em,
        AbstractAlgebra.exponent_vector_ref(r.poly, i),
        AbstractAlgebra.exponent_vector_ref(yi[j].poly, 1),
      )
      j += 1
    end
    if j > length(yi)
      i += 1
      continue
    end

    AbstractAlgebra.setcoeff!(m.poly, 1, coeff(r, i))
    _mul_m_p!(A, q.poly, em, yi[j].poly)
    AbstractAlgebra.divexact!(coeff(m, 1), coeff(q, 1))
    submul!(r, q, coeff(m, 1))
  end

  return r
end

###############################################################################
#
#   Groebner basis
#
###############################################################################

struct PBWAlgebraIdeal{T} <: AbstractAlgebra.Ideal{T}
  gens::Vector{PBWAlgebraElem{T}}
  gb::Vector{PBWAlgebraElem{T}}
end

function Base.show(io::IO, I::PBWAlgebraIdeal)
  io = pretty(io)
  print(io, LowercaseOff(), "PBW algebra ideal")
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

function short_s_polynomial(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  z = zeros(Int, ngens(parent(x)))
  monomial_lcm!(
    z, AbstractAlgebra.exponent_vector_ref(x.poly, 1),
    AbstractAlgebra.exponent_vector_ref(y.poly, 1),
  )
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

#=
function groebner_basis(I::PBWAlgebraIdeal)
  struct _lobject
    x::PBWAlgebraElem
    y::PBWAlgebraElem
  end

  struct _strat
    P::Any
    enterS::Any
    Ll::Int
    L::Vector{_lobject}
    Interpt::Bool
    red
    tail
  end

  strat = _strat()
  # init
  strat.red = function(l::_lobject, strat::_strat)
    j = 1
    if j > strat.ls
      return 0
    end
    if is_divisible(strat.S[j], l.p)
      l.p = reduce_s_polynomial(strat.S[j], l.p)
    end
    if l.p
      l.lcm = zero!(l.lcm)
      return 0
    end
  end

  while !isempty(strat.L)
    strat.P = pop!(strat.L) # strat.L[strat.Ll]
    strat.red(strat.P, strat)

    if !iszero(strat.P)
      if next(strat.P.p) == strat.tail
      end

      if strat.P.p == strat.tail
        strat.red(strat.P, strat)
      end
    end
  end
end
=#
