const pbwAlg_multGrow = 5

struct PBWAlg{T<:FieldElem}
  R::MPolyRing{T}
  rels::Vector{MPoly{T}}
  mult::Vector{Matrix{MPoly{T}}}

  function PBWAlg()
    A, q = laurent_polynomial_ring(ZZ, "q")
    QA = fraction_field(A)
    R, f = polynomial_ring(QA, ("F" => 1:3))

    rels = [q^-1 * f[1] * f[2], q * f[1] * f[3] + f[2], q^-1 * f[2] * f[3]]

    m1 = Array{MPoly}(undef, pbwAlg_multGrow, pbwAlg_multGrow)
    m1[1, 1] = rels[1]
    m2 = Array{MPoly}(undef, pbwAlg_multGrow, pbwAlg_multGrow)
    m2[1, 1] = rels[2]
    m3 = Array{MPoly}(undef, pbwAlg_multGrow, pbwAlg_multGrow)
    m3[1, 1] = rels[3]

    return new{elem_type(QA)}(
      R,
      rels,
      [m1, m2, m3]
    )
  end

  function PBWAlg(R::MPolyRing{T}, rels::Vector{MPoly{T}}) where {T<:FieldElem}
    mult = Vector{Matrix{MPoly{T}}}(undef, length(rels))
    for i in 1:length(rels)
      if length(rels[i]) == 1 # quasi-commuative case
        mult[i] = Matrix{MPoly{T}}(undef, 1, 1)
      else
        mult[i] = Matrix{MPoly{T}}(undef, pbwAlg_multGrow, pbwAlg_multGrow)
      end
      mult[i][1, 1] = rels[i]
    end

    return new{T}(R, rels, mult)
  end
end

function pbw_algebra(R::MPolyRing{T}, rels::Vector{MPoly{T}}) where {T<:FieldElem}
  return PBWAlg(R, rels)
end

function pbw_algebra(R::MPolyRing{T}, rels::Matrix{MPoly{T}}) where {T<:FieldElem}
  N = ngens(R)
  rels2 = [rels[i, j] for i in 1:N, j in i+1:N]
  return PBWAlg(R, rels2)
end

function _linear_index(A::PBWAlg, i::Int, j::Int)
  return div((i - 1) * (2 * ngens(A) - i - 2), 2) + j - 1
end

mutable struct PBWAlgebraElem{T<:FieldElem}
  parent::PBWAlg{T}
  poly::MPoly{T}
end

function elem_type(::Type{PBWAlg{T}}) where {T}
  return PBWAlgebraElem{T}
end

function parent_type(::Type{PBWAlgebraElem{T}}) where {T}
  return PBWAlgebraElem{T}
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

function mul!(z::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T<:FieldElem}
  _mul_p_p!(parent(z), z.poly, x.poly, y.poly)
  return z
end

function mul!(z::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}, y::T) where {T<:FieldElem}
  z.poly = mul!(z.poly, x.poly, y)
  return z
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

function swap!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T<:FieldElem}
  x.poly, y.poly = y.poly, x.poly
  return x
end

###############################################################################

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

function leading_exponent_vector(x::PBWAlgebraElem)
  return leading_exponent_vector(x.poly)
end

function tail(x::PBWAlgebraElem)
  return PBWAlgebraElem(parent(x), tail(x.poly))
end

###############################################################################
#
#   multiplication implementation
#
###############################################################################

# multiply polynomials x and y and store the result in z
function _mul_p_p!(A::PBWAlg{T}, z::MPoly{T}, x::MPoly{T}, y::MPoly{T}) where {T<:FieldElem}
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
      _mul_m_m!(A, res, AbstractAlgebra.exponent_vector_ref(x, i), AbstractAlgebra.exponent_vector_ref(y, j))
      if !is_one(coeff(x, i))
        mul!(res, coeff(x, i))
      end
      if !is_one(coeff(y, j))
        mul!(res, coeff(y, j))
      end
      println(res)
      add!(z, res)
    end
  end
end

function _mul_p_m!(A::PBWAlg{T}, z::MPoly{T}, x::MPoly{T}, y::AbstractVector{Int}) where {T<:FieldElem}
  z = zero!(z)
  res = zero(z)
  for i in 1:length(x)
    _mul_m_m!(A, res, AbstractAlgebra.exponent_vector_ref(x, i), y)
    if !is_one(coeff(x, i))
      mul!(res, coeff(x, i))
    end
    add!(z, res)
  end
end

function _mul_m_p!(A::PBWAlg, z::MPoly, x::AbstractVector{Int}, y::MPoly)
  z = zero!(z)
  res = zero(z)
  for i in 1:length(y)
    _mul_m_m!(A, res, x, AbstractAlgebra.exponent_vector_ref(y, i))
    if !is_one(coeff(y, i))
      mul!(res, coeff(y, i))
    end
    add!(z, res)
  end
end

# multiply monomials x and z and store result in z
function _mul_m_m!(A::PBWAlg, z::MPoly, x::AbstractVector{Int}, y::AbstractVector{Int})
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
  mon = zeros(Int, ngens(A) + 1)

  # apply exchange relation
  _mul_gens(A, z, xl, x[xl], yf, y[yf])
  xl = findprev(!iszero, x, xl - 1)
  yf = findnext(!iszero, y, yf + 1)

  if !isnothing(yf)
    copyto!(mon, yf, y, yf, length(y) - yf + 1)
    _mul_p_m!(A, tmp, z, mon)
  else
    AbstractAlgebra.swap!(tmp, z)
  end
  if !isnothing(xl)
    zero!(mon)
    copyto!(mon, 1, x, 1, xl)
    _mul_m_p!(A, z, mon, tmp)
  else
    AbstractAlgebra.swap!(z, tmp)
  end
  return z
end

#=
r1 = zero(z)
r2 = zero(z)
mon = zeros(Int, ngens(A))

# apply exchange relation
i = xl
while i > yf
  if x[i] == 0
    i -= 1
    continue
  end

  j = yf
  copyto!(mon, y)
  while i > j
    if y[j] != 0
      _mul_gens(A, r1, i, x[i], j, y[j])
      r2 = mul!(r2, r1)
      mon[j] = 0
      if length(r1) > 1
        break
      end
    end
    j += 1
  end
  if i > j # we stopped because length(res) > 1
    _mul_p_m!(A, r1, r2, mon)
    AbstractAlgebra.swap!(r2, r1)
  else
    AbstractAlgebra.add_exponent_vector!(r2, 1, mon)
  end

  zero!(mon)
  while x[i] == 0 && i > 0
    i -= 1
  end
  copyto!(mon, 1, x, 1, i)
  if i > yf
    _mul_m_p!(A, r1, mon, r2)
  else
    for k in 1:length(r2)
      add_exponent_vector!(r2, k, mon)
    end
    AbstractAlgebra.swap!(r1, r2)
  end
  z = add!(z, r1)
end
return z
=#

# i > j
function _mul_gens(A::PBWAlg{T}, z::MPoly, i::Int, n::Int, j::Int, m::Int) where {T<:FieldElem}
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
    mult = Matrix{MPoly{T}}(undef, newSize, newSize)
    copyto!(mult, A.mult[ind])
    A.mult[ind] = mult
  end

  mon = zeros(Int, ngens(A) + 1)
  mon[i] = 1
  for k in 2:n
    if !isassigned(A.mult[ind], k, 1)
      A.mult[ind][k, 1] = zero(A.R)
      _mul_m_p!(A, A.mult[ind][k, 1], mon, A.mult[ind][k-1, 1])
    end
  end

  mon[i], mon[j] = 0, 1
  for k in 2:m
    if !isassigned(A.mult[ind], n, k)
      A.mult[ind][n, k] = zero(A.R)
      _mul_p_m!(A, A.mult[ind][n, k], A.mult[ind][n, k-1], mon)
    end
  end

  return add!(z, A.mult[ind][n, m])
end
