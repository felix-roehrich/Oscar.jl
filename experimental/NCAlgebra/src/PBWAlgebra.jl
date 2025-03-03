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

###############################################################################
#
#   PBWAlgebraElem
#
###############################################################################

mutable struct PBWAlgebraElem{T<:FieldElem}
  parent::PBWAlg{T}
  poly::MPoly{T}
end

function parent(x::PBWAlgebraElem)
  return x.parent
end

function Base.show(io::IO, x::PBWAlgebraElem)
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

function mul!(z::PBWAlgebraElem, x::PBWAlgebraElem, y::PBWAlgebraElem)
  _mul_p_p!(parent(z), z.poly, x.poly, y.poly)
  return z
end

function Base.:+(x::PBWAlgebraElem, y::PBWAlgebraElem)
  @req parent(x) == parent(y) "parent mismatch"
  return add!(zero(x), x, y)
end

function Base.:*(x::PBWAlgebraElem, y::PBWAlgebraElem)
  @req parent(x) == parent(y) "parent mismatch"
  return mul!(zero(x), x, y)
end

function Base.:^(x::PBWAlgebraElem, n::Int)
  if n < 0
    throw(DomainError(n, "exponent must be >= 0"))
  elseif n == 0
    return one(x)
  end

  z1 = deepcopy(x)
  z2 = deepcopy(x)
  for i in 1:n
    if isodd(i)
      mul!(z1, z2, x)
    else
      mul!(z2, z1, x)
    end
  end

  return isodd(n) ? z2 : z1
end

###############################################################################

function is_one(x::PBWAlgebraElem)
  return is_one(x.poly)
end

function one(x::PBWAlgebraElem)
  return PBWAlgebraElem(parent(x), one(x.poly))
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

#define UPMATELEM(i,j,nVar) ( (nVar * ((i)-1) - ((i) * ((i)-1))/2 + (j)-1)-(i) )
function _multiplication(rels::Vector, C::ZZMatrix, D::ZZMatrix)
  N = 0 # number of vars
  MTsize = [] # size of the multiplication table
  R::Ring # polynomial ring
  MT = zeros(R, N) # multiplication table

  p = one(R)
  for i in 1:N
    for j in i+1:N
      n = div(i * (i - 1), 2) + j - 1
      # quasi-commuative case
      if is_zero_entry(D, n)
        MTsize[n] = 1
        MT[n] = identity_matrix(ZZ, 2)
      else
        MTsize[n] = pbwAlg_multGrow
        MT[n] = zero_matrix(ZZ, pbwAlg_multGrow, pbwAlg_multGrow)
      end

      # purely non-commutative case
      if !is_zero_entry(C, n)
        MT[n][1, 1] = rels[n]
      end
    end
  end
end

# multiply polynomials x and y and store the result in z
function _mul_p_p!(A::PBWAlg{T}, z::MPoly{T}, x::MPoly{T}, y::MPoly{T}) where {T<:FieldElem}
  z = zero!(z)

  mx = zeros(Int, ngens(A))
  my = zeros(Int, ngens(A))
  res = zero(z)
  for i in 1:length(x)
    AbstractAlgebra.exponent_vector!(mx, x, i)
    for j in 1:length(y)
      AbstractAlgebra.exponent_vector!(my, y, j)

      _mul_m_m!(A, res, mx, my)
      mul!(res, coeff(x, i))
      mul!(res, coeff(y, j))
      add!(z, res)
    end
  end
end

function _mul_p_m!(A::PBWAlg, z::MPoly, x::MPoly, y::Vector{Int})
  z = zero!(z)
  res = zero(z)
  mx = zeros(Int, ngens(A))
  for i in 1:length(x)
    AbstractAlgebra.exponent_vector!(mx, x, i)
    _mul_m_m!(A, res, mx, y)
    addmul!(z, res, coeff(x, i))
  end
end

function _mul_m_p!(A::PBWAlg, z::MPoly, x::Vector{Int}, y::MPoly)
  z = zero!(z)
  res = zero(z)
  my = zeros(Int, ngens(A))
  for i in 1:length(y)
    AbstractAlgebra.exponent_vector!(my, y, i)
    _mul_m_m!(A, res, x, my)
    addmul!(z, res, coeff(y, i))
  end
end

# multiply monomials x and z and store result in z
function _mul_m_m!(A::PBWAlg, z::MPoly, x::Vector{Int}, y::Vector{Int})
  xl = findlast(!iszero, x)
  if isnothing(xl)
    return zero!(z)
  end

  yf = findfirst(!iszero, y)
  if isnothing(yf)
    return zero!(z)
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
    AbstractAlgebra.swap!(tmp, z)
  end
  if !isnothing(xl)
    zero!(mon)
    copyto!(mon, 1, x, 1, xl)
    _mul_m_p!(A, z, mon, tmp)
  else
    AbstractAlgebra.swap!(z, tmp)
  end
end

# i > j
function _mul_gens(A::PBWAlg{T}, z::MPoly, i::Int, n::Int, j::Int, m::Int) where {T<:FieldElem}
  ind = (j - 1) * (ngens(A) - j) + i - 1

  # quasi-commutative case
  if length(A.rels[ind]) == 1
    one!(z)

    AbstractAlgebra.add_exponent!(z, 1, i, n)
    AbstractAlgebra.add_exponent!(z, 1, j, m)
    cf = coeff(A.mult[ind][1, 1], 1)
    if !is_one(cf)
      pow!(coeff(z, 1), cf, n * m)
    end
    return
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

  mon = zeros(Int, ngens(A))
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

  add!(z, A.mult[ind][n, m])
end
