const pbwAlg_multGrow = 5

using Base: deepcopy_internal
mutable struct PBWAlgElemData
  poly::MPoly

  exp::Vector{Int} # exponent vector of leading monomial
  nzf::Int # index of first non-zero entry in exp
  nzl::Int # index of last non-zero entry in exp
end

struct _pbw_multiplication_table
  N::Int
  T::Vector{Matrix{PBWAlgElemData}}
end

mutable struct _pbw_monomial
  exp::Vector{Int}
  nzf::Int
  nzl::Int
end

function Base.getindex(mult::_pbw_multiplication_table, i::Int, j::Int)
  return mult.T[(i-1)*mult.N+j-1]
end

struct PBWAlg{T<:FieldElem}
  R::MPolyRing{T}
  rels::Vector{MPolyRingElem{T}}
  mult::_pbw_multiplication_table

  function PBWAlg()
    A, q = laurent_polynomial_ring(ZZ, "q")
    QA = fraction_field(A)
    R, f = polynomial_ring(QA, ("F" => 1:3))

    rels = [q^-1 * f[1] * f[2], q * f[1] * f[3] + f[2], q^-1 * f[2] * f[3]]

    m1 = Array{PBWAlgElemData}(undef, pbwAlg_multGrow, pbwAlg_multGrow)
    m1[1, 1] = PBWAlgElemData(rels[1], [1, 1, 0], 1, 2)
    m2 = Array{PBWAlgElemData}(undef, pbwAlg_multGrow, pbwAlg_multGrow)
    m2[1, 1] = PBWAlgElemData(rels[2], [1, 0, 1], 1, 3)
    m3 = Array{PBWAlgElemData}(undef, pbwAlg_multGrow, pbwAlg_multGrow)
    m3[1, 1] = PBWAlgElemData(rels[3], [0, 1, 1], 2, 3)

    return new{elem_type(QA)}(
      R,
      rels,
      _pbw_multiplication_table(ngens(R), [m1, m2, m3])
    )
  end
end

# i must be less than j
function _offset(A::PBWAlg, i::Int, j::Int)
  return (i - 1) * ngens(A) + j - 1
end

function one(x::PBWAlgElemData)
  return PBWAlgElemData(one(x.poly), zeros(Int, length(x.exp)), 0, 0)
end

function one!(x::PBWAlgElemData)
  x.poly = one!(x.poly)
  x.exp = zero!(x.exp)
  x.nzf = 0
  x.nzl = 0
  return x
end

function add!(x::PBWAlgElemData, y::PBWAlgElemData)
  x.poly = add!(x.poly, y.poly)
  if x.exp < y.exp
    copy!(x.exp, y.exp)
    x.nzf = y.nzf
    x.nzl = y.nzl
  end

  return x
end

function mul!(x::PBWAlgElemData, y::PBWAlgElemData)
  x.poly = mul!(x.poly, y.poly)
  x.exp = add!(x.exp, y.exp)
  x.nzf = min(x.nzf, y.nzf)
  x.nzl = max(x.nzl, y.nzl)

  return x
end

function zero(x::PBWAlgElemData)
  return PBWAlgElemData(
    zero(x.poly),
    zeros(Int, length(x.exp)),
    0,
    0
  )
end

function zero!(x::PBWAlgElemData)
  x.poly = zero!(x.poly)
  x.exp = zero!(x.exp)
  x.nzf = 0
  x.nzl = 0
  return x
end

function leading_term(x::PBWAlgElemData)
  return PBWAlgElemData(
    leading_term(x.poly),
    copy(x.exp),
    x.nzf,
    x.nzl
  )
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

function coefficient_ring(A::PBWAlg)
  return coefficient_ring(A.R)
end

function gen(A::PBWAlg, i::Int)
  exp = zeros(Int, ngens(A))
  exp[i] = 1
  return PBWAlgebraElem(A, PBWAlgElemData(gen(A.R, i), exp, i, i))
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

struct PBWAlgebraElem
  parent::PBWAlg
  data::PBWAlgElemData
end

function parent(x::PBWAlgebraElem)
  return x.parent
end

function Base.show(io::IO, x::PBWAlgebraElem)
  show(io, x.data.poly)
end

function Base.deepcopy_internal(x::PBWAlgebraElem, dict::IdDict)
  return get!(dict, x, PBWAlgebraElem(parent(x), Base.deepcopy_internal(x.data, dict)))
end

function Base.hash(x::PBWAlgebraElem, h::UInt)
  b = 0xa080ada44dcea378 % UInt
  h = hash(parent(x), h)
  h = hash(x.data, h)

  return xor(h, b)
end

function Base.:*(x::PBWAlgebraElem, y::PBWAlgebraElem)
  return mul!(zero(x), x, y)
end

###############################################################################

function add!(z::PBWAlgebraElem, x::PBWAlgebraElem, y::PBWAlgebraElem)
  @req parent(z) == parent(x) == parent(y) "parent mismatch"
  z.data = add!(x.data, y.data)
  return z
end

function mul!(z::PBWAlgebraElem, x::PBWAlgebraElem, y::PBWAlgebraElem)
  @req parent(z) == parent(x) == parent(y) "parent mismatch"
  _mul_p_p!(parent(z), z.data, x.data, y.data)
  return z
end

###############################################################################

function is_one(x::PBWAlgebraElem)
  return is_one(x.data.poly)
end

function one(x::PBWAlgebraElem)
  return PBWAlgebraElem(parent(x), one(x.data))
end

function is_zero(x::PBWAlgebraElem)
  return is_zero(x.data.poly)
end

function zero(x::PBWAlgebraElem)
  return PBWAlgebraElem(parent(x), zero(x.data))
end

function zero!(x::PBWAlgebraElem)
  x.data = zero!(x.data)
  return x
end

function leading_exponent_vector(x::PBWAlgebraElem)
  return x.data.exp
end

function tail(x::PBWAlgElemData)
  return PBWAlgElemData(tail(x.poly), copy(x.exp), x.nzf, x.nzl)
end

function tail(x::PBWAlgebraElem)
  return PBWAlgebraElem(parent(x), tail(x.data))
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
function _mul_p_p!(A::PBWAlg, z::PBWAlgElemData, x::PBWAlgElemData, y::PBWAlgElemData)
  z = zero!(z)

  # temporary storage for monomials
  mx = _pbw_monomial(zeros(Int, length(x.exp)), 0, 0)
  my = _pbw_monomial(zeros(Int, length(y.exp)), 0, 0)
  res = zero(z)
  for py in terms(y.poly)
    for px in terms(x.poly)
      AbstractAlgebra.exponent_vector!(mx.exp, px, 1)
      mx.nzf = findfirst(!is_zero, mx.exp)
      mx.nzl = findlast(!is_zero, mx.exp)

      AbstractAlgebra.exponent_vector!(my.exp, py, 1)
      my.nzf = findfirst(!is_zero, my.exp)
      my.nzl = findlast(!is_zero, my.exp)

      _mul_m_m!(A, res, mx, my)
      mul!(res.poly, coeff(px, 1))
      mul!(res.poly, coeff(py, 1))
      add!(z, res)
    end
  end
end

function _mul_p_m!(A::PBWAlg, z::PBWAlgElemData, x::PBWAlgElemData, y::_pbw_monomial)
  z = zero!(z)
  t = zero(z)
  r = zero(z)
  for poly in terms(x.poly)
    t.poly = poly
    t.exp = leading_exponent_vector(poly)
    t.nzf = findfirst(!is_zero, t.exp)
    t.nzl = findlast(!is_zero, t.exp)
    z = add!(z, _mul_m_m!(A, r, t, y))
  end
end

function _mul_m_p!(A::PBWAlg, z::PBWAlgElemData, x::_pbw_monomial, y::PBWAlgElemData)
  z = zero!(z)
  t = zero(z)
  r = zero(z)
  for poly in terms(y.poly)
    t.poly = poly
    t.exp = leading_exponent_vector(poly)
    t.nzf = findfirst(!is_zero, t.exp)
    t.nzl = findlast(!is_zero, t.exp)
    z = add!(z, _mul_m_m!(A, r, x, t))
  end
end

# multiply monomials x and z and store result in z
function _mul_m_m!(A::PBWAlg, z::PBWAlgElemData, x::_pbw_monomial, y::_pbw_monomial)
  # monomials are ordered
  if x.nzl <= y.nzf
    z.poly = one!(z.poly)
    z.exp = add!(z.exp, x.exp, y.exp)
    AbstractAlgebra.set_exponent_vector!(z.poly, 1, z.exp)
    z.nzf = x.nzf
    z.nzl = y.nzl

    return z
  end

  # temporary storage for monomials
  mon = _pbw_monomial(copy(y.exp), y.nzf, y.nzl)

  # factors where we need to apply exchange relations
  r2 = zero(z)
  r3 = zero(z)
  for i in x.nzl:-1:y.nzf
    if is_zero(x.exp[i])
      continue
    end

    zero!(mon.exp[1:y.nzf-1])
    copy!(mon.exp[y.nzf:end], y.exp[y.nzf:end])
    for j in y.nzf:i-1
      if is_zero(y.exp[j])
        continue
      end
      mon.exp[j] = 0
      mon.nzf = j + 1

      # apply exchange relation
      _mul_gens(A, r2, i, x.exp[i], j, y.exp[j])
      for d in 1:length(r2.poly)
        AbstractAlgebra.add_exponent_vector!(r2.poly, d, mon.exp)
      end
      add!(r3.poly, r2.poly)
      r3.exp, r3.nzf, r3.nzl = r2.exp, r2.nzf, r2.nzl

      l = tail(r2)
      if !is_zero(l.poly)
        _mul_p_m!(A, r2, l, mon)
        add!(r3, r2)
      end
    end

    # first nonzero entry of LM(r1) is i
    copy!(mon.exp[1:i-1], x.exp[1:i-1])
    zero!(mon.exp[i:end])
    for d in 1:length(r3.poly)
      AbstractAlgebra.add_exponent_vector!(r3.poly, d, mon.exp)
    end

    add!(z, r3)
  end
end

# i > j
function _mul_gens(A::PBWAlg, z::PBWAlgElemData, i::Int, n::Int, j::Int, m::Int)
  # quasi-commutative case
  if length(A.rels[_offset(A, j, i)]) == 1
    z = one!(z)
    mul!(z.poly, A.mult[j, i][1, 1].poly)
    copy!(z.exp, A.mult[j, i][1, 1].exp)
    z.nzf, z.nzl = A.mult[j, i][1, 1].nzf, A.mult[j, i][1, 1].nzl
    return
  end

  # current and required multiplication table size
  curSize = size(A.mult[j, i], 1)
  reqSize = max(n, m)

  z = zero!(z)
  if (curSize >= reqSize)
    if !ismissing(A.mult[j, i][m, n])
      return add!(z, A.mult[j, i][m, n])
    end
  else
    newSize = reqSize + pbwAlg_multGrow
    mult = Matrix{MPoly}(undef, newSize, newSize)
    copyto!(mult, A.mult[j, i])
    A.mult = mult
  end

  mon = _pbw_monomial(zeros(Int, length(z.exp)), i, i)
  mon.exp[i] = 1
  for k in 2:n
    if ismissing(A.mult[j, i][k, 1])
      _mul_p_m!(A, A.mult[j, i][k, 1], A.mult[j, i][k-1, 1], mon)
    end
  end

  mon.exp[i], mon.exp[j] = 0, 1
  mon.nzf, mon.nzl = j, j
  for k in 2:m
    if ismissing(A.mult[j, i][n, k])
      _mul_m_p!(A, A.mult[j, i][n, k], mon, A.mult[j, i][n, k-1])
    end
  end

  return add!(z, A.mult[j, i][m, n])
end
