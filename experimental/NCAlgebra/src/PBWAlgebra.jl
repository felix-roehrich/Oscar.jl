mutable struct PBWAlgElemData
  poly

  exp::Vector{Int} # exponent vector of leading monomial
  nzf::Int # index of first non-zero entry in exp
  nzl::Int # index of last non-zero entry in exp
end

struct PBWAlg
  C
  D
  MT::Vector{Matrix{PBWAlgElemData}}
end

function ngens(A::PBWAlg)
  return 0
end

# i must be less than j
function _offset(A::PBWAlg, i::Int, j::Int)
  return (i-1)*gens(A)+j-1
end

function zero(x::PBWAlgElemData)
  return PBWAlgebraElemData(
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

struct PBWAlgElem
  parent::PBWAlg
  data::PBWAlgElemData
end

function parent(x::PBWAlgElem)
  return x.parent
end

#define UPMATELEM(i,j,nVar) ( (nVar * ((i)-1) - ((i) * ((i)-1))/2 + (j)-1)-(i) )
function _multiplication(rels::Vector, C::ZZMatrix, D::ZZMatrix)
  N = 0 # number of vars
  MTsize = [] # size of the multiplication table
  R::Ring # polynomial ring
  MT = zeros(R, N) # multiplication table
  
  p = one(R)
  for i in 1:N
    for j in i+1:N
      n = div(i*(i-1), 2)+j-1
      # quasi-commuative case
      if is_zero_entry(D, n)
        MTsize[n] = 1
        MT[n] = identity_matrix(ZZ, 2)
      else
        const defaultSize = 7
        MTsize[n] = defaultSize
        MT[n] = zero_matrix(ZZ, defaultSize, defaultSize)
      end
      
      # purely non-commutative case
      if !is_zero_entry(C, n)
        MT[n][1,1] = rels[n]
      end
    end
  end
end

function mul!(z::PBWAlgebraElem, x::PBWAlgebraElem, y::PBWAlgebraElem)
  @req parent(z) == parent(x) == parent(y) "parent mismatch"
  _mul_p_p!(z.data, x.data, y.data)
  return z
end

# multiply polynomials x and y and store the result in z
function _mul_p_p!(z::PBWAlgebraElemData, x::PBWAlgebraElemData, y::PBWAlgebraElemData)
  z = zero!(z)
  r = zero(z)
  t = zero(z)
  for poly in terms(y.poly)
    t.poly = poly
    t.exp = leading_exponent_vector(poly)
    t.nzf = findfirst(!is_zero, t.exp)
    t.nzl = findlast(!is_zero, t.exp)
    z = add!(z, _mul_p_m!(r, x, t))
  end
end

# multiply a polynomial with a monomial
function _mul_p_m!(z::PBWAlgebraElemData, x::PBWAlgebraElemData, y::PBWAlgebraElemData)
  z = zero!(z)
  r = zero(z)
  t = zero(z)
  for poly in terms(x.poly)
    t.poly = poly
    t.exp = leading_exponent_vector(poly)
    t.nzf = findfirst(!is_zero, t.exp)
    t.nzl = findlast(!is_zero, t.exp)
    z = add!(z, _mul_m_m!(r, t, y))
  end
end

function _mul_m_p!(z::PBWAlgebraElemData, x::PBWAlgebraElemData, y::PBWAlgebraElemData)
  z = zero!(z)
  r = zero(z)
  t = zero(z)
  for poly in terms(y.poly)
    t.poly = poly
    t.exp = leading_exponent_vector(poly)
    t.nzf = findfirst(!is_zero, t.exp)
    t.nzl = findlast(!is_zero, t.exp)
    z = add!(z, _mul_m_m!(r, x, t))
  end
end

# multiply monomials x and z and store result in z
function _mul_m_m!(z::PBWAlgebraElemData, x::PBWAlgebraElemData, y::PBWAlgebraElemData)
  # monomials are ordered
  if x.nzl <= y.nzf
    z.poly = mul!(z.poly, x.poly, y.poly)
    z.exp = add!(z.exp, x.exp, y.exp)
    z.nzf = x.nzf
    z.nzl = y.nzl
    return z
  end
  
  # y is a univariate monomial
  if y.nzf == y.nzl
    return _mul_m_u!(z, x, y)
  end
  
  nxt = zeros(Int, length(y_exp))
  next_nz = 0 # number of non-zero terms
  
  prv = zeros(Int, length(y_exp))
  prv_nz = 0 # number of non-zero terms
  
  for i in 1:xi
    prv[i] = y_exp[i]
    if !is_zero(y_exp[i])
      prv_nz += 1
    end
  end
  
  for i in xi+1:length(y_exp)
    nxt[i] = y_exp[i]
    if !is_zero(y_exp[i])
      next_nz += 1
    end
  end
end

# multiply a univaviate monomial with a monomial and store the result in z
function _mul_u_m!(z::PBWAlgebraElemData, x::PBWAlgebraElemData, y::PBWAlgebraElemData)
  if is_zero(y)
    return z
  end
  
  x_exp = leading_exponent(x.poly)
  y_exp = leading_exponent(y.poly)
  
  xi = findfirst(!is_zero, x_exp)
  y_first = findfirst(!is_zero, y_exp)
  y_last = findlast(!is_zero, y_exp)
  
  if xi <= y_first
    add!(z.poly, x.poly * y.poly)
    return z
  end
  
  # y is also univariate
  if y_first == y_last
    add!(z, parent(z).MT[xi, y_first])
    return z
  end
  
  nxt = zeros(Int, length(y_exp))
  next_nz = 0 # number of non-zero terms
  
  prv = zeros(Int, length(y_exp))
  prv_nz = 0 # number of non-zero terms
  
  for i in 1:xi
    prv[i] = y_exp[i]
    if !is_zero(y_exp[i])
      prv_nz += 1
    end
  end
  
  for i in xi+1:length(y_exp)
    nxt[i] = y_exp[i]
    if !is_zero(y_exp[i])
      next_nz += 1
    end
  end
  
  if nxt_nz == 1
  end
end

function _mul_gens(z, i::Int, n::Int, j::Int, m::Int)
  A = parent(z)
  g = gens(A)
  if i <= j
    return mul!(z, g[i], g[j])
  end
  
  # quasi-commutative case
  if iszero(A.D[j, i])
    z = mul!(z, g[j], g[i])
    return mul!(z, A.C[_offset(A, j, i)]^(n*m))
  end
  
  # current and required multiplication table size
  curSize = A.MTSize[_offset(A, j, i)]
  reqSize = max(a, b)
  
  idx = _offset(A, j, i)
  z = zero!(z)
  if(curSize >= reqSize)
    if !iszero(A.MT[idx][m, n])
      return add!(z, A.MT[idx][m, n])
    end
  else
    const defaultGrow = 5
    newSize = reqSize + defaultGrow
  end
  
  for k in 2:n
    if iszero(A.MT[idx][k, 1])
      _mul_p_m!(A, A.MT[idx][k, 1], A.MT[idx][k-1, 1], gen(A.poly, i))
    end
  end
  
  for k in 2:m
    if iszero(A.MT[idx][n, k])
      _mul_m_p!(A, A.MT[idx][n, k], gen(A.poly, j), A.MT[idx][n, k-1])
    end
  end
  
  return add!(z, A.MT[idx][m, n])
end
