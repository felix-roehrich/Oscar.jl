function Nemo.exponent_vector!(
  z::Vector{Int}, a::AbstractAlgebra.Generic.MPoly{T}, i::Int
) where {T<:RingElement}
  copy!(z, a.exps[:, i])
  reverse!(z)
  return z
end

function bar!(z::LaurentPolyWrap, x::LaurentPolyWrap)
  reverse!(z.poly, x.poly, length(x.poly))
  z.mindeg = -degree(x.poly) - x.mindeg
  return z
end

function bar!(z::FracFieldElem{T}, x::FracFieldElem{T}) where {T<:LaurentPolyWrap}
  bar!(z.num, x.num)
  bar!(z.den, x.den)
  return z
end

function Base.show(io::IO, U::QuantumGroup)
  #@show_name(io, W)
  #@show_special(io, W)
  io = pretty(io)
  if is_terse(io)
    print(io, LowercaseOff(), "Quantum group")
  else
    print(io, LowercaseOff(), "Quantum group over ", Lowercase(), root_system(U))
  end
end

function root_system(U::QuantumGroup)
  return U.root_system
end

function weyl_group(U::QuantumGroup)
  return weyl_group(root_system(U))
end

function (U::QuantumGroup)(x::PBWAlgElem)
  @req parent(x) == U.alg "parent mismatch"
  return QuantumGroupElem(U, x)
end

function one(U::QuantumGroup)
  return QuantumGroupElem(U, one(U.alg))
end

function zero(U::QuantumGroup)
  return QuantumGroupElem(U, zero(U.alg))
end

function coefficient_ring(U::QuantumGroup)
  return coefficient_ring(U.alg)
end

###############################################################################
# Generators

function gen(U::QuantumGroup{T,S}, i::Int) where {T,S}
  return QuantumGroupElem{T,S}(U, gen(U.alg, i))
end

function gens(U::QuantumGroup)
  return [gen(U, i) for i in 1:ngens(U)]
end

function chevalley_gens(U::QuantumGroup)
  return [gen(U, U.cvx[i]) for i in 1:rank(root_system(U))]
end

function ngens(U::QuantumGroup)
  return ngens(U.alg)
end

###############################################################################

function quantum_parameter(U::QuantumGroup)
  return U.q
end

@doc raw"""
    quantum_parameter(U::QuantumGroup, i::Int) -> Int
    
Return the quantum parameter for the `i`-th positive root (ordered in convex order).
"""
function quantum_parameter(U::QuantumGroup, i::Int)
  return U.qi[i]
end

function expressify(x::QuantumGroupElem; context=nothing)
  expr = Expr(:call, :+)
  for t in terms(x.elem)
    exp = leading_exponent_vector(t)
    coeff = leading_coefficient(t)
    for i in 1:length(exp)
      coeff = mul!(coeff, quantum_factorial(exp[i], parent(x).qi[i]))
    end
    push!(
      expr.args,
      Expr(:call, :*, expressify(coeff; context),
        (
          string("F[$i]", exp[i] > 1 ? "^($(exp[i]))" : "") for
          i in 1:length(exp) if exp[i] > 0
        )...,
      ),
    )
  end
  return expr
end

@enable_all_show_via_expressify QuantumGroupElem

function Base.deepcopy_internal(x::QuantumGroupElem, dict::IdDict)
  return get!(dict, x, QuantumGroupElem(parent(x), deepcopy_internal(x.elem, dict)))
end

function Base.hash(x::QuantumGroupElem, h::UInt)
  b = 0xe2cb30215ca391a1 % UInt
  h = hash(parent(x), h)
  h = hash(x.elem, h)

  return xor(h, b)
end

function parent(x::QuantumGroupElem)
  return x.U
end

###############################################################################

function iszero(x::QuantumGroupElem)
  return iszero(x.elem)
end

function zero(x::QuantumGroupElem)
  return zero(parent(x))
end

function zero!(x::QuantumGroupElem)
  x.elem = zero!(x.elem)
  return x
end

function one(x::QuantumGroupElem)
  return one(parent(x))
end

function length(x::QuantumGroupElem)
  return length(x.elem)
end

###############################################################################

function add!(z::QuantumGroupElem, x::QuantumGroupElem, y::QuantumGroupElem)
  z.elem = add!(z.elem, x.elem, y.elem)
  return z
end

function add!(x::QuantumGroupElem, y::QuantumGroupElem)
  return add!(x, x, y)
end

function addmul!(z::QuantumGroupElem{T}, x::QuantumGroupElem{T}, a::T) where {T}
  z.elem = addmul!(z.elem, x.elem, a)
  return z
end

function div!(z::QuantumGroupElem{T}, x::QuantumGroupElem{T}, a::T) where {T}
  z.elem = mul!(z.elem, x.elem, inv(a)) # TODO: use div!, when possible
  return z
end

function div!(x::QuantumGroupElem{T}, a::T) where {T}
  return div!(x, x, a)
end

function mul!(z::QuantumGroupElem, x::QuantumGroupElem, y::QuantumGroupElem)
  z.elem = mul!(z.elem, x.elem, y.elem)
  return z
end

function mul!(
  z::QuantumGroupElem{T}, x::QuantumGroupElem{T}, a::Union{Integer,Rational,AbstractFloat,T}
) where {T}
  z.elem = mul!(z.elem, x.elem, a)
  return z
end

function neg!(z::QuantumGroupElem, x::QuantumGroupElem)
  z.elem = neg!(z.elem, x.elem)
  return z
end

function neg!(x::QuantumGroupElem)
  return neg!(x, x)
end

function sub!(z::QuantumGroupElem, x::QuantumGroupElem, y::QuantumGroupElem)
  z.elem = sub!(z.elem, x.elem, y.elem)
  return z
end

function sub!(x::QuantumGroupElem, y::QuantumGroupElem)
  return sub!(x, x, y)
end

function submul!(z::QuantumGroupElem{T}, x::QuantumGroupElem{T}, a::T) where {T}
  z.elem = submul!(z.elem, x.elem, a)
  return z
end

###############################################################################

function Base.:*(x::QuantumGroupElem{T}, y::QuantumGroupElem{T}) where {T<:FieldElem}
  @req parent(x) == parent(y) "parent mismatch"
  return mul!(deepcopy(x), y)
end

function Base.:*(
  x::QuantumGroupElem{T}, a::Union{Integer,Rational,AbstractFloat,T}
) where {T<:FieldElem}
  return mul!(deepcopy(x), a)
end

function Base.:/(x::QuantumGroupElem{T}, a::T) where {T<:FieldElem}
  return div!(deepcopy(x), a)
end

function Base.://(x::QuantumGroupElem{T}, a::T) where {T<:FieldElem}
  return div!(deepcopy(x), a)
end

function Base.:*(
  a::Union{Integer,Rational,AbstractFloat,T}, x::QuantumGroupElem{T}
) where {T<:FieldElem}
  return mul!(deepcopy(x), a)
end

function Base.:+(x::QuantumGroupElem, y::QuantumGroupElem)
  return add!(deepcopy(x), y)
end

function Base.:-(x::QuantumGroupElem, y::QuantumGroupElem)
  return sub!(deepcopy(x), y)
end

function Base.:-(x::QuantumGroupElem)
  return neg!(deepcopy(x))
end

function Base.:^(x::QuantumGroupElem, n::Int)
  return QuantumGroupElem(parent(x), x.elem^n)
end

function Base.:(==)(x::QuantumGroupElem, y::QuantumGroupElem)
  return parent(x) == parent(y) && x.elem == y.elem
end

###############################################################################

function coeff(x::QuantumGroupElem, i::Int)
  return coeff(x.elem, i)
end

function setcoeff!(x::QuantumGroupElem{T}, i::Int, a::T) where {T}
  setcoeff!(x.elem, i, a)
  return x
end

function leading_coefficient(x::QuantumGroupElem)
  return leading_coefficient(x.elem)
end

function leading_exponent_vector(x::QuantumGroupElem)
  return leading_exponent_vector(x.elem)
end

function leading_monomial(x::QuantumGroupElem)
  return QuantumGroupElem(parent(x), leading_monomial(x.elem))
end

function trailing_coefficient(x::QuantumGroupElem)
  return trailing_coefficient(x.elem)
end

###############################################################################
#
#   Quantum Group constructor
#
###############################################################################

function quantum_group(R::RootSystem, w0=word(longest_element(weyl_group(R))))
  A, q = laurent_polynomial_ring(ZZ, "q")
  QA = fraction_field(A)
  P, theta = polynomial_ring(QA, :F => 1:length(w0))

  # for now we rely on the QuaGroup package to compute the PBW relations
  GAP.Packages.load("QuaGroup")
  gapR = GAP.Globals.Objectify(
    GAP.Globals.NewType(
      GAP.Globals.NewFamily(GAP.GapObj("RootSystemFam"), GAP.Globals.IsObject),
      GAP.evalstr("IsAttributeStoringRep and IsRootSystem"),
    ),
    GAP.evalstr("rec()"),
  )

  gapSim = GAP.GapObj(map(r -> GAP.GapObj(coefficients(r))[1], (simple_roots(R))))
  gapPos = GAP.GapObj(map(r -> GAP.GapObj(coefficients(r))[1], (positive_roots(R))))
  gapBil = GAP.GapObj(Oscar.bilinear_form(R))

  GAP.Globals.SetPositiveRoots(gapR, gapPos)
  GAP.Globals.SetNegativeRoots(gapR, -gapPos)
  GAP.Globals.SetSimpleSystem(gapR, gapSim)
  GAP.Globals.SetCartanMatrix(gapR, GAP.Obj(transpose(cartan_matrix(R))))
  GAP.Globals.SetBilinearFormMat(gapR, gapBil)
  GAP.Globals.SetPositiveRootsNF(gapR, gapPos)
  GAP.Globals.SetSimpleSystemNF(gapR, gapSim)
  GAP.Globals.SetBilinearFormMatNF(gapR, gapBil)
  GAP.Globals.SetTypeOfRootSystem(
    gapR, collect(Iterators.flatmap(t -> GAP.GapObj.(t), root_system_type(R)))
  )
  GAP.Globals.SetLongestWeylWord(gapR, GAP.GapObj(GAP.GapObj.(w0)))

  gapU = GAP.Globals.QuantizedUEA(gapR)
  gapF = GAP.Globals.GeneratorsOfAlgebra(gapU)

  # set rels
  npos = number_of_positive_roots(R)
  rels = zero_matrix(P, npos, npos)

  term = one(P)
  for i in 1:npos, j in (i + 1):npos
    rep = Oscar.GAPWrap.ExtRepOfObj(gapF[j] * gapF[i])
    for n in 1:2:length(rep)
      for m in 1:2:length(rep[n])
        for _ in 1:rep[n][m + 1]
          term = mul!(term, theta[rep[n][m]])
        end
        term = mul!(term, inv(quantum_factorial(rep[n][m + 1], QA(q))))
      end
      coeffRep = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[n + 1])
      coeff = zero(QA)
      for n in 1:length(coeffRep[1])
        coeff = addmul!(coeff, QA(q)^(n + coeffRep[2] - 1), coeffRep[1][n])
      end
      rels[i, j] = addmul!(rels[i, j], term, coeff)
      term = one!(term)
    end
  end

  cvx = zeros(Int, npos)
  for i in 1:npos
    beta = w0[i]
    for j in (i - 1):-1:1
      beta = weyl_group(R).refl[w0[j], beta]
    end
    cvx[beta] = i
  end

  alg, _ = pbw_algebra(P, rels, lex(theta))
  return QuantumGroup(
    A,
    alg,
    QA(q),
    [QA(q)^div(r * gapBil * r, 2) for r in GAP.Globals.PositiveRootsInConvexOrder(gapR)],
    R,
    w0,
    cvx,
    Dict{Vector{Int},PBWAlgElem}(),
  )
end

function quantum_integer(n::Int, q::RingElem)
  z = zero(q)
  for i in (-n + 1):2:(n - 1)
    z = add!(z, q^i)
  end
  return z
end

function quantum_factorial(n::Int, q::RingElem)
  z = one(q)
  for i in 1:n
    z = mul!(z, quantum_integer(i, q))
  end
  return z
end

#=
function quantum_binomial(U::QuantumGroup, n::Int, k::Int)
  @req k >= 0 "k must be non-negative"

  z = one(U.q)
  for i in 0:(k - 1)
    z = mul!(z, U.q^(n - i) - U.q^(i - n))
  end
  for i in 1:k
    z = div!(z, U.q^i - U.q^(-i))
  end
  return z
end

=#

###############################################################################
#
#   Printing
#
###############################################################################

function string_representation(x::QuantumGroupElem{T}) where {T}
  U = parent(x)
  rep = Tuple{T,Vector{Int}}[]
  y = deepcopy(x.elem.sdata)

  npos = number_of_positive_roots(root_system(U))
  refl = weyl_group(root_system(U)).refl
  cvx = zeros(Int, npos)
  for i in 1:npos
    beta = U.w0[i]
    for j in (i - 1):-1:1
      beta = refl[U.w0[j], beta]
    end
    cvx[beta] = i
  end

  F = gens(parent(y))
  f = zero(parent(y))
  while !iszero(y)
    t = Singular.trailing_term(y)
    exp = leading_exponent_vector(t)
    s = adapted_string(lusztig_datum(root_system(U), exp, U.w0))

    f = one!(f)
    for i in 1:length(U.w0)
      f = mul!(f, F[cvx[U.w0[i]]]^s[i])
      f = mul!(f, inv!(quantum_factorial(s[i], U.qi[cvx[U.w0[i]]])))
    end

    coeff = coefficient_ring(U)(trailing_coefficient(t)//trailing_coefficient(f))
    push!(rep, (coeff, s))
    y = submul!(y, f, coeff)
  end

  return rep
end

function canonical_basis_representation(x::QuantumGroupElem{T}) where {T<:FieldElem}
  U = parent(x)
  rep = Tuple{T,Vector{Int}}[]

  y = deepcopy(x.elem.sdata)
  while !iszero(y)
    t = Singular.trailing_term(y)
    exp = leading_exponent_vector(t)
    can = _canonical_basis_elem(U, exp).sdata

    coeff = coefficient_ring(U)(trailing_coefficient(t)//trailing_coefficient(can))
    push!(rep, (coeff, exp))
    y = submul!(y, can, coeff)
  end

  return rep
end

###############################################################################
#
#   Internal functions
#
###############################################################################
