using AbstractAlgebra.Generic: LaurentPolyWrap
import ..Oscar: add!, addmul!, mul!, neg!, sub!, submul!, one
import ..Oscar:
  coefficient_ring, leading_coefficient, leading_exponent_vector, leading_monomial

function Nemo.exponent_vector!(
  z::Vector{Int}, a::AbstractAlgebra.Generic.MPoly{T}, i::Int
) where {T<:RingElement}
  copy!(z, a.exps[:, i])
  reverse!(z)
  return z
end

struct PBWAlgModule
  alg::PBWAlgRing
  sdata::Singular.smodule
end

struct QuantumGroup
  A::Any#::LaurentPolynomialRing

  alg::PBWAlgRing
  gens::Vector{PBWAlgElem}
  q::Any #::RationalFunctionFieldElem{QQFieldElem,QQPolyRingElem}
  qi::Vector # {RationalFunctionFieldElem{QQFieldElem,QQPolyRingElem}}
  root_system::RootSystem
  w0::Vector{UInt8}

  # cache
  canonical_basis::Dict{Vector{Int},PBWAlgElem}
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

function gen(U::QuantumGroup, i::Int)
  return QuantumGroupElem(U, U.gens[i])
end

function gens(U::QuantumGroup)
  return [gen(U, i) for i in 1:length(U.gens)]
end

mutable struct QuantumGroupElem
  U::QuantumGroup
  elem::PBWAlgElem
end

function Base.show(io::IO, x::QuantumGroupElem)
  return show(io, x.elem)
end

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

###############################################################################

function add!(z::QuantumGroupElem, x::QuantumGroupElem, y::QuantumGroupElem)
  z.elem = add!(z.elem, x.elem, y.elem)
  return z
end

function add!(x::QuantumGroupElem, y::QuantumGroupElem)
  return add!(x, x, y)
end

function addmul!(z::QuantumGroupElem, x::QuantumGroupElem, a)
  z.elem = addmul!(z.elem, x.elem, a)
  return z
end

function mul!(z::QuantumGroupElem, x::QuantumGroupElem, y::QuantumGroupElem)
  z.elem = mul!(z.elem, x.elem, y.elem)
  return z
end

function mul!(x::QuantumGroupElem, y::QuantumGroupElem)
  return mul!(x, x, y)
end

function mul!(z::QuantumGroupElem, x::QuantumGroupElem, a)
  z.elem = mul!(z.elem, x.elem, a)
  return z
end

function mul!(x::QuantumGroupElem, a)
  return mul!(x, x, a)
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

function submul!(z::QuantumGroupElem, x::QuantumGroupElem, a)
  z.elem = submul!(z.elem, x.elem, a)
  return z
end

function Base.:*(x::QuantumGroupElem, y::QuantumGroupElem)
  return mul!(deepcopy(x), y)
end

function Base.:*(x::QuantumGroupElem, a)
  return mul!(deepcopy(x), a)
end

function Base.:*(a::LaurentPolyWrap, x::QuantumGroupElem)
  return mul!(deepcopy(x), coefficient_ring(parent(x))(a))
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

function leading_coefficient(x::QuantumGroupElem)
  return leading_coefficient(x.elem)
end

function leading_exponent_vector(x::QuantumGroupElem)
  return leading_exponent_vector(x.elem)
end

function leading_monomial(x::QuantumGroupElem)
  return QuantumGroupElem(parent(x), leading_monomial(x.elem))
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
  GAP.Globals.SetCartanMatrix(gapR, GAP.Obj(cartan_matrix(R)))
  GAP.Globals.SetBilinearFormMat(gapR, gapBil)
  GAP.Globals.SetPositiveRootsNF(gapR, gapPos)
  GAP.Globals.SetSimpleSystemNF(gapR, gapSim)
  GAP.Globals.SetBilinearFormMatNF(gapR, gapBil)
  GAP.Globals.SetTypeOfRootSystem(
    gapR, collect(Iterators.flatmap(t -> GAP.GapObj.(t), root_system_type(R)))
  )

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

  alg, F = pbw_algebra(P, rels, deglex(theta))
  return QuantumGroup(
    A,
    alg,
    F,
    q,
    [q for _ in 1:length(F)],
    R,
    w0,
    Dict{Vector{Int},PBWAlgElem}(),
  )
end

function _terms(x::QuantumGroupElem)
  return terms(x.elem)
end

###############################################################################
#
#   Standard automorphisms
#
###############################################################################

function bar_involution(U::QuantumGroup)
  R = root_system(U)
  refl = weyl_group(R).refl

  nsim = number_of_simple_roots(R)
  npos = number_of_positive_roots(R)

  gen = zeros(Int, npos)
  for i in 1:npos
    beta = U.w0[i]
    for j in (i - 1):-1:1
      beta = refl[U.w0[j], beta]
    end
    gen[beta] = i
  end

  img = zeros(U.alg, length(U.gens))
  for i in 1:nsim
    img[gen[i]] = add!(img[gen[i]], U.gens[gen[i]])
  end

  # we construct the images of the PBW generators inductively by height
  for m in (nsim + 1):npos
    n = 0
    s = 0
    # find positive root with smaller height
    for i in 1:nsim
      n = Int(refl[i, m])
      if n < m
        s = i
        break
      end
    end

    if gen[s] > gen[n]
      s, n = n, s
    end

    b = one(U.alg)
    img[gen[m]] = add!(img[gen[m]], U.gens[gen[n]] * U.gens[gen[s]])
    for term in Singular.terms(U.alg.relations[gen[s], gen[n]])
      exp = Singular.leading_exponent_vector(term)
      if exp[gen[m]] != 0
        continue
      end

      for n in eachindex(exp)
        for _ in 1:exp[n]
          b = mul!(b, img[n])
        end
      end

      img[gen[m]] = submul!(
        img[gen[m]], b,
        evaluate(coefficient_ring(U.alg)(Singular.leading_coefficient(term)), U.q^-1),
      )
      b = one!(b)
    end
  end

  return function (x::QuantumGroupElem)
    @req parent(x) == U "parent mismatch"

    val = zero(U.alg)
    t = one(U.alg)
    for term in terms(x.elem)
      exp = leading_exponent_vector(term)
      for i in eachindex(exp)
        for _ in 1:exp[i]
          t = mul!(t, img[i])
        end
      end
      val = addmul!(val, t, evaluate(leading_coefficient(term), U.q^-1))
      t = one!(t)
    end

    return QuantumGroupElem(U, val)
  end
end

###############################################################################
#
#   Canonical basis
#
###############################################################################

function _canonical_basis_elem(U::QuantumGroup, b::Vector{Int})
  datum = Int.(b)
  bar = bar_involution(U)
  return get!(U.canonical_basis, datum) do
    F = one(U)
    for i in 1:length(datum)
      for _ in 1:datum[i]
        F = mul!(F, gen(U, i))
      end
      mul!(F, inv(quantum_factorial(U, datum[i])))
    end

    G = F
    F = sub!(bar(F), F)
    n = 0
    while !iszero(F)
      lc = numerator(leading_coefficient(F))
      c = coefficient_ring(U)(
        parent(lc.poly)(
          [zero(ZZ); (coeff(lc.poly, n) for n in (-lc.mindeg + 1):degree(lc.poly))...]
        ),
      )

      elem = canonical_basis_elem(U, leading_exponent_vector(F))
      G = addmul!(G, elem, c)
      F = submul!(F, elem, div(leading_coefficient(F), leading_coefficient(elem)))
      n += 1
      if n > 2
        error("infinite loop")
      end
    end

    return G.elem
  end
end

function canonical_basis_elem(U::QuantumGroup, b::Vector{Int})
  return QuantumGroupElem(U, _canonical_basis_elem(U, b))
end

function quantum_integer(U::QuantumGroup, n::Int)
  field = coefficient_ring(U.alg)
  z = zero(field)
  for i in (-n + 1):2:(n - 1)
    z = add!(z, field(U.q)^i)
  end
  return z
end

function quantum_factorial(U::QuantumGroup, n::Int)
  z = one(coefficient_ring(U.alg))
  for i in 1:n
    z = mul!(z, quantum_integer(U, i))
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

function set_canonical_basis_elem!(U::QuantumGroup, b::Vector{Int}, x::QuantumGroupElem)
  return U.canonical_basis[b] = x.elem
end

=#

###############################################################################
#
#   Internal functions
#
###############################################################################
