struct QBosonAlgebra <: NCRing
  root_system::RootSystem
  pbw_algebra::PBWAlgebra{QuantumFieldElem}
  
  cvx::Vector{Int}
end

mutable struct QBosonAlgebraElem <: NCRingElem
  parent::QBosonAlgebra
  elem::PBWAlgebraElem{QuantumFieldElem}
end

function elem_type(::Type{QBosonAlgebra})
  return QBosonAlgebraElem
end

function parent_type(::Type{QBosonAlgebraElem})
  return QBosonAlgebra
end

function q_boson_algebra(R::RootSystem)
  w0 = Int.(word(longest_element(weyl_group(R))))
  inner_product = Oscar.bilinear_form(R)
  
  QF, q = quantum_field()
  vars = Symbol[]
  for i in 1:length(w0)
    push!(vars, Symbol("F$i"))
  end
  for i in 1:length(w0)
    push!(vars, Symbol("E$i'"))
  end

  P, theta = polynomial_ring(QF, vars)
  rels = Matrix{elem_type(P)}(undef, length(theta), length(theta))

  nsim = number_of_simple_roots(R)
  npos = number_of_positive_roots(R)

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
  gapBil = GAP.GapObj(inner_product)

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
  rels = Matrix{elem_type(P)}(undef, length(vars), length(vars))

  term = one(P)
  term2 = one(P)
  for i in 1:length(w0), j in i+1:length(w0)
    i2 = i + length(w0)
    j2 = j + length(w0)
    
    rep = Oscar.GAPWrap.ExtRepOfObj(gapF[j] * gapF[i])
    rels[i, j] = zero(P)
    rels[i2, j2] = zero(P)
    for n in 1:2:length(rep)
      for m in 1:2:length(rep[n])
        for _ in 1:rep[n][m + 1]
          term = mul!(term, theta[rep[n][m]])
          term2 = mul!(term2, theta[rep[n][m]+npos])
        end
        term = mul!(term, inv(q_factorial(rep[n][m + 1], q)))
        term2 = mul!(term2, inv(q_factorial(rep[n][m + 1], q)))
      end
      
      coeffRep = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[n + 1])
      coeff = zero(QF)
      for n in 1:length(coeffRep[1])
        coeff = addmul!(coeff, q^(n + coeffRep[2] - 1), coeffRep[1][n])
      end
      rels[i, j] = addmul!(rels[i, j], term, coeff)
      rels[i2, j2] = addmul!(rels[i2, j2], term2, coeff)
      term = one!(term)
      term2 = one!(term2)
    end
  end
  
  #=
  for i in length(w0)+1:2*length(w0), j in i+1:2*length(w0)
    rep = Oscar.GAPWrap.ExtRepOfObj(gapF[j+2*rank(R)] * gapF[i+2*rank(R)])
    rels[i, j] = zero(P)
    for n in 1:2:length(rep)
      for m in 1:2:length(rep[n])
        for _ in 1:rep[n][m + 1]
          term = mul!(term, theta[rep[n][m]-rank(R)])
        end
        term = mul!(term, inv(q_factorial(rep[n][m + 1], q)))
      end
      
      coeffRep = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[n + 1])
      coeff = zero(QF)
      for n in 1:length(coeffRep[1])
        coeff = addmul!(coeff, q^(n + coeffRep[2] - 1), coeffRep[1][n])
      end
      rels[i, j] = addmul!(rels[i, j], term, coeff)
      term = one!(term)
    end
  end
  =#
  
  cvx = zeros(Int, npos) # std -> cvx
  std = zeros(Int, npos) # cvx -> std
  for i in 1:npos
    beta = w0[i]
    for j in (i - 1):-1:1
      beta = weyl_group(R).refl[w0[j], beta]
    end
    cvx[beta] = i
    std[i] = beta
  end
  
  for i in 1:length(w0), j in length(w0)+1:2*length(w0)
    s = std[i]
    t = std[j-npos]
    if s <= nsim && t <= nsim
      e = Int(dot(coefficients(simple_root(R, s)) * inner_product, coefficients(simple_root(R, t))))
      rels[i, j] = q^-e * theta[i]*theta[j]
      if s == t
        rels[i, j] = add!(rels[i, j], one(P))
      end
    else
        rels[i, j] = zero(P)
    end
  end

  alg = pbw_algebra(P, rels)
  refl = weyl_group(R).refl
  for i in 1:npos, j in 1:npos
    if i <= nsim && j <= nsim
      continue
    end
    
    a = 0
    b = 0
    # find positive root with smaller height
    for k in 1:nsim
      a = Int(refl[k, i])
      if a < i
        b = k
        break
      end
    end
    
    c = 0
    d = 0
    for k in 1:nsim
      c = Int(refl[k, j])
      if c < j
        d = k
        break
      end
    end
    
    if i <= nsim
      c, d = cvx[c] < cvx[d] ? (cvx[c]+npos, cvx[d]+npos) : (cvx[d]+npos, cvx[c]+npos)
      
      cd = gen(alg, d) * gen(alg, c)
      cf = inv!(coeff(cd, length(cd)))
      cd = mul!(cd, cf)
      
      rel = zero(alg)
      rel = addmul!(rel, gen(alg, d) * (gen(alg, c) * gen(alg, cvx[i])), cf)
      rel = sub!(rel, (cd - gen(alg, cvx[j]+npos)) * gen(alg, cvx[i]))
      alg.mult[_linear_index(alg, cvx[i], cvx[j]+npos)][1, 1] = rel.poly
    elseif j <= nsim
      a, b = cvx[a] < cvx[b] ? (cvx[a], cvx[b]) : (cvx[b], cvx[a])
      
      ba = gen(alg, b) * gen(alg, a)
      cf = inv!(coeff(ba, length(ba)))
      ba = mul!(ba, cf)
      
      rel = zero(alg)
      rel = addmul!(rel, (gen(alg, cvx[j]+npos) * gen(alg, b)) * gen(alg, a), cf)
      rel = sub!(rel, gen(alg, cvx[j]+npos)*(ba - gen(alg, cvx[i])))

      alg.mult[_linear_index(alg, cvx[i], cvx[j]+npos)][1, 1] = rel.poly
    else
      a, b = cvx[a] < cvx[b] ? (cvx[a], cvx[b]) : (cvx[b], cvx[a])
      c, d = cvx[c] < cvx[d] ? (cvx[c]+npos, cvx[d]+npos) : (cvx[d]+npos, cvx[c]+npos)
      
      ba = gen(alg, b) * gen(alg, a)
      dc = gen(alg, d) * gen(alg, c)
      
      rel = zero(alg)
      rel = add!(rel, gen(alg, d) * (gen(alg, c) * gen(alg, b)) * gen(alg, a))
      rel = sub!(rel, gen(alg, d) * (gen(alg, c) * (ba - gen(alg, cvx[i]))))
      rel = sub!(rel, ((dc - gen(alg, cvx[j]+npos)) * gen(alg, b)) * gen(alg, a))
      rel = add!(rel, (dc - gen(alg, cvx[j]+npos)) * (ba - gen(alg, cvx[i])))
      
      alg.mult[_linear_index(alg, cvx[i], cvx[j]+npos)][1, 1] = rel.poly
    end
  end

  return QBosonAlgebra(R, alg, cvx)
end

function root_system(B::QBosonAlgebra)
  return B.root_system
end

function parent(x::QBosonAlgebraElem)
  return x.parent
end

function ngens(B::QBosonAlgebra)
  return ngens(B.pbw_algebra)
end

function Base.show(io::IO, x::QBosonAlgebraElem)
  print(io, x.elem)
end

function gen(B::QBosonAlgebra, i::Int)
  return QBosonAlgebraElem(B, gen(B.pbw_algebra, i))
end

function gens(B::QBosonAlgebra)
  return [gen(B, i) for i in 1:ngens(B)]
end

function f_gens(B::QBosonAlgebra)
  return [gen(B, i) for i in 1:div(ngens(B), 2)] # number_of_positive_roots(B.root_system)
end

function e_gens(B::QBosonAlgebra)
  R = root_system(B)
  return [gen(B, B.cvx[i]+number_of_positive_roots(R)) for i in 1:rank(R)] # number_of_positive_roots(B.root_system)
end

function one(B::QBosonAlgebra)
  return QBosonAlgebraElem(B, one(B.pbw_algebra))
end

function zero(B::QBosonAlgebra)
  return QBosonAlgebraElem(B, zero(B.pbw_algebra))
end

function zero(x::QBosonAlgebraElem)
  return zero(parent(x))
end

function monomial(B::QBosonAlgebra, i::Vector{Int}, n::Vector{Int})
  f = one(B)
  q = gen(coefficient_ring(B.pbw_algebra))
  
  for (i, n) in Iterators.zip(i, n)
    f = mul!(f, gen(B, B.cvx[i])^n)
    f = div!(f, q_factorial(n, q)) # B.qi[B.cvx[i]])
  end

  return f
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function add!(z::QBosonAlgebraElem, x::QBosonAlgebraElem, y::QBosonAlgebraElem)
  z.elem = add!(z.elem, x.elem, y.elem)
  return z
end

function div!(z::QBosonAlgebraElem, x::QBosonAlgebraElem, a::QuantumFieldElem)
  z.elem = div!(z.elem, x.elem, a)
  return z
end

function mul!(z::QBosonAlgebraElem, x::QBosonAlgebraElem, a::QuantumFieldElem)
  z.elem = mul!(z.elem, x.elem, a)
  return z
end

function mul!(z::QBosonAlgebraElem, x::QBosonAlgebraElem, y::QBosonAlgebraElem)
  z.elem = mul!(z.elem, x.elem, y.elem)
  return z
end

function neg!(z::QBosonAlgebraElem, x::QBosonAlgebraElem)
  z.elem = neg!(z.elem, x.elem)
  return z
end

function sub!(z::QBosonAlgebraElem, x::QBosonAlgebraElem, y::QBosonAlgebraElem)
  z.elem = sub!(z.elem, x.elem, y.elem)
  return z
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function Base.:+(x::QBosonAlgebraElem, y::QBosonAlgebraElem)
  check_parent(x, y)
  return add!(zero(x), x, y)
end

function Base.:-(x::QBosonAlgebraElem, y::QBosonAlgebraElem)
  check_parent(x, y)
  return sub!(zero(x), x, y)
end

function Base.:-(x::QBosonAlgebraElem)
  return neg!(zero(x), x)
end

function Base.:^(x::QBosonAlgebraElem, n::Int)
  return QBosonAlgebraElem(parent(x), x.elem^n)
end

function Base.:*(x::QBosonAlgebraElem, a::QuantumFieldElem)
  return mul!(zero(x), x, a)
end

function Base.:*(a::QuantumFieldElem, x::QBosonAlgebraElem)
  return mul!(zero(x), x, a)
end

function Base.:*(x::QBosonAlgebraElem, y::QBosonAlgebraElem)
  check_parent(x, y)
  return mul!(zero(x), x, y)
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

struct PBWAlgebraDerivation{T}
  domain::PBWAlgebra{T}
  codomain::PBWAlgebra{T}
  
  img::Vector{PBWAlgebraElem{T}}
  skew::Vector{T}
end

function (D::PBWAlgebraDerivation{T})(x::PBWAlgebraElem{T}) where {T}
  return image(D, x)
end

function image(D::PBWAlgebraDerivation{T}, x::PBWAlgebraElem{T}) where {T}
  @req parent(x) === D.domain "parent mismatch"
  return image!(zero(D.codomain), D, x)
end

function image!(z::PBWAlgebraElem{T}, D::PBWAlgebraDerivation{T}, x::PBWAlgebraElem{T}) where {T}
  z = zero!(z)
  m = one(D.codomain)
  tmp1 = zero(D.codomain)
  tmp2 = zero(D.codomain)

  exp = zeros(Int, ngens(D.domain))
  for i in 1:length(x)
    monomial_set!(m.poly, 1, x.poly, i)
    for j in 1:ngens(D.domain)
      if exponent(m, 1, j) > 0
        exp[j] = 1
        set_exponent!(m, 1, j, exponent(m, 1, j)-1)
        tmp1 = mul!(mul!(image!(tmp1, D, m), D.skew[j]), coeff(x, i))
        _mul_m_p!(D.domain, tmp2.poly, exp, tmp1.poly)
        exp[j] = 0
        
        z = add!(z, tmp2)
        z = addmul!(z, mul!(tmp1, D.img[j], m), coeff(x, i))
        break
      end
    end
  end

  return z
end

function q_boson_derivation(U::QuantumGroup, i::Int)
  R = root_system(U)
  nsim = number_of_simple_roots(R)
  npos = number_of_positive_roots(R)

  skew = Vector{QuantumFieldElem}(undef, npos)
  img = Vector{PBWAlgebraElem{QuantumFieldElem}}(undef, npos)
  
  for j in 1:nsim
    if i == j
      img[U.cvx[j]] = one(U.algebra)
    else
      img[U.cvx[j]] = zero(U.algebra)
    end
  end
  for j in 1:npos
    e = Int(dot(coefficients(simple_root(R, i)) * U.bilinear_form, coefficients(positive_root(R, j))))
    skew[U.cvx[j]] = gen(coefficient_ring(U))^-e
  end
  
  D = PBWAlgebraDerivation(U.algebra, U.algebra, img, skew)
  
  z = zero(U.algebra)
  m = one(U.algebra)
  for j in (nsim + 1):npos
    r = 0
    s = 0
    # find positive root with smaller height
    for k in 1:nsim
      r = Int(weyl_group(R).refl[k, j])
      if r < j
        s = k
        break
      end
    end
    
    e = Int(height(positive_root(R, j)) - height(positive_root(R, r)))
    if U.cvx[s] < U.cvx[r]
      rel = gen(U.algebra, U.cvx[r]) * gen(U.algebra, U.cvx[s])^e
      D.img[U.cvx[j]] = D.skew[U.cvx[r]] * gen(U.algebra, U.cvx[r]) * image(D, gen(U.algebra, U.cvx[s])^e) + D.img[U.cvx[r]] * gen(U.algebra, U.cvx[s])^e
    else
      rel = gen(U.algebra, U.cvx[s])^e * gen(U.algebra, U.cvx[r])
      D.img[U.cvx[j]] = D.skew[U.cvx[s]] * gen(U.algebra, U.cvx[s]) * image(D, gen(U.algebra, U.cvx[s])^(e-1) * gen(U.algebra, U.cvx[r])) + D.img[U.cvx[s]] * gen(U.algebra, U.cvx[s])^(e-1) * gen(U.algebra, U.cvx[r])
    end
    
    for i in 1:length(rel)-1
      setcoeff!(m, 1, coeff(rel, i))
      monomial_set!(m.poly, 1, rel.poly, i)
      D.img[U.cvx[j]] = sub!(D.img[U.cvx[j]], image!(z, D, m))
    end
    
    D.img[U.cvx[j]] = div!(D.img[U.cvx[j]], coeff(rel, length(rel)))
  end
  
  return D
end

function dot(U::QuantumGroup, i::Vector{Int}, n::Vector{Int})
  f = monomial(U, i, n)
  e = [Oscar.QuantumGroups.q_boson_derivation(U, j) for j in 1:rank(root_system(U))]
  
  z1 = f.elem
  z2 = zero(z1)
  for j in 1:length(i)
    println("$j")
    for _ in 1:n[j]
      z2 = Oscar.QuantumGroups.image!(z2, e[i[j]], z1)
      Oscar.QuantumGroups.swap!(z1, z2)
    end
  end
  
  cf = coeff(z1, length(z1))
  for j in 1:length(n)
    cf = div!(cf, q_factorial(n[j], U.qi[j]))
  end
  
  return cf
end
