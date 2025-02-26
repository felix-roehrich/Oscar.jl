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

struct QuantumGroupElem
  U::QuantumGroup
  elem::PBWAlgElem
end

function parent(x::QuantumGroupElem)
  return x.U
end

function quantum_group(R::RootSystem, w0=word(longest_element(weyl_group(R))))
  A, q = laurent_polynomial_ring(ZZ, "q")
  QA = fraction_field(A)
  P, theta = polynomial_ring(QA, :F => 1:length(w0))

  # compute the index of the positive root from w0 in R.positive_roots
  refl = weyl_group(R).refl
  w0_pos = zeros(Int, length(w0))
  pos_w0 = zeros(Int, length(w0))
  for i in 1:length(w0)
    beta = w0[i]
    for j in (i - 1):-1:1
      beta = refl[w0[j], beta]
    end
    w0_pos[i] = beta
    pos_w0[beta] = i
  end

  rels = zero_matrix(P, length(w0), length(w0))
  root = zero(RootSpaceElem, R)

  nsim = number_of_simple_roots(R)
  npos = number_of_positive_roots(R)

  for i in 1:npos, j in (i + 1):npos
    if pos_w0[i] > pos_w0[j]
      i, j = j, i
    end
    ij = Int(
      dot(
        coefficients(positive_root(R, i)) * cartan_matrix(R),
        coefficients(positive_root(R, j)),
      ),
    )

    if ij >= 0
      k = l = 0
      for n in 1:i
        root = sub!(root, positive_root(R, j), positive_root(R, n))
        t, l = is_positive_root_with_index(root)
        if t
          k = n
          break
        end
      end

      # l == 0 means that this a commutation relation
      if l == 0
        rels[pos_w0[i], pos_w0[j]] = theta[pos_w0[i]] * theta[pos_w0[j]]
        continue
      end

      # compute new relation from lower relations
      rel = one(P)
      if pos_w0[i] < pos_w0[k]
        rel = mul!(rel, rels[pos_w0[i], pos_w0[k]])
      else
        rel = mul!(rel, theta[pos_w0[k]] * theta[pos_w0[i]])
      end
      if pos_w0[i] < pos_w0[l]
        rel = mul!(rel, rels[pos_w0[i], pos_w0[l]])
      else
        rel = mul!(rel, theta[pos_w0[l]] * theta[pos_w0[i]])
      end
      println(rel)

      s = one(P)
      for t in terms(rels[pos_w0[k], pos_w0[l]])
        exp = leading_exponent_vector(t)
        if !iszero(exp[pos_w0[j]])
          continue
        end

        for n in eachindex(exp)
          for _ in 1:exp[n]
            if pos_w0[i] < n
              s = mul!(s, rels[pos_w0[i], n])
            else
              s = mul!(s, theta[n] * theta[pos_w0[i]])
            end
          end
        end
        rel = sub!(rel, s)
        s = one!(s)
      end
      println("ij = $i$j")
      println(rel)
      rel = add!(rel, theta[pos_w0[i]] * theta[pos_w0[j]])
      rels[pos_w0[i], pos_w0[j]] = rel
    elseif ij == -1
      root = add!(vec, positive_root(R, i), positive_root(R, j))
      _, l = is_positive_root_with_index(root)
      rels[pos_w0[i], pos_w0[j]] =
        theta[pos_w0[l]] + q^-ij * theta[pos_w0[i]] * theta[pos_w0[j]]
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
      val = addmul!(val, t, leading_coefficient(term))
    end

    return QuantumGroupElem(U, val)
  end
end

#=
function quantum_integer(U::QuantumGroup, n::Int)
  z = zero(U.q)
  for i in (-n + 1):2:(n - 1)
    z = add!(z, U.q^i)
  end
  return z
end

function quantum_factorial(U::QuantumGroup, n::Int)
  z = one(U.q)
  for i in 1:n
    z = mul!(z, quantum_integer(U, i))
  end
  return z
end

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

function canonical_basis_elem(U::QuantumGroup, b::Vector{Int})
  datum = Int.(b)
  return get!(U.canonical_basis, datum) do
    F = one(U.alg)
    bar = bar_involution(U)
    for i in 1:length(datum)
      for _ in 1:datum[i]
        F = mul!(F, U.gens[i])
      end
    end

    G = deepcopy(F)
    F = sub!(F, bar(F))
    while !iszero(F)
      lc = leading_coefficient(F)
      c = parent(lc.poly)(
        [zero(ZZ); (coeff(lc.poly, n) for n in (-lc.mindeg + 1):degree(lc.poly))...]
      )
      elem = canonical_basis_elem(U)

      G = addmul!(G, elem, c)
      F = sub!(F, elem)
    end

    return F
  end
end
=#
