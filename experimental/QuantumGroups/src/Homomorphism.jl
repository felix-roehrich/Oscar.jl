function (f::QuantumGroupHomomorphism)(x::QuantumGroupElem)
  z = zero(parent(x))
  for i in 1:length(x)
    exp = exponent_vector(x, i)
    for j in 1:length(exp)
      for _ in 1:exp[j]
        t = mul!(t, f.img[i])
      end
    end
  end

  return z
end

function bar_involution(U::QuantumGroup)
  R = root_system(U)
  refl = weyl_group(R).refl

  nsim = number_of_simple_roots(R)
  npos = number_of_positive_roots(R)

  # order of the positive roots in the convex order
  cvx = zeros(Int, npos)
  for i in 1:npos
    beta = U.w0[i]
    for j in (i - 1):-1:1
      beta = refl[U.w0[j], beta]
    end
    cvx[beta] = i
  end

  img = zeros(U.alg, ngens(U.alg))
  for i in 1:nsim
    img[cvx[i]] = add!(img[cvx[i]], gen(U.alg, cvx[i]))
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

    pow = Int(height(positive_root(R, m)) - height(positive_root(R, n)))
    if cvx[s] < cvx[n]
      rel = gen(U.alg, cvx[n]) * gen(U.alg, cvx[s])^pow
      img[cvx[m]] = img[cvx[n]] * img[cvx[s]]^pow
    else
      rel = gen(U.alg, cvx[s])^pow * gen(U.alg, cvx[n])
      img[cvx[m]] = img[cvx[s]]^pow * img[cvx[n]]
    end

    b = one(U.alg)
    barred = zero(coefficient_ring(U))
    coeff = zero(coefficient_ring(U))
    for term in terms(rel)
      exp = leading_exponent_vector(term)
      if exp[cvx[m]] != 0
        coeff = bar!(coeff, inv((leading_coefficient(term))))
        continue
      end

      for n in eachindex(exp)
        for _ in 1:exp[n]
          b = mul!(b, img[n])
        end
      end

      bar!(barred, leading_coefficient(term))
      img[cvx[m]] = submul!(img[cvx[m]], b, barred)
      b = one!(b)
    end
    img[cvx[m]] = mul!(img[cvx[m]], coeff)
  end

  return function (x::QuantumGroupElem)
    @req parent(x) == U "parent mismatch"

    val = zero(U.alg)
    t = one(U.alg)
    barred = zero(coefficient_ring(U))
    for term in terms(x.elem)
      exp = leading_exponent_vector(term)
      for i in eachindex(exp)
        for _ in 1:exp[i]
          t = mul!(t, img[i])
        end
      end

      bar!(barred, leading_coefficient(term))
      val = addmul!(val, t, barred)
      t = one!(t)
    end

    return QuantumGroupElem(U, val)
  end
end

function star_involution(U::QuantumGroup)
  return function (x::QuantumGroupElem)
    z = zero(U)
    t = one(U)

    i = 1
    for exp in Singular.exponent_vectors(x.elem)
      for j in length(exp):-1:1
        t = mul!(t, gen(U, j)^exp[j])
      end
      setcoeff!(t, 1, coeff(x.elem, i))
      z = add!(z, t)
      i += 1
    end

    return z
  end
end

function braid_automorphism(U::QuantumGroup, i::Int)
  F = chevalley_gens(U)
  img0 = QuantumGroupElem[]
  for j in 1:(i - 1)
    push!(img0, F[j] * F[i] - U.qi[i] * F[i] * F[j])
  end
  push!(img0, F[i])
  for j in (i + 1):rank(root_system(U))
    push!(img0, F[j] * F[i] - U.qi[i] * F[i] * F[j])
  end

  img = _image(U, img0)
  return QuantumGroupHomomorphism(img)
end

function frobenius_splitting(U::QuantumGroup, l::Int)
  img = _image(U, [f^l for f in chevalley_gens(U)])
  return QuantumGroupHomomorphism(img)
end

###############################################################################
#
#   Internals
#
###############################################################################

# Accepts as input the image of Chevalley generators and returns the image of
# the PBW generators. Treated as QQ(q)-linear map.
function _image(U::QuantumGroup, image::Vector{QuantumGroupElem})
  R = root_system(U)
  cvx = U.cvx
  refl = weyl_group(R).refl

  nsim = number_of_simple_roots(R)
  npos = number_of_positive_roots(R)

  img = zeros(U.alg, ngens(U.alg))
  for i in 1:nsim
    img[U.cvx[i]] = add!(img[U.cvx[i]], image[i].elem)
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

    pow = Int(height(positive_root(R, m)) - height(positive_root(R, n)))
    if cvx[s] < cvx[n]
      rel = gen(U.alg, cvx[n]) * gen(U.alg, cvx[s])^pow
      img[cvx[m]] = img[cvx[n]] * img[cvx[s]]^pow
    else
      rel = gen(U.alg, cvx[s])^pow * gen(U.alg, cvx[n])
      img[cvx[m]] = img[cvx[s]]^pow * img[cvx[n]]
    end

    b = one(U.alg)
    coeff = zero(coefficient_ring(U))
    for term in terms(rel)
      exp = leading_exponent_vector(term)
      if exp[cvx[m]] != 0
        coeff = inv(leading_coefficient(term))
        continue
      end

      for n in eachindex(exp)
        for _ in 1:exp[n]
          b = mul!(b, img[n])
        end
      end

      img[cvx[m]] = submul!(img[cvx[m]], b, leading_coefficient(term))
      b = one!(b)
    end
    img[cvx[m]] = mul!(img[cvx[m]], coeff)
  end
end
