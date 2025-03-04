function bar_involution(U::QuantumGroup)
  R = root_system(U)
  refl = weyl_group(R).refl

  nsim = number_of_simple_roots(R)
  npos = number_of_positive_roots(R)

  # order of the positive roots in the convex order
  cvx = zeros(Int, npos)
  for i in 1:npos
    beta = U.w0[i]
    for j in (i-1):-1:1
      beta = refl[U.w0[j], beta]
    end
    cvx[beta] = i
  end

  img = zeros(U.alg, ngens(U.alg))
  for i in 1:nsim
    img[cvx[i]] = add!(img[cvx[i]], gen(U.alg, cvx[i]))
  end

  # we construct the images of the PBW generators inductively by height
  for m in (nsim+1):npos
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
        coeff = inv(leading_coefficient(term))
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
    for i in 1:length(x)
      exp = exponent_vector(x, i)
      for j in length(exp):-1:1
        t = mul!(t, gen(U, j)^exp[j])
      end
      set_coeff!(t, coeff(x, i))
      z = add!(z, t)
    end

    return z
  end
end
