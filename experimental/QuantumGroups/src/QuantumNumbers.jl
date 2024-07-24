#A, v = LaurentPolynomialRing(ZZ, :v)

function q_binomial(a::IntegerUnion, t::IntegerUnion, q)
  @req t >= 0 "$t must be a positive integer"
  return evaluate(prod(v^(a-s) - v^(-a+s) for s in 0:t-1; init=A(1)) / prod(v^s - v^-s for s in 1:t; init=A(1)), q)
end

function check(v::Vector{LSFanElem})
  for a in v
    i, n = first(sequence(a))
    if !isnothing(ei!(deepcopy(a), i, n+1))
      return a
    end
  end
  return nothing
end

function verify()
  A, v = cyclotomic_field(24)
  b = [2,2,2,2,2,2,2,2]
  qbinom = (a, t) -> prod(v^(a-s) - v^(-a+s) for s in 0:t-1; init=A(1)) / prod(v^s - v^-s for s in 1:t; init=A(1))

  s = zero(A)
  r = [2,2,0,0,0,0,0,0]
  n = 0
  while !isnothing(r)
    q = [2+r[1],2+r[2],2+r[3],2+r[4],2+r[5],0,0,0]
    qt = 0
    for i in 1:length(q)
      qt += q[i]
      if qt >= 8
        q[i] -= qt - 8
        q[i+1:end] .= 0
        break
      end
    end
    
    while !isnothing(q)
      n += 1
      s += prod(v^(-(12-sum(r[t] for t in i:8))*(2-r[i]))*
      v^(-(24-sum(q[t] for t in i:8))*(2+r[i]-q[i]))*
      v^(-(12-sum(2-r[t] for t in i:8))*(q[i]-r[i]))*
      v^(-(12-sum(4-q[t] for t in i:8))*(-q[i]))*
      qbinom(2, r[i])*qbinom(2+r[i], q[i]) for i in 1:8)
      q = next_partition(q, [2+r[1],2+r[2],2+r[3],2+r[4],2+r[5],0,0,0])
    end
    r = next_partition(r, b)
  end
  return (n, s)
end
