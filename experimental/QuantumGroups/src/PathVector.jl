struct PathVector
  a::LSFanElem
end

function (p::PathVector)(a::LSFanElem)
  if length(max(a)) < 3
    if p.a == a
      return 1
    else
      return 0
    end
  end

  b = p.a
  LS = parent(p.a)
  W = weyl_group(LS)
  s = gens(W)
  seq = sequence(a)
  return binomial(trunc(Int, Float64(LS.dual_heighst_weight[1]*(b[s[1]*s[2]]+b[s[1]]))), trunc(Int, Float64(seq[3][2] - LS.dual_heighst_weight[1]*(b[s[1]*s[2]*s[1]]+b[s[2]*s[1]]))))
end
