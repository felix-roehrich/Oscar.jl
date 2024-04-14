struct LSLattice
  bonds::Vector{Int}
  chain::Vector{Int}
  poset::Poset
end

struct LSLatticeElem
  parent::LSLattice
  point::Vector{QQFieldElem}
  support::Vector{Int} # rename to chain
end

function Base.:(==)(a::LSLatticeElem, b::LSLatticeElem)
  return support(a) == support(b) && a.point == b.point
end

function Base.hash(a::LSLatticeElem, h::UInt)
  b = 0x4bba5b0a76b3db59 % UInt
  h = hash(a.point, h)
  h = hash(a.support, h)

  return xor(h, b)
end

function parent(a::LSLatticeElem)
  return a.parent
end

function support(a::LSLatticeElem)
  return a.support
end

struct LSMonoid
  degree::Int
  parent::LSLattice
  points::Vector{LSLatticeElem}
end

struct LSFanOfMonoids
  degree::Int
  points::Set{LSLatticeElem}
  poset::Poset
end

function Base.show(io::IO, LS::LSFanOfMonoids)
  io = pretty(io)
  print(io, "LS fan of monoids of degree $(LS.degree) with $(length(LS.points)) elements")
end

function expressify(a::LSLatticeElem, s=:e; context=nothing)
  sum = Expr(:call, :+)
  for i in 1:length(a.point)
    p = poset_elem(parent(a).poset, a.support[i])
    push!(sum.args, Expr(:call, :*, expressify(a.point[i]; context=context), "$s($p)"))
  end
  return sum
end
@enable_all_show_via_expressify LSLatticeElem

function ls_lattice(ps::Poset, chain::Vector{Int})
  # compute bonds
  b = zeros(Int, length(chain))
  b[1] = 1
  for i in 2:length(chain)
    b[i] = ps.cov[chain[i - 1], chain[i]]
  end

  return LSLattice(b, chain, ps)
end

function (LS::LSLattice)(degree::Integer)
  l = lcm(LS.bonds)
  d = [div(l, b) for b in LS.bonds]

  pts = LSLatticeElem[]
  p = zero(LS.bonds)

  p[end] = l * degree
  while true
    push!(
      pts,
      LSLatticeElem(
        LS,
        [QQ(x, l) for x in p if !iszero(x)],
        [LS.chain[i] for i in 1:length(p) if !iszero(p[i])],
      ),
    )

    # search
    s = 0
    j = 2
    while p[j] < d[j]
      s += p[j - 1]
      p[j - 1] = 0
      j += 1
      if j > length(LS.chain)
        return LSMonoid(degree, LS, pts)
      end
    end

    # go to next point
    p[j] -= d[j]
    p[j - 1] += s + d[j]
  end
end

function points(LS::LSMonoid)
  return LS.points
end

function points(LS::LSFanOfMonoids)
  return LS.points
end

function ls_fan_of_monoids(P::Poset, degree::Integer)
  fan = LSFanOfMonoids(degree, Set{LSLatticeElem}(), P)

  chains = maximal_chains(P.cov)
  for chain in chains
    monoid = ls_lattice(P, chain)(degree)
    union!(fan.points, monoid.points)
  end

  return fan
end
