struct Poset
  cov::Matrix{Int} # covering relations of the poset
  rel::BitVector # general relations
  set::BitVector # indicates whether a general relation was computed
end

struct PosetElem
  i::Int
  parent::Poset
end

function Base.:(<)(x::PosetElem, y::PosetElem)
  @req parent(x) === parent(y) "elements must belong to the same poset"
  ps = parent(x)
  
  # upper triangularity
  if y.i < x.i
    return false
  end
  
  # linearised index
  len = ncols(ps.cov)
  xy = div((x.i-1)*(2*len-x.i), 2) + (y.i - x.i)
  
  # fast path
  if ps.set[xy]
    return ps.rel[xy]
  end
  
  # slow path using covering relations
  q = Int[x.i]
  while !isempty(q)
    n = last(q)
    
    # because of upper triangularity we only need to go to y.i
    for k in n+1:y.i
      if !iszero(ps.cov[n, k])
        # check if we are done or can take the fast path
        ky = div((k-1)*(2*len-k), 2) + (y.i - k)
        if k == y.i || ps.set[ky]
          # set relation for elements in the stack
          rel = k == y.i || ps.rel[ky]
          for m in q
            my = div((m-1)*(2*len-m), 2) + (y.i - m)
            ps.set[my] = true
            ps.rel[my] = rel
          end
          return ps.rel[xy]
        end
        
        # add k to the stack
        push!(q, k)
        continue
      end
    end
    
    # we now know that n is not comparable y
    ny = div((n-1)*(2*len-n), 2) + (y.i - n)
    ps.set[ny] = true
    ps.rel[ny] = false
    
    # continue with previous element
    pop!(q)
  end
  
  return false
end

function parent(x::PosetElem)
  return x.parent
end

function poset_elem(P::Poset, i::Int)
  return PosetElem(i, P)
end

function (P::Poset)(i::Int)
  return poset_elem(P, i)
end

function expressify(a::PosetElem; context=nothing)
  return a.i
end
@enable_all_show_via_expressify PosetElem

function poset(cov::Matrix{Int})
  @req is_upper_triangular(cov, 1) "matrix must be upper triangular"
  
  d = nrows(cov)
  rel = BitVector(!iszero(cov[i, j]) for i in 1:d for j in i+1:d)
  return Poset(cov, rel, copy(rel))
end

struct LSLattice
  bonds::Vector{Int}
  chain::Vector{Int}
  poset::Poset
end

struct LSLatticeElem
  parent::LSLattice
  point::Vector{Int}
end

function expressify(a::LSLatticeElem, s=:e; context=nothing)
  sum = Expr(:call, :+)
  for i in 1:length(a.point)
    push!(
      sum.args,
      Expr(:call, :*, expressify(a.point[i]; context=context), "$i"),
    )
  end
  return sum
end
@enable_all_show_via_expressify LSLatticeElem

struct LSMonoid
  degree::Int
  parent::LSLattice
  points::Vector{LSLatticeElem}
end

function ls_lattice(ps::Poset, chain::Vector{Int})
  # compute bonds
  b = zeros(Int, length(chain))
  b[1] = 1
  for i in 2:length(chain)
    b[i] = ps.covers[chain[i - 1], chain[i]]
  end

  return LSLattice(b, chain, ps)
end

function maximal_chains(cov::Matrix{Int})
  chains = Vector{Int}[]
  c = [1]
  j = 1
  while true
    j = findnext(!=(0), cov[last(c), :], j)
    if !isnothing(j)
      push!(c, j)
      continue
    end

    push!(chains, deepcopy(c))
    while !isempty(c)
      s = pop!(c) + 1
      if isempty(c)
        return chains
      end

      j = findnext(!=(0), cov[last(c), :], s)
      if !isnothing(j)
        break
      end
    end
  end
end

function (LS::LSLattice)(degree::Integer)
  l = lcm(LS.bonds)
  d = [div(l, b) for b in LS.bonds]

  pts = LSLatticeElem[]
  p = zero(LS.bonds)

  p[end] = l*degree
  while true
    push!(pts, LSLatticeElem(LS, deepcopy(p)))
  
    # search
    s = 0
    j = 2
    while p[j] < d[j]
      s += p[j-1]
      p[j-1] = 0
      j += 1
      if j > length(LS.chain)
        return LSMonoid(degree, LS, pts)
      end
    end

    # go to next point
    p[j] -= d[j]
    p[j - 1] += d[j]
  end
end

function points(LS::LSMonoid)
  return LS.points
end

#=
m = [0  1  1  0  0  0
0  0  0  1  1  0
0  0  0  1  1  0
0  0  0  0  0  1
0  0  0  0  0  1
0  0  0  0  0  0
]
=#
