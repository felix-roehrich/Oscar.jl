struct Poset
  cov::Matrix{Int} # covering relations of the poset
  rel::BitVector # general relations
  set::BitVector # indicates whether a general relation was computed
  elems::Vector{Symbol} # symbols to use for the elements
end

struct PosetElem
  i::Int
  parent::Poset
end

function poset(cov::Matrix{Int}, elems::Vector{<:VarName}=["x_$i" for i in 1:ncols(cov)])
  @req is_upper_triangular(cov, 1) "matrix must be strictly upper triangular"
  @req nrows(cov) == ncols(cov) "must be a square matrix"
  @req ncols(cov) == length(elems) "size of matrix must match number of elements"

  d = nrows(cov)
  rel = BitVector(!iszero(cov[i, j]) for i in 1:d for j in (i + 1):d)
  return Poset(cov, rel, copy(rel), Symbol.(elems))
end

function poset_elem(P::Poset, i::Int)
  @req 1 <= i <= length(P.elems) "index out of range"
  return PosetElem(i, P)
end

function poset_elem(P::Poset, elem::VarName)
  i = findfirst(==(Symbol(elem)), P.elems)
  if isnothing(i)
    error("unknown element")
  end

  return poset_elem(P, i)
end

function Base.show(io::IO, P::Poset)
  print(io, "poset with $(length(P)) elements")
end

function Base.length(P::Poset)
  return length(P.elems)
end

function Base.:(<)(x::PosetElem, y::PosetElem)
  @req parent(x) === parent(y) "elements must belong to the same poset"
  ps = parent(x)

  # upper triangularity
  if y.i <= x.i
    return false
  end

  # linearised index
  len = ncols(ps.cov)
  xy = div((x.i - 1) * (2 * len - x.i), 2) + (y.i - x.i)

  # fast path
  if ps.set[xy]
    return ps.rel[xy]
  end

  # slow path using covering relations
  q = Int[x.i]

  while !isempty(q)
    @label outer
    n = last(q)

    # because of upper triangularity we only need to go to y.i
    for k in (n + 1):(y.i)
      if !iszero(ps.cov[n, k])
        # set the relation for all previous elements
        for m in q
          mk = div((m - 1) * (2 * len - m), 2) + (k - m)
          ps.set[mk] = true
          ps.rel[mk] = true
        end

        # we are done
        if k == y.i
          return ps.rel[xy]
        end

        ky = div((k - 1) * (2 * len - k), 2) + (y.i - k)
        if ps.set[ky]
          # check if can take the fast path
          if ps.rel[ky]
            # set relation for elements in the stack
            for m in q
              my = div((m - 1) * (2 * len - m), 2) + (y.i - m)
              ps.set[my] = true
              ps.rel[my] = true
            end
            return ps.rel[xy]
          else
            # k is not comparable to y
            continue
          end
        end

        # add k to the stack and continue from k
        push!(q, k)
        @goto outer
      end
    end

    # we now know that n is not comparable y
    ny = div((n - 1) * (2 * len - n), 2) + (y.i - n)
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

function (P::Poset)(i::Int)
  return poset_elem(P, i)
end

function (P::Poset)(elem::VarName)
  return poset_elem(P, elem)
end

function expressify(a::PosetElem; context=nothing)
  return parent(a).elems[a.i]
end
@enable_all_show_via_expressify PosetElem

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

function rel_index(P::Poset, x::Int, y::Int)
  return div((x - 1) * (2 * length(P) - x), 2) + (y - x)
end
