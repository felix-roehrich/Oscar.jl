struct NZPolytope
  iota::Vector{Int}
  R::RootSystem
  
  function NZPolytope(R::RootSystem)
    iota = Int.(word(longest_element(weyl_group(R))))
    return new(iota, R)
  end
end

struct NZPolytopePoint
  coords::Vector{Int}
  parent::NZPolytope
end

function Base.show(io::IO, p::NZPolytopePoint)
  print(io, p.coords)
end

function Base.getindex(p::NZPolytopePoint, i::Int)
  return p.coords[i]
end

function Base.setindex!(p::NZPolytopePoint, n::Int, i::Int)
  return p.coords[i] = n
end

function Base.zero(P::NZPolytope)
  return NZPolytopePoint(zeros(Int, length(P.iota)), P)
end

function Base.parent(p::NZPolytopePoint)
  return p.parent
end

function root_system(P::NZPolytope)
  return P.R
end

function sigma(p::NZPolytopePoint, k::Int)
  P = parent(p)
  R = root_system(P)
  A = cartan_matrix(R)

  s = p[k]
  for i in k-1:-1:1
    s += A[P.iota[k], P.iota[i]] * p[i]
  end

  return s
end

function mi_sigma(p::NZPolytopePoint, i::Int)
  P = parent(p)
  k0 = findlast(==(i), P.iota)
  
  m = sigma(p, k0)
  mk = k0
  for k in k0-1:-1:1
    if P.iota[k] != i
      continue
    end
    m2 = max(m, sigma(p, k))
    if m < m2
      m = m2
      mk = k
    end
  end
  
  return mk
end

#=
function sigma0(p::NZPolytopePoint, k::Int)
  P = parent(P)
  R = root_system(P)
  A = cartan_matrix(R)

  x = -highest_weight(LS)[k]
  for i in 1:length(seq)
    x += A[k, Int(w0[i])] * seq[i]
  end

  return x
end
=#

function fi!(p::NZPolytopePoint, i::Int)
  mk = mi_sigma(p, i)
  p[mk] += 1
  return p
end