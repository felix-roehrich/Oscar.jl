module Crystal

function e end

function e! end

function f end

function f! end

function eps end

function phi end

function weight end

end

function Base.deepcopy_internal(::AbstractCrystalElem)
  error("not implemented")
end

function Crystal.e!(::AbstractCrystalElem, ::Int)
  error("not implemented")
end

function Crystal.f!(::AbstractCrystalElem, ::Int)
  error("not implemented")
end

@doc raw"""
    Crystal.eps(b::AbstractCrystalElem, i::Int) -> Int
"""
function Crystal.eps(::AbstractCrystalElem, ::Int)
  error("not implemented")
end

@doc raw"""
    Crystal.phi(b::AbstractCrystalElem, i::Int) -> Int
"""
function Crystal.phi(::AbstractCrystalElem, ::Int)
  error("not implemented")
end

@doc raw"""
    weight(b::AbstractCrystalElem) -> WeightLatticeElem
"""
function Crystal.weight(::AbstractCrystalElem)
  error("not implemented")
end

# AbstractCrystalElem provided methods

# e_i variants

@doc raw"""
    Crystal.e(b::AbstractCrystalElem, i::Int) -> AbstractCrystalElem
"""
function Crystal.e(b::AbstractCrystalElem, i::Int)
  return Crystal.e!(deepcopy(b), i)
end

@doc raw"""
    Crystal.e(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem
"""
function Crystal.e(b::AbstractCrystalElem, i::Int, n::Int)
  return Crystal.e!(deepcopy(b), i, n)
end

@doc raw"""
    Crystal.e!(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem
"""
function Crystal.e!(b::AbstractCrystalElem, i::Int, n::Int)
  for _ in 1:n
    Crystal.e!(b, i)
  end
  return b
end

@doc raw"""
    Crystal.e!(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem
"""
function Crystal.e(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  Crystal.e!(deepcopy(b), i, n)
end

@doc raw"""
    Crystal.e!(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int}) -> AbstractCrystalElem
"""
function Crystal.e!(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  @req length(i) == length(n) "length mismatch"

  for l in 1:length(i)
    Crystal.e!(b, i[l], n[l])
  end
  return b
end

# f_i variants

@doc raw"""
    Crystal.f(b::AbstractCrystalElem, i::Int) -> AbstractCrystalElem
"""
function Crystal.f(b::AbstractCrystalElem, i::Int)
  return Crystal.f!(deepcopy(b), i)
end

@doc raw"""
    Crystal.e(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem
"""
function Crystal.f(b::AbstractCrystalElem, i::Int, n::Int)
  return Crystal.f!(deepcopy(b), i, n)
end

@doc raw"""
    Crystal.f!(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem
"""
function Crystal.f!(b::AbstractCrystalElem, i::Int, n::Int)
  for _ in 1:n
    Crystal.f!(b, i)
  end
  return b
end

@doc raw"""
    Crystal.f(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int}) -> AbstractCrystalElem
"""
function Crystal.f(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  Crystal.f!(deepcopy(b), i, n)
end

@doc raw"""
    Crystal.f!(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int}) -> AbstractCrystalElem
"""
function Crystal.f!(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  @req length(i) == length(n) "i and n must have same length"

  for l in 1:length(i)
    Crystal.f!(b, i[l], n[l])
  end
  return b
end

struct CrystalTensorProduct <: AbstractCrystal
end

struct CrystalTensorProductElem <: AbstractCrystalElem
  b::Vector{<:AbstractCrystalElem}
  parent::CrystalTensorProduct
end

function Base.show(io::IO, b::CrystalTensorProductElem)
  print(io, "$(b.b[1]) (x) $(b.b[2])")
end

function Crystal.f!(b::CrystalTensorProductElem, i::Int)
  if Crystal.phi(b.b[1], i) > Crystal.eps(b.b[2], i)
    Crystal.f!(b.b[1], i)
  else
    Crystal.f!(b.b[2], i)
  end
  return b
end

function Crystal.e!(b::CrystalTensorProductElem, i::Int)
  if Crystal.phi(b.b[1], i) >= Crystal.eps(b.b[2], i)
    Crystal.e!(b.b[1], i)
  else
    Crystal.e!(b.b[2], i)
  end
end

function Crystal.eps(b::CrystalTensorProductElem, i::Int)
  return max(
    Crystal.eps(b.b[1], i), Crystal.eps(b.b[2], i) - Int(Crystal.weight(b.b[1])[i])
  )
end

function tensor_product(b1::AbstractCrystalElem, b2::AbstractCrystalElem)
  return CrystalTensorProductElem([b1, b2], CrystalTensorProduct())
end
