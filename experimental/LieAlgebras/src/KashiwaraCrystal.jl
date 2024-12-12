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

Return $\epsilon_i(b)$.
"""
function Crystal.eps(::AbstractCrystalElem, ::Int)
  error("not implemented")
end

@doc raw"""
    Crystal.phi(b::AbstractCrystalElem, i::Int) -> Int

Return $\varphi_i(b)$.
"""
function Crystal.phi(::AbstractCrystalElem, ::Int)
  error("not implemented")
end

@doc raw"""
    weight(b::AbstractCrystalElem) -> WeightLatticeElem

Return the weight of `b`.
"""
function Crystal.weight(::AbstractCrystalElem)
  error("not implemented")
end

# AbstractCrystalElem provided methods

# e_i variants

@doc raw"""
    Crystal.e(b::AbstractCrystalElem, i::Int) -> AbstractCrystalElem

Return the result of applying $\tilde e_i$ to `b`.
"""
function Crystal.e(b::AbstractCrystalElem, i::Int)
  return Crystal.e!(deepcopy(b), i)
end

@doc raw"""
    Crystal.e(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem

Shorthand for `n`-fold appliction of `Crystal.e(b, i)`.
"""
function Crystal.e(b::AbstractCrystalElem, i::Int, n::Int)
  return Crystal.e!(deepcopy(b), i, n)
end

@doc raw"""
    Crystal.e!(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem
    
Shorthand for `n`fold appliction of `Crystal.e!(b, i)`.
"""
function Crystal.e!(b::AbstractCrystalElem, i::Int, n::Int)
  for _ in 1:n
    Crystal.e!(b, i)
  end
  return b
end

@doc raw"""
    Crystal.e!(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem
    
Shorthand for `n`fold appliction of `Crystal.e!(b, i)`.
"""
function Crystal.e(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  Crystal.e!(deepcopy(b), i, n)
end

@doc raw"""
    Crystal.e!(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int}) -> AbstractCrystalElem

Return the result of $\tilde e_{i_1} \dots \tilde e_{i_r}b$.
"""
function Crystal.e!(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  @req length(i) == length(n) "length mismatch"

  for l in length(i):-1:1
    Crystal.e!(b, i[l], n[l])
  end
  return b
end

# f_i variants

@doc raw"""
    Crystal.f(b::AbstractCrystalElem, i::Int) -> AbstractCrystalElem

Return the result of applying $\tilde f_i$ to `b`.
"""
function Crystal.f(b::AbstractCrystalElem, i::Int)
  return Crystal.f!(deepcopy(b), i)
end

@doc raw"""
    Crystal.e(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem

Shorthand for `n`fold appliction of `Crystal.f(b, i)`.
"""
function Crystal.f(b::AbstractCrystalElem, i::Int, n::Int)
  return Crystal.f!(deepcopy(b), i, n)
end

@doc raw"""
    Crystal.f!(b::AbstractCrystalElem, i::Int, n::Int) -> AbstractCrystalElem
    
Shorthand for `n`fold appliction of `Crystal.f!(b, i)`.
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

Return the result of $\tilde e_{i_1} \dots \tilde e_{i_r}b$.
"""
function Crystal.f!(b::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
  @req length(i) == length(n) "i and n must have same length"

  for l in length(i):-1:1
    Crystal.f!(b, i[l], n[l])
  end
  return b
end
