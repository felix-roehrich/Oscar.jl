struct MPolyRing{T}
  N::Int
  ordering::Matrix{Int}
end

function Base.one(R::MPolyRing{T}) where {T}
  return MPolyRingElem(R, 3, zeros(Int, cld(ngens(R), sizeof(Int)), 1), [1], 1)
end

function Base.zero(R::MPolyRing)
  return MPolyRingElem(R)
end

function gen(R::MPolyRing{Int}, i::Int)
  x = one(R)
  monomial_set!(x, 1, i, 1)
  return x
end

function ngens(R::MPolyRing)
  return R.N
end

mutable struct MPolyRingElem{T}
  parent::MPolyRing{T}

  bit::Int
  exps::Matrix{Int}
  coeffs::Vector{T}
  len::Int

  function MPolyRingElem(R::MPolyRing{T}) where {T}
    return new{T}(R)
  end

  function MPolyRingElem(
    R::MPolyRing{T}, bit::Int, exps::Matrix{Int}, coeffs::Vector{T}, len::Int
  ) where {T}
    return new{T}(R, bit, exps, coeffs, len)
  end
end

function length(x::MPolyRingElem)
  return x.len
end

function parent(x::MPolyRingElem{T}) where {T}
  return x.parent
end

function zero(x::MPolyRingElem)
  return zero(parent(x))
end

function zero!(x::MPolyRingElem)
  x.len = 0
  return x
end

function swap!(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  x.bit, y.bit = y.bit, x.bit
  x.exps, y.exps = y.exps, x.exps
  x.coeffs, y.coeffs = y.coeffs, x.coeffs
  x.len, y.len = y.len, x.len
end

function coeff(x::MPolyRingElem, i::Int)
  return x.coeffs[i]
end

function ngens(x::MPolyRingElem)
  return ngens(parent(x))
end

function Base.:+(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  check_parent(x, y)
  return add!(zero(x), x, y)
end

function monomial_cmp(R::MPolyRing, m::Vector{Int}, n::Vector{Int})
  cmp = R.ordering * m - R.ordering * n
  for c in cmp
    if c != 0
      return c
    end
  end
end

function monomial_set!(x::MPolyRingElem, i::Int, y::MPolyRingElem, j::Int)
  for k in 1:size(x.exps, 1)
    x.exps[k, i] = y.exps[k, j]
  end
end

function monomial_set!(x::MPolyRingElem, i::Int, k::Int, e::Int)
  bits = 1 << x.bit
  j = k >> x.bit + 1
  @inbounds x.exps[j, i] =
    (x.exps[j, i] & ~((1 << bits - 1) << ((k - 1) * bits))) | (e << ((k - 1) * bits))
end

function unpack_exponent(p::Int, i::Int, bits::Int)
  return (p >> ((i - 1) * bits)) & (1 << bits - 1)
end

function pack_exponent!(p::Ptr{Int}, k::Int, e::Int, bits::Int)
  unsafe_modify!(p, &, ~((1 << bits - 1) << ((k - 1) * bits)))
  unsafe_modify!(p, |, e << ((k - 1) * bits))
end

function unpack_monomial!(m::Vector{Int}, x::MPolyRingElem, i::Int)
  bits = 1 << x.bit
  for k in 1:ngens(x)
    @inbounds m[k] = unpack_exponent(x.exps[k >> x.bit + 1, i], k, bits)
  end
end

function pack_monomial!(x::MPolyRingElem, i::Int, m::Vector{Int})
  bits = 1 << x.bit
  p = pointer(x.exps) + (i - 1) * sizeof(Int)
  for k in 1:ngens(x)
    @inbounds pack_exponent!(p + (k >> x.bit) * sizeof(Int), k, m[k], bits)
  end
end

function fit!(x::MPolyRingElem, len::Int, bit::Int)
  resize!(x.coeffs, len)
  if bit > x.bit || len > size(x.exps, 2)
    bits = 1 << bit
    xbits = 1 << x.bit
    exps = Matrix{Int}(
      undef, cld(ngens(x), (sizeof(Int) * 8) >> bit), max(len, 2 * size(x.exps, 2))
    )

    for i in 1:size(x.exps, 2), k in 1:ngens(x)
      @inbounds pack_exponent!(
        pointer(exps, k >> bit + i),
        k,
        unpack_exponent(x.exps[k >> x.bit + 1, i], k, xbits),
        bits,
      )
    end
    x.exps = exps
    x.bit = bit
  end
end

function add!(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  R = parent(x)

  len = length(x) + length(y)
  fit!(x, len, max(x.bit, y.bit))

  for i in 1:length(x)
    x.coeffs[length(y) + i] = x.coeffs[i]
    monomial_set!(x, length(y) + i, x, i)
  end

  i = length(y) + 1
  j = k = 1

  mx = Vector{Int}(undef, ngens(R))
  unpack_monomial!(mx, x, i)

  my = Vector{Int}(undef, ngens(R))
  unpack_monomial!(my, y, j)

  while i <= len && j <= length(y)
    cmp = monomial_cmp(R, mx, my)
    if cmp > 0
      x.coeffs[k] = x.coeffs[i]
      pack_monomial!(x, k, mx)
      i += 1
      unpack_monomial!(mx, x, i)
    elseif cmp == 0
      x.coeffs[k] = add!(x.coeffs[i], y.coeffs[j])
      if !iszero(x.coeffs[k])
        pack_monomial!(x, k, i)
      else
        k -= 1
      end

      i += 1
      j += 1
      unpack_monomial!(mx, x, i)
      unpack_monomial!(my, y, j)
    else
      x.coeffs[k] = deepcopy(y.coeffs[j])
      pack_monomial!(x, k, j)
      j += 1
      unpack_monomial!(my, y, j)
    end
    k += 1
  end
  while i <= len
    x.coeffs[k] = x.coeffs[i]
    monomial_set!(x, k, x, i)
    i += 1
    k += 1
  end
  while j <= length(y)
    x.coeffs[k] = deepcopy(y.coeffs[j])
    monomial_set!(x, k, y, j)
    j += 1
    k += 1
  end

  x.len = k - 1
  resize!(x.coeffs, x.len)
  return x
end

function add!(z::MPolyRingElem{T}, x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  R = parent(x)

  len = length(x) + length(y)
  fit!(z, len, max(x.bit, y.bit))

  i = j = k = 1

  mx = Vector{Int}(undef, ngens(R))
  unpack_monomial!(mx, x, i)

  my = Vector{Int}(undef, ngens(R))
  unpack_monomial!(my, y, j)

  while i <= length(x) && j <= length(y)
    cmp = monomial_cmp(R, mx, my)
    if cmp > 0
      z.coeffs[k] = x.coeffs[i]
      pack_monomial!(x, k, mx)
      i += 1
      unpack_monomial!(mx, x, i)
    elseif cmp == 0
      x.coeffs[k] = add!(x.coeffs[i], y.coeffs[j])
      if !iszero(x.coeffs[k])
        pack_monomial!(x, k, i)
      else
        k -= 1
      end

      i += 1
      j += 1
      unpack_monomial!(mx, x, i)
      unpack_monomial!(my, y, j)
    else
      x.coeffs[k] = deepcopy(y.coeffs[j])
      pack_monomial!(x, k, j)
      j += 1
      unpack_monomial!(my, y, j)
    end
    k += 1
  end
  while i <= len
    x.coeffs[k] = x.coeffs[i]
    monomial_set!(x, k, x, i)
    i += 1
    k += 1
  end
  while j <= length(y)
    x.coeffs[k] = deepcopy(y.coeffs[j])
    monomial_set!(x, k, y, j)
    j += 1
    k += 1
  end

  x.len = k - 1
  resize!(x.coeffs, x.len)
  return x
end

function mul!(x::MPolyRingElem{T}, a::T) where {T}
  if iszero(a)
    return zero!(x)
  elseif isone(a)
    return x
  end

  for i in 1:length(x)
    x.coeffs[i] = mul!(x.coeffs[i], a)
  end
  return x
end
