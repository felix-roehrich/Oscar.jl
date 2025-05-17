struct LusztigDatum <: AbstractCrystalElem
  # datum wrt. to _w0
  datum::Vector{Int}
  _w0::Vector{UInt8}
  root_system::RootSystem

  # w0 used for display purposes
  w0::Vector{UInt8}
end

function Base.deepcopy_internal(d::LusztigDatum, id::IdDict)
  return get!(
    id,
    d,
    LusztigDatum(
      deepcopy_internal(d.datum, id),
      deepcopy_internal(d._w0, id),
      d.root_system,
      deepcopy_internal(d.w0, id),
    ),
  )
end

function Base.hash(d::LusztigDatum, h::UInt)
  b = 0xbf09cbc877180763 % UInt
  h = hash(d.datum, h)
  h = hash(root_system(d), h)
  h = hash(d.w0, h)

  return xor(h, b)
end

function lusztig_datum(R::RootSystem, w0::Vector{UInt8})
  return lusztig_datum(R, zeros(Int, number_of_positive_roots(R)), w0)
end

function lusztig_datum(R::RootSystem, datum::Vector{Int}, w0::Vector{UInt8})
  @req number_of_positive_roots(R) == length(datum) == length(w0) "length does not match number of positive roots"
  @req weyl_group(R)(w0) == longest_element(weyl_group(R)) "w0 not a reduced decomposition"
  return LusztigDatum(datum, copy(w0), R, w0)
end

function lusztig_datum(d::LusztigDatum, w0::Vector{<:Integer})
  datum = copy(d.datum)
  w = UInt8.(w0)
  for mv in braid_moves(weyl_group(d), w, d._w0)
    _move!(datum, mv)
  end
  return LusztigDatum(datum, copy(w), d.root_system, w)
end

function Base.show(io::IO, d::LusztigDatum)
  for mv in braid_moves(weyl_group(d), d.w0, d._w0)
    _move!(d.datum, mv)
  end
  copyto!(d._w0, d.w0)

  print(io, "Lusztig datum (")
  join(io, d.datum, ", ")
  print(io, ") for w0 = (")
  join(io, Iterators.map(Int, d.w0), ", ")
  print(io, ")")
end

function root_system(datum::LusztigDatum)
  return datum.root_system
end

function weyl_group(datum::LusztigDatum)
  return weyl_group(datum.root_system)
end

function adapted_string(d::LusztigDatum)
  s = zeros(Int, length(d.w0))

  wn = copy(d._w0)
  wo = copy(d._w0)
  datum = copy(d.datum)
  for i in 1:length(s)
    exchange_left!(weyl_group(d), wn, d.w0[i])
    for mv in braid_moves(weyl_group(d), wn, wo)
      _move!(datum, mv)
    end
    copyto!(wo, wn)

    s[i] = datum[1]
    datum[1] = 0
  end

  return s
end

function Crystal.e!(d::LusztigDatum, i::Int)
  _state!(d, UInt8(i))

  d.datum[1] -= 1
  if d.datum[1] < 0
    return nothing
  end
  return d
end

function Crystal.f!(d::LusztigDatum, i::Int)
  _state!(d, UInt8(i))
  d.datum[1] += 1
  return d
end

function _state!(d::LusztigDatum, i::UInt8)
  if d._w0[1] == i
    return nothing
  end

  w0 = copy(d._w0)
  exchange_left!(weyl_group(d), w0, i)
  for mv in braid_moves(weyl_group(d), w0, d._w0)
    _move!(d.datum, mv)
  end
  copyto!(d._w0, w0)
end

function _move!(d::Vector{Int}, mv::Tuple{Int,Int,Int})
  n, len, dir = mv
  if len == 2
    # A1xA1 move
    d[n], d[n + 1] = d[n + 1], d[n]
  elseif len == 3
    # A2 move
    m = min(d[n], d[n + 2])
    d[n], d[n + 1], d[n + 2] = d[n + 1] + d[n + 2] - m, m, d[n] + d[n + 1] - m
  elseif len == 4
    # B2/C2 move
    if dir == -1
      reverse!(d, n, n + 3)
    end
    m1 = min(d[n], d[n + 2])
    m2 = min(d[n] + d[n + 1], m1 + d[n + 3])
    m3 = min(2 * d[n] + d[n + 1], 2 * m1 + d[n + 3])

    d[n], d[n + 1], d[n + 2], d[n + 3] = d[n + 1] + 2 * d[n + 2] + d[n + 3] - m3,
    m3 - m2, 2 * m2 - m3,
    d[n] + d[n + 1] + d[n + 2] - m2

    if dir == -1
      reverse!(d, n, n + 3)
    end
  elseif len == 6
    if dir == -1
      reverse!(d, n, n + 5)
    end
    m1 = min(d[n], d[n + 2])
    m2 = min(d[n + 2], d[n + 4])
    m3 = min(d[n] + d[n + 2], 2 * d[n + 2], d[n + 2] + d[n + 4], d[n] + d[n + 4])
    m4 = min(
      d[n] + d[n + 1] + 2 * d[n + 2] + d[n + 3],
      d[n] + d[n + 1] + 2 * m2 + d[n + 5],
      m1 + d[n + 3] + 2 * d[n + 4] + d[n + 5],
    )
    m5 = min(
      2 * d[n] + 2 * d[n + 1] + 3 * d[n + 2] + d[n + 3],
      2 * d[n] + 2 * d[n + 1] + 3 * m2 + d[n + 5],
      2 * m1 + 2 * d[n + 3] + 3 * d[n + 4] + d[n + 5],
      d[n] + d[n + 1] + d[n + 3] + 2 * d[n + 4] + d[n + 5] + m3,
    )
    m6 = min(
      3 * d[n] + 2 * d[n + 1] + 3 * d[n + 2] + d[n + 3],
      3 * d[n] + 2 * d[n + 1] + 3 * m2 + d[n + 5],
      3 * m1 + 2 * d[n + 3] + 3 * d[n + 4] + d[n + 5],
      2 * d[n] + d[n + 1] + d[n + 3] + 2 * d[n + 4] + d[n + 5] + m3,
    )
    m7 = min(
      2 * d[n] + 2 * d[n + 1] + 3 * d[n + 2] + d[n + 3] +
      min(
        d[n] + d[n + 1] + 3 * d[n + 2] + d[n + 3],
        d[n] + d[n + 1] + 3 * m2 + d[n + 5],
        m3 + d[n + 3] + 2 * d[n + 4] + d[n + 5],
      ),
      2 * d[n + 5] + 3 * min(d[n] + d[n + 1] + 2 * m2, m1 + d[n + 3] + 2 * d[n + 4]),
    )

    d[n], d[n + 1], d[n + 2], d[n + 3], d[n + 4], d[n + 5] = d[n + 1] + 3 * d[n + 2] +
                                                             2 * d[n + 3] + 3 * d[n + 4] +
                                                             d[n + 5] - m6,
    m6 - m5, 3 * m5 - m6 - m7, m7 - m4 - m5, 3 * m4 - m7,
    d[n] + d[n + 1] + 2 * d[n + 2] + d[n + 3] + d[n + 4] - m4

    if dir == -1
      reverse!(d, n, n + 5)
    end
  end
end
