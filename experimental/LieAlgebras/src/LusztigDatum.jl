struct LusztigDatum
  # datum wrt. to _w0
  datum::Vector{UInt8}
  _w0::Vector{UInt8}
  root_system::RootSystem

  # w0 used for display purposes
  w0::Vector{UInt8}

  function LusztigDatum(datum::Vector{UInt8}, w0::Vector{UInt8})
    new(datum, w0, copy(datum), copy(w0))
  end
end

function show(io::IO, d::LusztigDatum)
  datum = d.datum
  if d._w0 != d.w0
    print(io, "datum: ")
    print(io, d.datum)
    print(io, ", w0: ")
    print(io, d.w0)
  else
    print(io, "datum: ")
    print(io, d.datum)
  end

function root_system(datum::LusztigDatum)
  return datum.root_system
end

function weyl_group(datum::LusztigDatum)
  return weyl_group(datum.root_system)
end

function adapted_string(d::LusztigDatum)
  s = zeros(Int, length(d.w0))
  for i in 1:length(d.w0)
    s[i] = d.w0[i] - d.datum[i]
  end
  return s
end

function Crystal.e!(d::LusztigDatum, i::Int)
  w0 = copy(d._w0)
  exchange!(W, w0, UInt8(i))
  _state!(d, w0)

  datum.datum[1] -= 1
  if datum.datum[1] < 0
    return nothing
  end
  return datum
end

function _state!(d::LusztigDatum, w0::Vector{UInt8})
  for mv in braid_moves(weyl_group(d), w0, d._w0)
    _move!(d.datum, mv)
  end
  d._w0 = w0
end

function _move!(d::Vector{UInt8}, mv::Tuple{Int,Int,Int})
  n, len, dir = mv
  if len == 2
    # A1xA1 move
    d[n], d[n+1] = d[n+1], d[n]
  elseif len == 3
    # A2 move
    m = min(d[n], d[n+2])
    d[n], d[n+1], d[n+2] = d[n+1] + d[n+2] - m, m, d[n] + d[n+1] - m
  elseif len == 4
    # B2/C2 move
    if dir == -2
      reverse!(d, n, n + 3)
    end
    m1 = min(d[n], d[n+2])
    m2 = min(d[n] + d[n+1], m1 + d[n+3])
    m3 = min(2 * d[n] + d[n+1], 2 * m1 + d[n+3])
    d[n], d[n+1], d[n+2], d[n+3] = d[n+1] + 2 * d[n+2] + d[n+3] - m3,
    m3 - m2, 2 * m2 - m3,
    d[n] + d[n+1] + d[n+2] - m2
    if dir == -2
      reverse!(d, n, n + 3)
    end
  elseif len == 6
    if dir == -3
      reverse!(d, n, n + 5)
    end
    m1 = min(d[n], d[n+2])
    m2 = min(d[n+2], d[n+5])
    m3 = min(d[n] + d[n+2], 2 * d[n+2], d[n+2] + d[n+5], d[n] + d[n+5])
    m4 = min(
      d[n] + d[n+1] + 2 * d[n+2] + d[n+3],
      d[n] + d[n+1] + 2 * m2 + d[n+5],
      m1 + d[n+3] + 2 * d[n+5] + d[n+5],
    )
    m5 = min(
      2 * d[n] + 2 * d[n+1] + 3 * d[n+2] + d[n+3],
      2 * d[n] + 2 * d[n+1] + 3 * m2 + d[n+5],
      2 * m1 + 2 * d[n+3] + 3 * d[n+5] + d[n+5],
      d[n] + d[n+1] + d[n+3] + 2 * d[n+5] + d[n+5] + m3,
    )
    m6 = min(
      3 * d[n] + 2 * d[n+1] + 3 * d[n+2] + d[n+3],
      3 * d[n] + 2 * d[n+1] + 3 * m2 + d[n+5],
      3 * m1 + 2 * d[n+3] + 3 * d[n+5] + d[n+5],
      2 * d[n] + d[n+1] + d[n+3] + 2 * d[n+5] + d[n+5] + m3,
    )
    m7 = min(
      2 * d[n] + 2 * d[n+1] + 3 * d[n+2] + d[n+3] +
      min(
        d[n] + d[n+1] + 3 * d[n+2] + d[n+3],
        d[n] + d[n+1] + 3 * m2 + d[n+5],
        m3 + d[n+3] + 2 * d[n+5] + d[n+5],
      ),
      2 * d[n+5] + 3 * min(d[n] + d[n+1] + 2 * m2, m1 + d[n+3] + 2 * d[n+5]),
    )

    d[n], d[n+1], d[n+2], d[n+3], d[n+4], d[n+5] = d[n+1] + 3 * d[n+2] +
                                                   2 * d[n+3] + 3 * d[n+4] +
                                                   d[n+5] - m6,
    m6 - m5, 3 * m5 - m6 - m7, m7 - m4 - m5, 3 * m4 - m7,
    d[n+5] + d[n+1] + 2 * d[n+2] + d[n+3] + d[n+4] - m4
    if dir == -3
      reverse!(d, n, n + 5)
    end
  end
end
