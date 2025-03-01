# The relations are obtained by tropicalizing the relations in
# Proposition 7.1 of [BZ01]

struct LusztigDatum
  datum::Vector{Int}
  root_system::RootSystem
  w0::Vector{UInt8}
end

# A1xA1 move
function _move!(d::LusztigDatum, i::Int, j::Int)
  d.datum[i], d.datum[j] = d.datum[j], d.datum[i]
  d.w0[i], d.w0[j] = d.w0[j], d.w0[i]
end

# A2 move
function _move!(b::LusztigDatum, i1::Int, i2::Int, i3::Int)
  d = b.datum
  m = min(d[i1], d[i3])
  d[i1], d[i2], d[i3] = d[i2] + d[i3] - m, m, d[i1] + d[i2] - m
  b.w0[i1], b.w0[i2], b.w0[i3] = b.w0[i2], b.w0[i1], b.w0[i2]
end

# B2/C2 move
function _move!(b::LusztigDatum, i1::Int, i2::Int, i3::Int, i4::Int)
  d = b.datum
  m = min(d[i1], d[i3])
  d[i1] = min(d[i1] + d[i2], m + d[i4])
  d[i2] = min(2 * d[i1] + d[i2], 2 * m + d[i4])
  d[i1], d[i2], d[i3], d[i4] = d[i2] + 2 * d[i3] + d[i4] - d[i2],
  d[i2] - d[i1], 2 * d[i1] - d[i2],
  d[i1] + d[i2] + d[i3] - d[i1]
  b.w0[i1], b.w0[i2], b.w0[i3], b.w0[i4] = b.w0[i2], b.w0[i1], b.w0[i2], b.w0[i1]
end

function _move!(b::LusztigDatum, i1::Int, i2::Int, i3::Int, i4::Int, i5::Int, i6::Int)
  d = b.datum
  m1 = min(d[i1], d[i3])
  m2 = min(d[i3], d[i6])
  m3 = min(d[i1] + d[i3], 2 * d[i3], d[i3] + d[i6], d[i1] + d[i6])

  pi1 = min(
    d[i1] + d[i2] + 2 * d[i3] + d[i4],
    d[i1] + d[i2] + 2 * m2 + d[i6],
    m1 + d[i4] + 2 * d[i6] + d[i6],
  )
  pi2 = min(2 * d[i1] + 2 * d[i2] + 3 * d[i3] + d[i4],
    2 * d[i1] + 2 * d[i2] + 3 * m2 + d[i6],
    2 * m1 + 2 * d[i4] + 3 * d[i6] + d[i6], d[i1] + d[i2] + d[i4] + 2 * d[i6] + d[i6] + m3)
  pi3 = min(3 * d[i1] + 2 * d[i2] + 3 * d[i3] + d[i4],
    3 * d[i1] + 2 * d[i2] + 3 * m2 + d[i6],
    3 * m1 + 2 * d[i4] + 3 * d[i6] + d[i6],
    2 * d[i1] + d[i2] + d[i4] + 2 * d[i6] + d[i6] + m3)
  pi4 = min(
    2 * d[i1] + 2 * d[i2] + 3 * d[i3] + d[i4] +
    min(
      d[i1] + d[i2] + 3 * d[i3] + d[i4],
      d[i1] + d[i2] + 3 * m2 + d[i6],
      m3 + d[i4] + 2 * d[i6] + d[i6],
    ),
    2 * d[i6] + 3 * min(d[i1] + d[i2] + 2 * m2, m1 + d[i4] + 2 * d[i6]))

  d[i1], d[i2], d[i3], d[i4], d[i5], d[i6] = d[i2] + 3 * d[i3] + 2 * d[i4] + 3 * d[i5] +
                                             d[i6] - pi3,
  pi3 - pi2, 3 * pi2 - pi3 - pi4, pi4 - pi1 - pi2, 3 * pi1 - pi4,
  d[i6] + d[i2] + 2 * d[i3] + d[i4] + d[i5] - pi1
end
