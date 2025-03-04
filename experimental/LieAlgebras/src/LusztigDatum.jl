struct LusztigDatum
  datum::Vector{UInt8}
  w0::Vector{UInt8}
  root_system::RootSystem

  # internal state for computations
  _datum::Vector{UInt8}
  _w0::Vector{UInt8}

  function LusztigDatum(datum::Vector{UInt8}, w0::Vector{UInt8})
    new(datum, w0, copy(datum), copy(w0))
  end
end

function root_system(datum::LusztigDatum)
  return datum.root_system
end

function weyl_group(datum::LusztigDatum)
  return weyl_group(datum.root_system)
end

function Crystal.e!(datum::LusztigDatum, i::Int)
  ui = UInt8(i)
  W = weyl_group(datum)
  if datum._w0[1] != ui
    w0 = copy(datum._w0)
    exchange!(W, w0, ui)
    braid_moves(W, w0, datum._w0)
  end

  datum._datum[1] -= 1
  if datum._datum[1] < 0
    return nothing
  end
  return datum
end
