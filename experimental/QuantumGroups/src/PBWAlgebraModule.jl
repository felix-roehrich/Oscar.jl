function zero(M::PBWAlgebraModule)
  return PBWAlgebraModuleElem(zero(M.I.gens[1]), M)
end

function Base.show(io::IO, x::PBWAlgebraModuleElem)
  print(io, x.elem)
end

function parent(m::PBWAlgebraModuleElem)
  return m.parent
end

function zero(m::PBWAlgebraModuleElem)
  return zero(parent(m))
end

function zero!(m::PBWAlgebraModuleElem)
  zero!(m.elem)
  return m
end

function add!(
  z::PBWAlgebraModuleElem{T}, x::PBWAlgebraModuleElem{T}, y::PBWAlgebraModuleElem{T}
) where {T}
  add!(z.elem, x.elem, y.elem)
  return z
end

function mul!(
  z::PBWAlgebraModuleElem{T}, x::PBWAlgebraElem{T}, m::PBWAlgebraModuleElem{T}
) where {T}
  rem!(mul!(z.elem, x, m.elem), parent(z).I)
  return z
end

function Base.:+(x::PBWAlgebraModuleElem, y::PBWAlgebraModuleElem)
  return add!(zero(x), x, y)
end

function Base.:*(x::PBWAlgebraElem{T}, m::PBWAlgebraModuleElem{T}) where {T}
  return mul!(zero(m), x, m)
end

function Base.:(==)(x::PBWAlgebraModuleElem, y::PBWAlgebraModuleElem)
  return x.elem == y.elem
end
