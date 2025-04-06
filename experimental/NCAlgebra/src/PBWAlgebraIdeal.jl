
function left_ideal(x::Vector{PBWAlgebraElem{T}}) where {T}
  return PBWAlgebraIdeal(x, PBWAlgebraElem{T}[])
end

function Base.in(x::PBWAlgebraElem{T}, I::PBWAlgebraIdeal{T}) where {T}
  return iszero(rem(x, groebner_basis(I)))
end

function rem!(z::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}, I::PBWAlgebraIdeal{T}) where {T}
  z = rem!(z, x, groebner_basis(I))
  return z
end