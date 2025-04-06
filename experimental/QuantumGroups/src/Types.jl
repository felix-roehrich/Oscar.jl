struct QuantumGroup{T<:FieldElem,S} <: NCRing # S needed for Singular data
  A::Any#::LaurentPolynomialRing

  alg::PBWAlgRing{T,S}
  q::T #::RationalFunctionFieldElem{QQFieldElem,QQPolyRingElem}
  qi::Vector{T} # {RationalFunctionFieldElem{QQFieldElem,QQPolyRingElem}}
  root_system::RootSystem
  w0::Vector{UInt8}

  cvx::Vector{Int}

  # cache
  canonical_basis::Dict{Vector{Int},PBWAlgElem}
end

mutable struct QuantumGroupElem{T<:FieldElem,S} <: NCRingElem
  U::QuantumGroup{T,S}
  elem::PBWAlgElem{T,S}
end

struct QuantumGroupHomomorphism
  img::Vector{QuantumGroupElem}
end
