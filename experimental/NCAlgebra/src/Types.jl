
struct PBWAlgebraModule{T} <: AbstractAlgebra.Module{T}
  I::PBWAlgebraIdeal{T}
end

struct PBWAlgebraModuleElem{T} <: AbstractAlgebra.ModuleElem{T}
  elem::PBWAlgebraElem{T}
  parent::PBWAlgebraModule{T}
end
