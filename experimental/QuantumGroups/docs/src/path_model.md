```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Path model

!!! note
    Since groups in OSCAR act from the right, the action of the root operators is likewise from the right.

```@docs
PathModel.f(::AbstractPathModelElem, ::Int) -> AbstractPathModelElem
PathModel.f!(::AbstractPathModelElem, ::Int) -> AbstractPathModelElem

PathModel.f(::AbstractPathModelElem, ::Int, ::Int) -> AbstractPathModelElem
PathModel.f!(::AbstractPathModelElem, ::Int, ::Int) -> AbstractPathModelElem

PathModel.f(::AbstractPathModelElem, ::Vector{Int}, ::Vector{Int}) -> AbstractPathModelElem
PathModel.f!(::AbstractPathModelElem, ::Vector{Int}, ::Vector{Int}) -> AbstractPathModelElem

PathModel.e(::AbstractPathModelElem, ::Int) -> AbstractPathModelElem
PathModel.e!(::AbstractPathModelElem, ::Int) -> AbstractPathModelElem

PathModel.e(::AbstractPathModelElem, ::Int, ::Int) -> AbstractPathModelElem
PathModel.e!(::AbstractPathModelElem, ::Int, ::Int) -> AbstractPathModelElem

PathModel.e(::AbstractPathModelElem, ::Vector{Int}, ::Vector{Int}) -> AbstractPathModelElem
PathModel.e!(::AbstractPathModelElem, ::Vector{Int}, ::Vector{Int}) -> AbstractPathModelElem

PathModel.h(::AbstractPathModelElem, ::Int) -> Int
```

```@docs
PathModel.eps(::AbstractPathModelElem, ::Int) -> Int
PathModel.phi(::AbstractPathModelElem, ::Int) -> Int
PathModel.weight(::AbstractPathModelElem) -> WeightLatticeElem
```

```@docs
dominant_path(::AbstractPathModel) -> AbstractPathModelElem
highest_weight(::AbstractPathModel) -> WeightLatticeElem
```

```@docs
is_dominant_path(::PathModelElem) -> Bool
```

# LS path model

```@docs
ls_path_model(::RootSystem, ::WeightLatticeElem)
```
