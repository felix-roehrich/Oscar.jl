```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Crystals

According to Kashiwara a crystal is a set $B$ together with the following maps
```
\begin{aligned}
  \tilde{f}_i, \tilde{e}_i: B \to B \cup \{0\} \\
  \varphi_i, \epsilon: B\to \mathbb{Z} \\
  weight: B \to P
\end{aligned}
```

This is realised in OSCAR through the `AbstractCrystal` and `AbstractCrystalElem` interfaces and the following methods

```@docs
Crystal.e(::AbstractCrystalElem, i::Int)
Crystal.f(::AbstractCrystalElem, i::Int)
Crystal.eps(::AbstractCrystalElem, i::Int)
Crystal.phi(::AbstractCrystalElem, i::Int)
Crystal.weight(::AbstractCrystalElem)
```

Additionally, the following methods are provided to act on a crystal element
```@docs
Crystal.e!(::AbstractCrystalElem, i::Int)
Crystal.e(::AbstractCrystalElem, i::Int, n::Int)
Crystal.e!(::AbstractCrystalElem, i::Int, n::Int)
Crystal.e(::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
Crystal.e!(::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
```

```@docs
Crystal.f!(::AbstractCrystalElem, i::Int)
Crystal.f(::AbstractCrystalElem, i::Int, n::Int)
Crystal.f!(::AbstractCrystalElem, i::Int, n::Int)
Crystal.f(::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
Crystal.f!(::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
```
