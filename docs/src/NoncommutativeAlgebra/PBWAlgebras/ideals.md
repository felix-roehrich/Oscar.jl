```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Ideals in PBW-algebras

## Types

The OSCAR type for ideals in PBW-algebras is of parametrized form
`PBWAlgIdeal{D, T, S}`, where `D` encodes the direction left, right,
or two-sided, and `T` is the element type of the field over which
the PBW-algebra is defined (the type `S` is added for internal use).

## Constructors

```@docs
left_ideal(g::Vector{<:PBWAlgElem})
right_ideal(g::Vector{<:PBWAlgElem})
two_sided_ideal(g::Vector{<:PBWAlgElem})
```

## Gröbner bases

## Data Associated to Ideals

If `I` is an ideal of a PBW-algebra  `A`, then

- `base_ring(I)` refers to `A`,
- `gens(I)` to the generators of `I`,
- `number_of_generators(I)` / `ngens(I)` to the number of these generators, and
- `gen(I, k)` as well as `I[k]` to the `k`-th such generator.

###### Examples

```jldoctest
julia> D, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"])
(Weyl-algebra over rational field in variables (x, y), PBWAlgElem{QQFieldElem, Singular.n_Q}[x, y, dx, dy])

julia> I = left_ideal(D, [x, dx])
left_ideal(x, dx)

julia> base_ring(I)
Weyl-algebra over rational field in variables (x, y)

julia> gens(I)
2-element Vector{PBWAlgElem{QQFieldElem, Singular.n_Q}}:
 x
 dx

julia> number_of_generators(I)
2

julia> gen(I, 2)
dx

```

## Operations on Ideals

### Simple Ideal Operations

#### Powers of Ideal

```@docs
^(I::PBWAlgIdeal{D, T, S}, k::Int) where {D, T, S}
```

#### Sum of Ideals

```@docs
+(I::PBWAlgIdeal{D, T, S}, J::PBWAlgIdeal{D, T, S}) where {D, T, S}
```

#### Product of Ideals

```@docs
*(I::PBWAlgIdeal{DI, T, S}, J::PBWAlgIdeal{DJ, T, S}) where {DI, DJ, T, S}
```

### Intersection of Ideals

```@docs
intersect(I::PBWAlgIdeal{D, T, S}, Js::PBWAlgIdeal{D, T, S}...) where {D, T, S}
intersect(V::Vector{PBWAlgIdeal{D, T, S}}) where {D, T, S}
```

### Elimination

Let

```math
A = K\langle x_1, \dots , x_n \mid x_jx_i = c_{ij}x_ix_j+d_{ij},  \ 1\leq i<j \leq n \rangle,
```
be a PBW-algebra. Fix a subset $\sigma\subset \{1,\dots, n\}$, write $x_\sigma$ 
for the set of variables $x_i$ with $i\in\sigma$, and let $A_\sigma$ be the $K$-linear 
subspace of $A$ which is generated by the standard monomials in $\langle x_\sigma \rangle$.
Suppose there exists a global monomial ordering $>$ on $A$ which is both admissible for $A$
and an elimination ordering for $x\smallsetminus x_\sigma$. Then $A_\sigma$ is a subalgebra
of $A$ with $d_{ij}\in A_\sigma$ for each pair of indices $1\leq i<j \leq n$ with $i,j\in\sigma$.
In particular, $A_\sigma$  is a PBW-algebra with admissible ordering $>_\sigma$, where $>_\sigma$
is the restriction of $>$ to the set of standard monomials in  $\langle x_\sigma\rangle$. Moreover,
if $I\subset A$ is a nonzero (left, right, two-sided) ideal, and $\mathcal G$ is a (left, right, two-sided)
Gröbner basis for $I$ with respect to $>$, then $\mathcal G\cap A_\sigma$ is a (left, right, two-sided)
Gröbner basis for $I\cap A_\sigma$ with respect to $>_\sigma$. We refer to computing $I\cap A_\sigma$
as *eliminating the variables in $x\smallsetminus x_\sigma$ from $I.$*

!!! note
    If the relevant $d_{ij}$ are all contained in $A_\sigma$, we also say that $A_\sigma$ is *admissible for elimination*.

!!! note
    A monomial ordering which is both admissible for $A$ and an elimination ordering for $x\smallsetminus x_\sigma$ may not exist. 

```@docs
eliminate(I::PBWAlgIdeal, V::Vector{<:PBWAlgElem}; ordering = nothing)
```

## Tests on Ideals

```@docs
is_zero(I:: PBWAlgIdeal)
```

```@docs
is_one(I:: PBWAlgIdeal)
```

```@docs
is_subset(I::PBWAlgIdeal{D, T, S}, J::PBWAlgIdeal{D, T, S}) where {D, T, S}
```

```@docs
==(I::PBWAlgIdeal{D, T, S}, J::PBWAlgIdeal{D, T, S}) where {D, T, S}
```

```@docs
ideal_membership(f::PBWAlgElem{T, S}, I::PBWAlgIdeal{D, T, S}) where {D, T, S}
```

