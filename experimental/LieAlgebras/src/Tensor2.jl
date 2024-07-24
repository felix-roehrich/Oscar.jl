import Nemo: pow

# p partition, b bound
function next_partition(p::Vector{Int}, b::Vector{Int})
  np = copy(p)
  len = length(np)

  r = np[len]+1
  np[len] = 0
  for i in len-1:-1:1
    if np[i] != 0
      # overflow over the bound
      of = np[i+1]+r - b[i+1]
      if of <= 0
        np[i] -= 1
        np[i+1] += r
        return np
      elseif of > sum(b[j] for j in i+2:len; init=0)
        r += np[i]
        np[i] = 0
        continue
      end

      np[i] -= 1
      np[i+1] += r
      for j in i+2:len
        np[j-1] -= of
        np[j] += of
        of -= b[j]
        if of <= 0
          return np
        end
      end
    end
  end

  return nothing
end

A, v = cyclotomic_field(12)
RR = root_system(:B, 2)
sr = gens(weyl_group(RR))

B = [4 -2; -2 2]
C = [2 -1; -2 2]

lam = [2,2]
f = [2,1,2,1]
p = [12, 24, 36, 12]
pbw = [4, 1, 4, 1]

GAP.Globals.LoadPackage(GAP.Obj("QuaGroup"))
R = GAP.Globals.RootSystem(GAP.Obj("B"), 2)
U = GAP.Globals.QuantizedUEA(R)
g = GAP.Globals.GeneratorsOfAlgebra(U)
V = GAP.Globals.HighestWeightModule(U, GAP.Obj([2,2]))
v0 = GAP.Globals.Basis(V)[1]

function pow(f::GAP.GapObj, n::Int)
  if n == 0
    return GAP.Globals.One(U)
  end

  fam = GAP.Globals.ElementsFamily(GAP.Globals.FamilyObj(U))
  rep = GAP.Globals.List(GAP.Globals.ExtRepOfObj(f), GAP.Globals.ShallowCopy)
  rep[1][2] = n
  return GAP.Globals.ObjByExtRep(fam, rep)
end

pbw_212_12_2 = [v^-4, v^-4, v^-4, v^-4, v^-20, v^-20, v^-20, v^-20, v^-20, A(1), A(1), A(1)]
m_212_12_2 = [
2 2 2 2 2 2 0 0 0 0 0 0
4 4 4 4 4 4 0 0 0 0 0 0
6 6 6 6 0 0 2 2 2 2 2 2
0 0 0 0 0 0 4 4 4 0 0 0
]
n_212_12_2 = [
2 2 2 2 2 2 0 0 0 0 0 0
6 6 6 6 2 2 2 2 2 2 0 0
2 2 2 2 2 2 4 4 4 0 0 0
2 2 2 2 0 0 0 0 0 0 2 2
]

m_212_1 = [
2 2 2 2 2 2 0 0 0 0 0 0
4 4 4 4 4 4 0 0 0 0 0 0
6 6 6 6 6 6 0 0 0 0 0 0
0 0 0 0 0 0 2 2 2 2 2 2
]

pbw_121_21_1 = [v^-8, v^-8, v^-8, A(1), A(1), A(1), A(1), A(1), A(1), A(1), A(1), A(1)]
m_121_21_1 = [
2 2 2 2 2 2 0 0 0 0 0 0
4 4 4 2 2 2 2 2 2 0 0 0
4 4 4 4 4 4 6 6 0 0 0 0
2 2 2 0 0 0 0 0 0 2 2 2
]
[
0 0 2 2 2 2 2 2 0 0 0 0
2 2 2 2 2 2 2 2 2 2 2 2
6 6 4 4 4 4 4 4 0 0 0 0
4 4 4 0 0 0 0 0 0 0 0 0
]
[
0 0 2 2 2 2 2 2 0 0 0 0
2 2 3 2 2 2 2 2 1 2 2 2
6 6 4 4 4 4 4 4 0 0 0 0
4 4 3 0 0 0 0 0 1 0 0 0
]
n_121_21_1 = [
2 2 2 2 2 2 0 0 0 0 0 0
6 6 6 6 6 6 0 0 0 0 0 0
4 4 4 0 0 0 2 2 2 2 2 2
0 0 0 0 0 0 6 6 0 0 0 0
]

m_121_2 = [
2 2 2 2 2 2 0 0 0 0 0 0
4 4 4 4 4 4 0 0 0 0 0 0
4 4 4 4 4 4 2 2 2 2 2 2
2 2 2 2 2 2 0 0 0 0 0 0
]

struct TensorIterator
  col::Vector{Vector{Int}}
  highest_weight::Vector{Int}
  pbw::Vector{nf_elem}
  row::Vector{Int}
end

function Base.iterate(iter::TensorIterator)
  sh = collect(1:iter.p+iter.q)
  return sh, (sh, iter.p+1)
end

function Base.iterate(iter::TensorIterator, t::Matrix{Int})
  n, m = size(t)
  nt = copy(t)

  @label loop
  nr = next_partition(nt[2, :], compute_bound(iter.highest_weight, nt, 2, iter.col))
  if nr != nothing
    nt[2, :] = nr
  else
    nr = next_partition(nt[1, :], compute_bound(iter.highest_weight, nt, 1, iter.col))
    if nr == nothing
      return nothing
    end

    nt[1, :] = nr

    b = compute_bound(lam, nt, 2, col)
    r = iter.row[2]
    for j in 1:m
      _min = min(r, b[j])
      nt[2, j] = _min
      r -= _min
    end
  end

  b3 = compute_bound(lam, nt, 3, col)
  for j in 1:m
    nt[3, j] = col[j][f[3]] - nt[1, j]
    nt[4, j] = col[j][f[4]] - nt[2, j]

    if nt[3, j] > b3[j]
      @goto loop
    end
  end

  return nt, nt
end

function tensor_iterator(t::Matrix{Int})
  n, m = size(t)

  col = [[sum(t[i, j] for i in 1:n if f[i] == k) for k in 1:length(lam)] for j in 1:m]
  row = [sum(t[i,j] for j in 1:m) for i in 1:n]
  pbw_coeffs = zeros(A, m)
  for j in 1:m
    wv = v0
    for i in 1:n
      wv = pow(g[pbw[i]], t[i,j])^wv
    end
    
    rep = GAP.Globals.ExtRepOfObj(GAP.Globals.ExtRepOfObj(wv))
    if iszero(rep)
      pbw_coeffs[j] = 0
      continue
    end
    
    coeffs = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[2])
    pbw_coeffs[j] = sum(l -> coeffs[1][l]*v^(l-1+coeffs[2]), 1:length(coeffs[1]))
  end
  
  return TensorIterator(col, pbw_coeffs, row)
end

function compute_bound(hw::Vector{Int}, t::Matrix{Int}, s::Int, col::Vector{Vector{Int}})
  n, m = size(t)

  b = zeros(Int, m)
  if s == 1
    for j in 1:m
      b[j] = min(hw[f[s]], col[j][f[s]])
    end
    return b
  elseif s == 2
    for j in 1:m
      b[j] = min(hw[f[s]]-t[1, j]*C[f[s], f[1]], col[j][f[s]])
    end
    return b
  elseif s == 3
    for j in 1:m
      b[j] = min(hw[f[s]]-t[1, j]-t[2, j]*C[f[s], f[2]], col[j][f[s]])
    end
    return b
  end
  
  error("not implemented")
end

function next_tensor(mat::Matrix{Int})
  n, m = size(mat)
  col = [[sum(mat[i, j] for i in 1:n if f[i] == k) for k in 1:length(lam)] for j in 1:m]
  row = [sum(mat[i,j] for j in 1:m) for i in 1:n]
  nt = copy(mat)

  @label loop
  nr = next_partition(nt[2, :], compute_bound(lam, nt, 2, col))
  if nr != nothing
    nt[2, :] = nr
  else
    nr = next_partition(nt[1, :], compute_bound(lam, nt, 1, col))
    if nr == nothing
      return nothing
    end

    nt[1, :] = nr

    b = compute_bound(lam, nt, 2, col)
    r = row[2]
    for j in 1:m
      _min = min(r, b[j])
      nt[2, j] = _min
      r -= _min
    end
    if r > 0
      @goto loop
    end
  end

  b3 = compute_bound(lam, nt, 3, col)
  for j in 1:m
    nt[3, j] = col[j][f[3]] - nt[1, j]
    nt[4, j] = col[j][f[4]] - nt[2, j]

    if nt[3, j] > b3[j]
      @goto loop
    end
  end
  
  return nt
end

# lam = coeff vector of simple roots
function tensor_coefficient(lam, mat)
  q = 1
  w = [copy(lam) for _ in 1:12]

  for i in 1:nrows(mat)
    s = p[i]
    for j in ncols(mat):-1:2
      s -= mat[i,j]
      w[j][f[i]] -= mat[i,j]
      q *= v^-(s*mat[i,j]+s*dot(B[f[i], :], w[j]))
    end
    w[1][f[i]] -= mat[i,1]
  end

  for i in 1:ncols(mat)
    obj = foldl((v, l) -> pow(g[pbw[l]], mat[l, i])^v, 1:nrows(mat); init=v0)
    rep = GAP.Globals.ExtRepOfObj(GAP.Globals.ExtRepOfObj(obj))

    if iszero(rep)
      q = 0
      break
    end
    coeffs = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[2])
    q *= sum(l -> coeffs[1][l]*v^(l-1+coeffs[2]), 1:length(coeffs[1])) / pbw_coeffs[i]
  end
  
  return q
end

function coefficient(iter::TensorIterator, t::Matrix{Int})
  q = 1
  
  w = [C*iter.highest_weight for _ in 1:12]
  for i in 1:nrows(t)
    s = iter.rows[i]
    for j in ncols(t):-1:2
      s -= t[i,j]
      w[j][f[i]] -= t[i,j]
      q *= v^-(s*t[i,j]+s*dot(B[f[i], :], w[j]))
    end
    w[1][f[i]] -= t[i,1]
  end

  for i in 1:ncols(t)
    obj = foldl((v, l) -> pow(g[pbw[l]], t[l, i])^v, 1:nrows(t); init=v0)
    rep = GAP.Globals.ExtRepOfObj(GAP.Globals.ExtRepOfObj(obj))

    if iszero(rep)
      q = 0
      break
    end
    coeffs = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[2])
    q *= sum(l -> coeffs[1][l]*v^(l-1+coeffs[2]), 1:length(coeffs[1])) / pbw_coeffs[i]
  end

  return q
end
