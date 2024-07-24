function shuffles(p, q)
  return ShuffleIterator(p, q)
end

struct qparam
end

struct ShuffleIterator
  p::Int
  q::Int
end

function Base.eltype(::Type{ShuffleIterator})
  return Vector{Int}
end

function Base.iterate(iter::ShuffleIterator)
  sh = collect(1:iter.p+iter.q)
  return sh, (sh, iter.p+1)
end

function Base.iterate(iter::ShuffleIterator, last::Tuple{Vector{Int}, Int})
  sh, n = last
  for i in max(1, n-1):iter.p
    if i == iter.p || sh[i]+1 < sh[i+1]
      if sh[i] == length(sh)
        return nothing
      end

      sh[1:i-1] = 1:i-1
      sh[i] += 1
      sh[iter.p+1:iter.p+sh[i]-i] = i:sh[i]-1
      return sh, (sh, i)
    end
  end
end

function Base.length(iter::ShuffleIterator)
  return binomial(iter.p+iter.q, iter.p)
end

struct QuantumElem
  word::Vector{UInt8}
end

struct GenericQuantumElem
  hom::Vector{QuantumElem}
end

struct Derivation
  word::Vector{UInt8}
end

function degree(w::QuantumElem)
  return
end

function derivation(w::QuantumElem)
  return Derivation(w.word)
end

function Base.(*)(w1::QuantumElem, w2::QuantumElem)
  p = zero()
  for sh in shuffles(length(w1.word), length(w2.word))
    # TODO: replace cartan with bilinear form
    e = sum(cartan[w1.word[sh[i]], w1.word[sh[l]]] for k in 1:length(w1.word), l in length(w1.word)+1:length(w1.word)+length(w2.word))
    p += qparam^e * []
  end
end

A, v = LaurentPolynomialRing(ZZ, "v")
QF = fraction_field(A)

B = [
  4 -2
  -2 2
]

lam = [3,4]

f = [2,1,2,1]
p = [12, 24, 36, 12]

mat = [
  2 2 2 2 2 0 0 0 0 1 1 0
  4 4 4 4 4 2 2 0 0 0 0 0
  6 6 6 6 0 2 2 2 2 1 1 2
  0 0 0 0 0 2 2 4 4 0 0 0
]

mat = [
  2 2 2 2 2 2 0 0 0 0 0 0
  4 4 4 4 4 4 0 0 0 0 0 0
  6 6 6 6 6 6 0 0 0 0 0 0
  2 2 2 2 2 2 0 0 0 0 0 0
]

mat = [
2 2 2 2 2 0 0 0 0 1 0 1
4 4 4 4 0 2 2 2 2 0 0 0
6 6 6 6 0 2 2 2 2 1 2 1
0 0 0 0 4 2 2 2 2 0 0 0
]

C = [2 -1; -2 2]

function dirm()
  n, m = size(mat)
  dir_mat = zeros(Int, n, m)
  for j in 1:m
    len = [lam[f[a]] for a in 1:length(lam)]
    for i in 1:n
      dir_mat[i, j] = len[f[i]] - mat[i, j]
      for a in 1:length(lam)
        len[a] -= mat[i, j]*C[a, f[i]]
      end
    end
  end
  return dir_mat
end

function next(mat)
  row_sum = [sum(r) for r in eachrow(mat)]
  col_sum = [[sum(c[j] for j in 1:length(c) if f[j] == i) for i in 1:length(lam)] for c in eachcol(mat)]

  w = [copy(lam) for _ in 1:12]
  n, m = size(mat)
  weight_mat = [lam for i in 1:n, j in 1:m]

  for j in 1:m
    weight_mat[1, j] = copy(lam)
    weight_mat[1, j][f[1]] -= mat[1,j]
  end
  for i in 2:n, j in 1:m
    weight_mat[i, j] = copy(weight_mat[i-1, j])
    weight_mat[i, j][f[i]] -= mat[i,j]
  end
  
  dir_mat = zeros(Int, n, m)
  for j in 1:m
    len = [lam[f[a]] for a in 1:length(lam)]
    for i in 1:n
      dir_mat[i, j] = len[f[i]] - mat[i, j]
      for a in 1:length(lam)
        len[a] -= mat[i, j]*C[a, f[i]]
      end
    end
  end

  for i in n:-1:1, j in m-1:-1:1
    if mat[i,j] == 0
      continue
    end

    next_row = findnext(x -> x == f[i], f, i+1)
    if next_row == nothing
      continue
    end

    next_col = nothing
    for k in j+1:m
      if dir_mat[i, k] > 0 && col_sum[k][f[i]] > 0 && mat[next_row, k] > 0
        next_col = k
        break
      end
    end

    if next_col === nothing
      continue
    end
    
    skip = false
    for l in i+1:next_row-1
      if dir_mat[l, j]+C[f[l], f[i]] < 0
        a = findnext(x -> x > 0, dir_mat[l, :], j+1)
        if a == nothing
          skip = true
          break
        end
        mat[l, j] -= 1
        mat[l, a] += 1
      end
    end
    if skip
      continue
    end
    
    mat[i, j] -= 1
    mat[i, next_col] += 1
    mat[next_row, j] += 1
    mat[next_row, next_col] -= 1
    
    return mat
  end
end

function compute()
  q = 1
  w = [copy(lam) for _ in 1:12]
  for i in 1:nrows(mat)
    s = p[i]
    for j in ncols(mat):-1:2
      s -= mat[i,j]
      w[j][f[i]] -= mat[i,j]
      q *= v^-(s*mat[i,j]+s*dot(B[f[i], :], w[j]))
      #print(s*mat[i,j]+s*dot(B[f[i], :], w[j]), " ")
    end
    w[1][f[i]] -= mat[i,1]
    #println("")
  end
  
  return q
end