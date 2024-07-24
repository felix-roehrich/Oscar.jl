#A, v = LaurentPolynomialRing(ZZ, "v")
A, v = cyclotomic_field(12)
B = [4 -2; -2 2]
C = [2 -1; -2 2]

GAP.Globals.LoadPackage(GAP.Obj("QuaGroup"))
R = GAP.Globals.RootSystem(GAP.Obj("B"), 2)
U = GAP.Globals.QuantizedUEA(R)
g = GAP.Globals.GeneratorsOfAlgebra(U)
V = GAP.Globals.HighestWeightModule(U, GAP.Obj([2,2]))
v0 = GAP.Globals.Basis(V)[1]

# number of partitions of n into m nonnegative numbers leq to bound
function partition(n::Int, m::Int, bound::Int = n)
  d = min(n, bound)
  
  if (m-1)*bound > n-d
    return sum(partition(n-i, m-1, min(n-i, bound)) for i in 0:d)
  elseif (m-1)*bound == n-d
    return 1
  else
    return 0
  end
end


function pow(f::GAP.GapObj, n::Int)
  if n == 0
    return GAP.Globals.One(U)
  end

  fam = GAP.Globals.ElementsFamily(GAP.Globals.FamilyObj(U))
  rep = GAP.Globals.List(GAP.Globals.ExtRepOfObj(f), GAP.Globals.ShallowCopy)
  rep[1][2] = n
  return GAP.Globals.ObjByExtRep(fam, rep)
end

lam = [2, 2]
f = [2,1,2,1]
pbw = [4, 1, 4, 1]
pbw_coeffs = [v^-2, v^-2, v^-2, v^-2, v^-12, v^-12, v^-12, v^-12, v^-12, A(1), A(1), A(1)]
pbw_coeffs = [v^-8, v^-8, v^-8, v^-8, v^-8, v^-8, A(1), A(1), A(1), A(1), A(1), A(1)]
p = [12, 24, 36, 12]

lex_max = [
2 2 2 2 2 2 0 0 0 0 0 0
4 4 4 4 4 4 0 0 0 0 0 0
6 6 6 6 0 0 2 2 2 2 2 2
0 0 0 0 0 0 4 4 4 0 0 0
]

lex_min = [
2 2 2 2 0 0 0 0 0 0 2 2
4 4 4 4 0 2 2 2 2 0 0 0
6 6 6 6 2 2 2 2 2 2 0 0
0 0 0 0 4 2 2 2 2 0 0 0
]

critical1 = [
2  2  2  2  2  2  0  0  0  0  0  0
4  4  4  4  0  2  2  2  2  0  0  0
6  6  6  6  0  0  2  2  2  2  2  2
0  0  0  0  4  2  2  2  2  0  0  0
]

critical2 = [
2  2  2  2  0  0  0  2  2  0  0  0
4  4  4  4  0  0  0  4  4  0  0  0
6  6  6  6  2  2  2  0  0  2  2  2
0  0  0  0  4  4  4  0  0  0  0  0
]

m_212_1 = [
2  2  2  2  2  2  0  0  0  0  0  0
4  4  4  4  4  4  0  0  0  0  0  0
6  6  6  6  6  6  0  0  0  0  0  0
0  0  0  0  0  0  2  2  2  2  2  2
]

m_121_2 = [
2 2 2 2 2 2 0 0 0 0 0 0
4 4 4 4 4 4 0 0 0 0 0 0
4 4 4 4 4 4 2 2 2 2 2 2
2 2 2 2 2 2 0 0 0 0 0 0
]

pbw_121_21_1 = [v^-8, v^-8, v^-8, A(1), A(1), A(1), A(1), A(1), A(1), A(1), A(1), A(1)]
m_121_21_1 = [
2 2 2 2 2 2 0 0 0 0 0 0
4 4 4 2 2 2 2 2 2 0 0 0
4 4 4 4 4 4 6 6 0 0 0 0
2 2 2 0 0 0 0 0 0 2 2 2
]

function direction_matrix(mat::Matrix{Int})
  n,m = size(mat)
  dir_mat = zeros(Int, n, m)
  
  for j in 1:m
    len = copy(lam)
    for i in 1:n
      for a in 1:length(lam)
        len[a] -= mat[i, j]*C[a, f[i]]
      end
      dir_mat[i, j] = mat[i,j] + len[f[i]]
    end
  end
  
  return dir_mat
end

function next_tensor(mat)
  n,m = size(mat)
  dir = direction_matrix(mat)
  
  next = copy(mat)
  dir_next = copy(dir)
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
      if dir[i, k] > 0 && mat[next_row, k] > 0
        next_col = k
        break
      end
    end
    if next_col == nothing
      continue
    end
    
    copy!(next, mat)
    copy!(dir_next, dir)
    
    dir_next[i, j] += 1
    dir_next[i, next_col] -= 1
    
    next[i, j] -= 1
    next[i, next_col] += 1
    next[next_row, j] += 1
    next[next_row, next_col] -= 1
    
    for ll in i+1:next_row-1
      dir_next[ll, j] += C[f[ll], f[i]]
      dir_next[ll, next_col] -= C[f[ll], f[i]]
    end
    
    c = next_col
    for l in i:next_row-1
      nr = findnext(x -> x == f[l], f, l+1)
      while true
        while c < m && dir_next[l, c] == 0
          c += 1
        end
        if dir_next[l, c] < 0
          @goto skip
        elseif c == m
          break
        end
        
        if nr == nothing || next[nr, c] == 0
          c += 1
          continue
        end
        
        k = c+1
        while k <= m && next[l, k] == 0
          k += 1
        end
        if k == m+1
          break
        end
        
        _min = min(next[nr, c], dir_next[l, c], next[l,k])
        
        dir_next[l, c] -= _min
        dir_next[l, k] += _min

        next[l, c] += _min
        next[l, k] -= _min
        next[nr, c] -= _min
        next[nr, k] += _min

        for ll in l+1:nr-1
          dir_next[ll, c] -= _min*C[f[ll], f[l]]
          dir_next[ll, k] += _min*C[f[ll], f[l]]
        end
      end
      c = j
    end
    return next
    
    @label skip
  end
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
    coeffs = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[2])
    q *= sum(l -> coeffs[1][l]*v^(l-1+coeffs[2]), 1:length(coeffs[1])) / pbw_coeffs[i]
  end
  
  return q
end

#=
partition_121_21_1 =[[2, 2, 2, 0, 0], [2, 2, 1, 1, 0], [2, 2, 1, 0, 1], [2, 2, 0, 2, 0], [2, 2, 0, 1, 1], [2, 2, 0, 0, 2], [2, 1, 2, 1, 0], [2, 1, 2, 0, 1], [2, 1, 1, 2, 0], [2, 1, 1, 1, 1], [2, 1, 1, 0, 2], [2, 1, 0, 2, 1], [2, 1, 0, 1, 2], [2, 0, 2, 2, 0], [2, 0, 2, 1, 1], [2, 0, 2, 0, 2], [2, 0, 1, 2, 1], [2, 0, 1, 1, 2], [2, 0, 0, 2, 2], [1, 2, 2, 1, 0], [1, 2, 2, 0, 1], [1, 2, 1, 2, 0], [1, 2, 1, 1, 1], [1, 2, 1, 0, 2], [1, 2, 0, 2, 1], [1, 2, 0, 1, 2], [1, 1, 2, 2, 0], [1, 1, 2, 1, 1], [1, 1, 2, 0, 2], [1, 1, 1, 2, 1], [1, 1, 1, 1, 2], [1, 1, 0, 2, 2], [1, 0, 2, 2, 1], [1, 0, 2, 1, 2], [1, 0, 1, 2, 2], [0, 2, 2, 2, 0], [0, 2, 2, 1, 1], [0, 2, 2, 0, 2], [0, 2, 1, 2, 1], [0, 2, 1, 1, 2], [0, 2, 0, 2, 2], [0, 1, 2, 2, 1], [0, 1, 2, 1, 2], [0, 1, 1, 2, 2], [0, 0, 2, 2, 2]]
=#