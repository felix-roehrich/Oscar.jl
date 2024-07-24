GAP.Globals.LoadPackage(GAP.Obj("QuaGroup"))

#A, v = LaurentPolynomialRing(ZZ, "v")
#P, x = polynomial_ring(A, 6)
#M = [0 q^-1*x[1]*x[2]+x[3] q^-1*x[1]*x[3]; 0 0 q*x[2]*x[3]; 0 0 0]
#U, f = pbw_algebra(P, M, deglex(gens(P)))
#I = left_ideal(U, [f[1]^2, f[2]^2])

A, v = QQ["v"]
Qv = fraction_field(A)

function quantum_integer(n::Int, d::Int = 1)
  return (v^(d*n) - v^-(d*n))/(v - v^-1)
end


struct QuantumGroup
  pbw::PBWAlg
end

struct QuantumGroupPBW
  power::Int

  # root is the positive root associated to the PBW basis element
  root::Vector{Int}
end

struct QuantumGroupElem
end

function quantum_group(fam::Symbol, rk::Int)
  R = GAP.Globals.RootSystem(GAP.Obj(fam), rk)
  U = GAP.Globals.QuantizedUEA(R)
  g = GAP.Globals.GeneratorsOfAlgebra(U)
  
  poly_ring, x = polynomial_ring(Qv, length(R.positive_roots))
  rels = zero_matrix(poly_ring, length(g), length(g))
  for i in 1:length(g)
    for j in i+1:length(g)
      rels[i,j] = g[j]*g[i]
    end
  end
  
  return QuantumGroup(pbw_algebra(poly_ring, rels, deglex(poly_ring)))
end

function chevallay_gens(U::QuantumGroup)
end

function pbw_gens(U::QuantumGroup)
  return gens(U.pbw_algebra)
end

#=
function quantum_group(fam::Symbol, rk::Int)
  R = root_system(fam, rk)
  
  cm = cartan_matrix(R)
  symmetrization = zeros(ZZ, rank(R))
  for i in 1:rank(R)
    symmetrization[i] = 1
    for j in i+1:rank(R)
      if cm[i, j] != 0
        symmetrization[i] = cm[i, j] * cm[j, i]
        break
      end
    end
  end
  
  bil = diagonal_matrix(symmetrization)*cm
  
  # TODO: positive roots -> accessor
  poly, x = polynomial_ring(Qv, length(R.positive_roots))
  rels = zero_matrix(poly, length(R.positive_roots), length(R.positive_roots))
  
  for i in 1:length(R.positive_roots)
    for j in i+1:length(R.positive_roots)
      r = dot(ZZ.(R.positive_roots[i]), bil, ZZ.(R.positive_roots[j]))
      if r >= 0
        continue
      end
      
      ij = findfirst(x -> x == R.positive_roots[i]+R.positive_roots[j], R.positive_roots)
      rels[i, j] = x[ij] + v^(-r)*x[i]*x[j]
      if r == -2
        if dot(ZZ.(R.positive_roots[i]), bil, ZZ.(R.positive_roots[i])) == 1
          iij = findfirst(x -> x == R.positive_roots[i]+R.positive_roots[ij], R.positive_roots)
          rels[i, ij] = quantum_integer(2)*x[iij] + x[i]*x[ij]
        else
          ijj = findfirst(x -> x == R.positive_roots[ij]+R.positive_roots[j], R.positive_roots)
          rels[j, ij] = x[j]*x[ij] - quantum_integer(2)*x[ijj]
        end
      end
    end
  end
  
end
=#