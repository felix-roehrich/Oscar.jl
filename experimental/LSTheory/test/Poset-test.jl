@testset "Poset" begin
  using LSTheory: rel_index

  # Weyl group for A2 (covering relations)
  a2_cov = [
    0 1 1 0 0 0
    0 0 0 2 1 0
    0 0 0 1 2 0
    0 0 0 0 0 1
    0 0 0 0 0 1
    0 0 0 0 0 0
  ]

  # Weyl group for A2 (general relations)
  a2_gen = BitMatrix([
    0 1 1 1 1 1
    0 0 0 1 1 1
    0 0 0 1 1 1
    0 0 0 0 0 1
    0 0 0 0 0 1
    0 0 0 0 0 0
  ])

  # covering relations, B2 adjoint rep
  b2_adj_cov = [
    0 1 1 0 0 0 0 0
    0 0 0 3 1 0 0 0
    0 0 0 1 2 0 0 0
    0 0 0 0 0 2 1 0
    0 0 0 0 0 1 3 0
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 0
  ]

  # Weyl group for B2 (general relations)
  b2_adj_gen = BitMatrix(
    [
      0 1 1 1 1 1 1 1
      0 0 0 1 1 1 1 1
      0 0 0 1 1 1 1 1
      0 0 0 0 0 1 1 1
      0 0 0 0 0 1 1 1
      0 0 0 0 0 0 0 1
      0 0 0 0 0 0 0 1
      0 0 0 0 0 0 0 0
    ],
  )

  @testset "<(x::PosetElem, y::PosetElem)" begin
    for _ in 1:10
      ps = poset(a2_cov)
      for _ in 1:5
        sz = length(ps.elems)
        x = rand(1:6)
        y = rand(1:6)

        @test a2_gen[x, y] == (ps(x) < ps(y))
        @test all(
          ps.set[rel_index(ps, i, j)] == false ||
          ps.rel[rel_index(ps, i, j)] == a2_gen[i, j] for i in 1:sz for j in (i + 1):sz
        )
      end
    end
  end
end
