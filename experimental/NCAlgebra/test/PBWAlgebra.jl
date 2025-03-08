@testset "PBWAlgebra" begin
  using NCAlgebra: _linear_index

  @testset "constructor" begin
    A, q = laurent_polynomial_ring(ZZ, "q")
    QA = fraction_field(A)
    R, F = polynomial_ring(QA, ("F" => 1:9))

    rels = [
      q^-2 * F[1] * F[2],
      q^-2 * F[1] * F[3],
      F[1] * F[4] + (q + q^-1) * F[2],
      F[1] * F[5] + (-q^2 + q^-2) * F[3] * F[4],
      F[1] * F[6] + (q + q^-1) * F[3],
      q^2 * F[1] * F[7] + F[4],
      q^2 * F[1] * F[8] + F[6],
      F[1] * F[9],
      q^-2 * F[2] * F[3],
      q^-2 * F[2] * F[4],
      q^-2 * F[2] * F[5] + (q - 2 * q^-1 + q^-3) * F[3] * F[4]^2,
      F[2] * F[6] + (-q^2 + q^-2) * F[3] * F[4],
      F[2] * F[7] + (-q + q^-1) // (q^2 + 1) * F[4]^2,
      q^2 * F[2] * F[8] + (-q^4 + 1) * F[3] * F[7] + (-q + q^-1) * F[4] * F[6] + (-q^2 - 1 + q^-2) * F[5],
      q^2 * F[2] * F[9] + F[3],
      F[3] * F[4],
      q^-2 * F[3] * F[5],
      q^-2 * F[3] * F[6],
      q^2 * F[3] * F[7] + F[5],
      F[3] * F[8] + (-q + q^-1) // (q^2 + 1) * F[6]^2,
      q^-2 * F[3] * F[9],
      q^-2 * F[4] * F[5],
      F[4] * F[6] + (q + q^-1) * F[5],
      q^-2 * F[4] * F[7],
      F[4] * F[8] + (-q^2 + q^-2) * F[6] * F[7],
      q^2 * F[4] * F[9] + F[6],
      q^-2 * F[5] * F[6],
      q^-2 * F[5] * F[7],
      q^-2 * F[5] * F[8] + (q - 2 * q^-1 + q^-3) * F[6]^2 * F[7],
      F[5] * F[9] + (-q + q^-1) // (q^2 + 1) * F[6]^2,
      F[6] * F[7],
      q^-2 * F[6] * F[8],
      q^-2 * F[6] * F[9],
      q^-2 * F[7] * F[8],
      q^2 * F[7] * F[9] + F[8],
      q^-2 * F[8] * F[9],
    ]

    for i in 1:9
      for j in i+1:9
        ind = _linear_index(U, i, j)
        @test (F[j] * F[i]).poly == rels[ind]
      end
    end
    @test (F[2] + F[3])^2 == F[2]^2 + (1 + q^-2) * F[2] * F[3] + F[3]^2
    @test F[1] * F[7] * F[1] == (q^3 + q) * F[1]^(2) * F[7] + F[1] * F[4]
    @test (q^2 + q^-2) * F[4]^2 * F[1]
    @test F[4]^2 * F[1]^2
    @test (F[7] + F[9]) * (F[1] + F[7]) == q^2 * F[1] * F[7] + F[1] * F[9] + F[4] + F[7]^2 + q^2 * F[7] * F[9] + F[8]
  end
end
