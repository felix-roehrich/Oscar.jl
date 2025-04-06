@testset "LieAlgebras.QuantumGroup" begin
  @testset "canonical_basis_elem" begin
    R = root_system(:B, 3)
    U = quantum_group(R)
    F = gens(U)

    f1 =
      F[1]^2 * F[3]^4 * F[9]^2 / (
        quantum_factorial(2, quantum_parameter(U, 1)) *
        quantum_factorial(4, quantum_parameter(U, 3)) *
        quantum_factorial(2, quantum_parameter(U, 9))
      )
    @test canonical_basis_elem(U, [2, 0, 4, 0, 0, 0, 0, 0, 2]) == f1
    @test canonical_basis_elem(U, [0, 1, 0, 0, 1, 0, 0, 1, 0]) == one(U)

    # canonical_basis_elem(U, [0, 1, 0, 0, 0, 0, 0, 1, 0])
  end
end
