@testset "LSFan" begin
  @testset "points_of_weight(LS::LSFan, w::WeightLatticeElem)" begin
    R = root_system(:A, 3)
    LS = ls_fan(R, [3, 0, 3])

    pts = points_of_weight(LS)
    @test length(pts) == 8

    R = root_system(:D, 4)
    LS = ls_fan(R, [1, 1, 1, 0])

    wt = WeightLatticeElem(R, [-1, 0, -1, 0])
    pts = points_of_weight(LS, wt)
    @test length(pts) == 9
    @test length(unique(pts)) == 9
    @test all(p -> weight(p) == wt, pts) == true
  end
end

@testset "LSFanElem" begin
  @testset "ei!" begin
    R = root_system(:B, 2)
    W = weyl_group(R)
    s = gens(W)
    LS = ls_fan(2*weyl_vector(R))
    
    a = LS((1//2, s[2]*s[1]), (1//2, s[2]))
    @test ei!(a, 2) == LS((1//2, s[2]*s[1]), (1//2, one(W)))
  end

  @testset "sequence(a::LSFanElem)" begin
    @testset "B2" begin
      R = root_system(:B, 2)
      W = weyl_group(R)
      s = gens(W)
      LS = ls_fan(2*weyl_vector(R))

      a = LS((1//2, longest_element(W)), (1//2, one(W)))
      @test sequence(a) == [(2, 1), (1, 2), (2, 3), (1, 1)]
      @test sequence(a, [2, 1, 2, 1]) == [(1, 1), (2, 3), (1, 2), (2, 1)]

      a = LS((1//4, s[1]*s[2]*s[1]), (5//12, s[2]*s[1]), (1//3, s[1]))
      @test sequence(a) == [(1, 1), (2, 4), (1, 2)]

      a = LS((1//3, s[2]*s[1]*s[2]), (5//12, s[1]*s[2]), (1//4, s[2]))
      @test sequence(a) == [(2, 2), (1, 3), (2, 2)]
    end
  end

  @testset "leq_points_of_weight(a::LSFanElem)" begin
    R = root_system(:B, 2)
    W = weyl_group(R)
    s = gens(W)
    LS = ls_fan(2*weyl_vector(R))

    a = LS((1//2, longest_element(W)), (1//2, one(W)))
    pts = leq_points_of_weight(a)
    @test length(pts) == 5

    R = root_system(:A, 3)
    W = weyl_group(R)
    s = gens(W)

    a = LS((1//3, s[1]*s[2]*s[3]*s[2]*s[1]), (1//6, s[3]*s[2]*s[1]), (1//6, s[2]*s[1]), (1//3, one(W)))
    pts = leq_points_of_weight(a)
    @test length(pts) == 8
  end
end

@testset "TensorIterator" begin
  @testset "B2" begin
    R = root_system(:B, 2)
    W = weyl_group(R)
    s = gens(W)
    LS = ls_fan(2*weyl_vector(R))

    a = LS((1//2, longest_element(W)), (1//2, one(W)))
    b = LS((1//4, s[1]*s[2]*s[1]), (5//12, s[2]*s[1]), (1//3, s[1]))

    iter = Oscar.QuantumGroups.TensorIterator(a, b; rdec=[1, 2, 1, 2])
    mat, _ = iterate(iter)
    mat == [
    2 2 2 2 2 2 0 0 0 0 0 0
    4 4 4 2 2 2 2 2 2 0 0 0
    4 4 4 4 4 4 6 6 0 0 0 0
    2 2 2 0 0 0 0 0 0 2 2 2
    ]
  end
end

@testset "tensor_coefficient" begin
  @testset "A3" begin
    R = root_system(:A, 3)
    W = weyl_group(R)
    s = gens(W)
    LS = ls_fan(R, [3, 0, 3])

    a = LS((1//3, s[1]*s[2]*s[3]*s[2]*s[1]), (1//6, s[3]*s[2]*s[1]), (1//6, s[2]*s[1]), (1//3, one(W)))
    pts = leq_points_of_weight(a)

    iter = Oscar.QuantumGroups.TensorIterator(pts[7], pts[8])
    tensor_coefficient(pts[7], pts[8]) == 3

    LS = ls_fan(R, 2*weyl_vector(R))
    pts = points_of_weight(LS, [0, 0, 0])
    a = LS((1//2, longest_element(W)), (1//2, one(W)))
    @test tensor_coefficient(a, pts[7]) == (3985980, 4)
  end

  @testset "fi!" begin
    R = root_system(:A, 3)
    W = weyl_group(R)
    s = gens(W)
    LS = ls_fan(R, 2*weyl_vector(R))

    a = LS()
    @test fi!(a, 1, 1) == LS((1//2, s[1]), (1//2, one(W)))
    @test fi!(a, 2, 2) == LS((1//2, s[2]*s[1]), (1//2, one(W)))
    @test fi!(a, 1, 2) == LS((1//2, s[1]*s[2]*s[1]), (1//2, s[1]))
    @test fi!(a, 3, 2) == LS((1//3, s[3]*s[1]*s[2]*s[1]), (1//6, s[1]*s[2]*s[1]), (1//2, s[1]))
    @test fi!(a, 2, 1) == LS((1//4, s[2]*s[3]*s[1]*s[2]*s[1]), (1//12, s[3]*s[1]*s[2]*s[1]), (1//6, s[1]*s[2]*s[1]), (1//2, s[1]))
    @test fi!(a, 2, 1) == LS((1//4, s[2]*s[3]*s[1]*s[2]*s[1]), (1//12, s[3]*s[1]*s[2]*s[1]), (1//6, s[1]*s[2]*s[1]), (1//4, s[2]*s[1]), (1//4, s[1]))

    # [0, 1, 2, 2, 2, 1] with fi!
    # [0, 2, 2, 2, 2, 1] with ei!

    # Oscar.QuantumGroups.fi!(LS(), [1, 2, 3, 1, 2, 1], [1, 3, 2, 1, 0, 0]) == 1//4*e(s1 * s2 * s3 * s1) + 1//4*e(s2 * s3 * s1) + 1//2*e(s3)
  end

  @testset "ei!" begin
    R = root_system(:A, 3)
    LS = ls_fan(R, 2*weyl_vector(R))

    a = Oscar.QuantumGroups.apply_sequence(LS, [0, 2, 2, 2, 2, 1])
    @test ei!(a, 2) == LS((1//4, s[2]*s[3]*s[1]*s[2]*s[1]), (1//12, s[3]*s[1]*s[2]*s[1]), (1//6, s[1]*s[2]*s[1]), (1//2, s[1]))

    # 1//4*e(s2 * s3 * s1 * s2 * s1) + 1//4*e(s3 * s1 * s2 * s1) + 1//4*e(s2 * s1) + 1//4*e(s1)
  end
end
