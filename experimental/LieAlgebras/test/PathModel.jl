@testset "LieAlgebras.PathModel" begin
  #GAP.Packages.load("QuaGroup")

  @testset "LSPathModel" begin
    @testset "dominant_path" begin
      R = root_system(:A, 3)

      wt = WeightLatticeElem(R, [1, 1, 1])
      P = LSPathModel(wt)

      p0 = dominant_path(P)
      @test weight(p0) == wt
      @test all(PathModel.eps(p0, i) == zero(QQ) for i in 1:rank(R))
      @test all(PathModel.phi(p0, i) == wt[i] for i in 1:rank(R))
    end

    @testset "PathModel.e" begin
      R = root_system(:A, 3)
      p = PathModel.f(p0, 1)

      for i in 1:rank(R)
        p2 = PathModel.f(p, i)
        @test PathModel.eps(p2, i) == PathModel.eps(p, i) + 1
        @test PathModel.phi(p2, i) == PathModel.phi(p, i) - 1
        @test PathModel.weight(p2) == PathModel.weight(p) - simple_root(R, i)
      end
    end

    @testset "PathModel.f" begin
      R = root_system(:A, 3)
      p = PathModel.f(p0, 1)

      for i in 1:rank(R)
        p2 = PathModel.e(p, i)
        @test PathModel.eps(p2, i) == PathModel.eps(p, i) - 1
        @test PathModel.phi(p2, i) == PathModel.phi(p, i) + 1
        @test PathModel.weight(p2) == PathModel.weight(p) + simple_root(R, i)
      end
    end
  end
end
