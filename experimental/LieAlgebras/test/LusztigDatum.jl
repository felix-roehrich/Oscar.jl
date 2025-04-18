@testset "LieAlgebras.LusztigDatum" begin
  @testset "Crystal.f!" begin
    R = root_system(:B, 3)
    w0 = UInt8[1, 2, 1, 3, 2, 3, 1, 2, 3]
    d = lusztig_datum(R, UInt8[1, 2, 1, 3, 2, 3, 1, 2, 3])
    Crystal.f!(d, [2, 1, 3, 2, 3, 1, 2], [2, 2, 4, 4, 4, 2, 2])
    LieAlgebras._state!(d, w0)
    @test d.datum == [0, 2, 0, 0, 2, 0, 0, 2, 0]
  end
end
