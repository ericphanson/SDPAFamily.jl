@testset "TemporaryDirectory" begin

    optimizer = SDPAFamily.Optimizer{Float64}()
    tmp = MOI.get(optimizer, SDPAFamily.TemporaryDirectory())
    @test tmp == optimizer.tempdir
    @test isdir(tmp)

    tmp_arg = mktempdir()
    optimizer = SDPAFamily.Optimizer{Float64}(TemporaryDirectory=tmp_arg)
    @test tmp_arg == MOI.get(optimizer, SDPAFamily.TemporaryDirectory())
    

    new_tmp = mktempdir()
    MOI.set(optimizer, SDPAFamily.TemporaryDirectory(), new_tmp)
    @test new_tmp == MOI.get(optimizer, SDPAFamily.TemporaryDirectory())

    # check that it survives `empty`
    MOI.empty!(optimizer)
    @test new_tmp == MOI.get(optimizer, SDPAFamily.TemporaryDirectory())

end
