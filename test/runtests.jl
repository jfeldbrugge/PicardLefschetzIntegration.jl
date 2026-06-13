using PicardLefschetzIntegration
using Test

@testset "PicardLefschetzIntegration.jl" begin
    pars = parameters(δ = 0.5, τ = -10., ϵ = 0.1, N = 20, n = 5, dim = 1)

    S(p) = p[1]^2
    thim = initialGrid([-4], [4], pars)
    flow!(S, thim, pars)
    @test abs.(PL_integrate(S, thim, pars) - (1+im) * sqrt(π / 2)) < 1e-4

    @test find_closest([1., 2., 3.], 1.8) == 2
end