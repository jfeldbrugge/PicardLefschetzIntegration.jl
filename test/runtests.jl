using PicardLefschetzIntegration
using Test

@testset "PicardLefschetzIntegration.jl" begin
    # Fresnel integral
    pars = parameters(δ = 0.5, τ = -20., ϵ = 0.1, N = 40, n = 5, dim = 1)
    S(p) = p[1]^2
    thim = initialGrid([-4], [4], pars)
    flow!(S, thim, pars)
    @test abs.(PL_integrate(S, thim, pars) - (1+im) * sqrt(π / 2)) < 1e-7


    # Two-dimensional Fresnel integral
    pars = parameters(δ = 0.5, τ = -20., ϵ = 0.1, N = 40, n = 5, dim = 2)
    S(p) = p[1]^2 + p[2]^2
    thim = initialGrid([-4, -4], [4, 4], pars)
    flow!(S, thim, pars)
    @test abs.(PL_integrate(S, thim, pars) - im * π) < 1e-7

    # Test find_closest
    @test find_closest([1., 2., 3.], 1.8) == 2
end