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

    # Quartic integral ∫exp(i t⁴)dt = 2Γ(5/4)exp(iπ/8)
    pars = parameters(δ = 0.25, τ = -20., ϵ = 0.1, N = 40, n = 8, dim = 1)
    S₄(p) = p[1]^4
    thim = initialGrid([-4], [4], pars)
    flow!(S₄, thim, pars)
    @test abs(PL_integrate(S₄, thim, pars) - 2 * 0.9064024770554771 * exp(im * π / 8)) < 1e-7

    # Reuse of the quadrature rule of a thimble
    pars = parameters(δ = 0.5, τ = -20., ϵ = 0.1, N = 40, n = 5, dim = 2)
    S₂(p) = p[1]^2 + p[2]^2
    thim = initialGrid([-4, -4], [4, 4], pars)
    flow!(S₂, thim, pars)
    Q = quadrature(thim, pars)
    @test length(Q.weights) == size(Q.nodes, 2) == length(thim.simplices) * Q.npts
    @test PL_integrate(S₂, Q, pars) == PL_integrate(S₂, thim, pars)
    for ω in (1., 2., 5.)
        @test abs(PL_integrate(p -> ω * S₂(p), Q, pars) - im * π / ω) < 1e-7
    end

    # Test find_closest
    @test find_closest([1.], 0.8) == 1
    @test find_closest([1., 2., 3.], 0.8) == 1
    @test find_closest([1., 2., 3.], 1.8) == 2
    @test find_closest([1., 2., 3.], 2.2) == 2
    @test find_closest([1., 2., 3.], 3.8) == 3
end