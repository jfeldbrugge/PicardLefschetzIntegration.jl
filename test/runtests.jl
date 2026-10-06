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

    # Integration threshold: a thimble flowed for S and reused for ω S. The default threshold pars.τ
    # skips most of the thimble; a lower threshold keeps the contributions that matter.
    pars = parameters(δ = 0.5, τ = -10., ϵ = 0.1, N = 40, n = 5, dim = 2)
    thim = initialGrid([-4, -4], [4, 4], pars)
    flow!(S₂, thim, pars)
    Q = quadrature(thim, pars)
    @test PL_integrate(S₂, Q, pars) == PL_integrate(S₂, Q, pars; τ = pars.τ)
    @test abs(PL_integrate(p -> 5 * S₂(p), Q, pars) - im * π / 5) > 1e-6
    @test abs(PL_integrate(p -> 5 * S₂(p), Q, pars; τ = -40.) - im * π / 5) < 1e-8
    @test PL_integrate(p -> 5 * S₂(p), thim, pars; τ = -40.) == PL_integrate(p -> 5 * S₂(p), Q, pars; τ = -40.)

    # The subdivision keeps the triangulation conforming: the edges that belong to a single triangle are
    # exactly the edges on the boundary of the initial grid
    function boundary_edges(thim)
        count = Dict{Tuple{Int, Int}, Int}()
        for sim in thim.simplices, (a, b) in ((1, 2), (2, 3), (1, 3))
            edge = minmax(sim.coord[a], sim.coord[b])
            count[edge] = get(count, edge, 0) + 1
        end
        return Base.count(==(1), values(count))
    end

    pars = parameters(δ = 0.35, τ = -20., ϵ = 0.1, N = 50, n = 5, dim = 2)
    S₂(p) = p[1]^2 + p[2]^2
    thim = initialGrid([-4, -4], [4, 4], pars)
    n_boundary = boundary_edges(thim)
    flow!(S₂, thim, pars)
    @test boundary_edges(thim) == n_boundary

    # Regression test: the subdivision used to loop forever at the 27th flow step of this exponent.
    # By the scaling x → 4^(-1/4) x the integral equals 2 I(c; ω = 4) for the controls c below.
    c = (0.3, -0.2, -0.5, 0.2, -0.4, 0.3, 0.1)
    S_c(p, c) = p[1]^4 - 6p[1]^2 * p[2]^2 + p[2]^4 + c[1] * p[1]^2 * p[2] + c[2] * p[1] * p[2]^2 +
                c[3] * p[1]^2 + c[4] * p[1] * p[2] + c[5] * p[2]^2 + c[6] * p[1] + c[7] * p[2]
    c_scaled = c .* (4^0.25, 4^0.25, 4^0.5, 4^0.5, 4^0.5, 4^0.75, 4^0.75)
    pars = parameters(δ = 0.35, τ = -10., ϵ = 0.1, N = 50, n = 8, dim = 2)
    thim_scaled = initialGrid([-4, -4], [4, 4], pars)
    flow!(p -> S_c(p, c_scaled), thim_scaled, pars)
    thim = initialGrid([-4, -4], [4, 4], pars)
    flow!(p -> S_c(p, c), thim, pars)
    I_scaled = PL_integrate(p -> S_c(p, c_scaled), thim_scaled, pars) / 2
    I = PL_integrate(p -> 4 * S_c(p, c), thim, pars)
    @test abs(I_scaled - I) / abs(I) < 1e-3

    # Test find_closest
    @test find_closest([1.], 0.8) == 1
    @test find_closest([1., 2., 3.], 0.8) == 1
    @test find_closest([1., 2., 3.], 1.8) == 2
    @test find_closest([1., 2., 3.], 2.2) == 2
    @test find_closest([1., 2., 3.], 3.8) == 3
end