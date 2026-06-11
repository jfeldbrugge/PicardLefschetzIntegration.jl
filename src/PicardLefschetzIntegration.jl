module PicardLefschetzIntegration

    using LinearAlgebra, FastGaussQuadrature, Interpolations, Combinatorics, SimplexQuad
    using CairoMakie, Makie.GeometryBasics

    include("PL.jl")
    include("PL_sim.jl")
    include("Plot.jl")
    
    export parameters, initialGrid, flow, find_closest, PL, point, index
    export initialGrid_sim
    export triPlot
end
