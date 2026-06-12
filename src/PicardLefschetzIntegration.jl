module PicardLefschetzIntegration

    using LinearAlgebra, FastGaussQuadrature, Interpolations, Combinatorics, SimplexQuad
    using CairoMakie, Makie.GeometryBasics

    include("PL.jl")
    include("PL_sim.jl")
    include("Plot.jl")
    
    export parameters, initialGrid, flow, find_closest, PL, point, index, thimble
    export initialGrid_sim, flow_sim, PL_sim
    export triPlot, linePlot
end
