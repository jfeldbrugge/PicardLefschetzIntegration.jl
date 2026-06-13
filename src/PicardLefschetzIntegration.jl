module PicardLefschetzIntegration

    using LinearAlgebra, FastGaussQuadrature, Interpolations, Combinatorics, SimplexQuad
    using CairoMakie, Makie.GeometryBasics

    include("PicardLefschetz.jl")
    include("Plot.jl")
    
    export parameters, point, index, thimble, find_closest, initialGrid, flow, PL
    export triPlot, linePlot
end
