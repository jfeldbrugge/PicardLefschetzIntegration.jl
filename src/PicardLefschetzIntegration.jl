module PicardLefschetzIntegration

    using LinearAlgebra, FastGaussQuadrature, Interpolations

    include("PL.jl")
    
    export parameters, initialGrid, flow, find_closest, PL

end
