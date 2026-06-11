
function triPlot(points, simplices)
    ps = [Polygon([Point2f(real.(points[i].coord)) for i in sim.coord]) for sim in simplices]
    
    f = Figure()
    Axis(f[1, 1]; backgroundcolor = :gray15)
    poly!(ps, color = rand(RGBf, length(ps)))
    
    return f
end