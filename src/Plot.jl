
function triPlot(points, simplices)
    ps = [Polygon([Point2f(real.(points[i].coord)) for i in sim.coord]) for sim in simplices]
    
    f = Figure()
    ax = Axis(f[1, 1]; aspect=1)
    poly!(ax, ps, color = rand(RGBf, length(ps)))

    vertices = stack(real.(map(p -> p.coord, points)), dims=1)
    scatter!(ax, vertices[:, 1], vertices[:, 2])

    for i in 1:length(points)
        text!(ax, vertices[i, 1], vertices[i, 2] + 0.2, text = string(i))
        # text!(ax, position = (1, 1), text = string(i), textsize = 12)
        # vertices[i, 1], vertices[i, 2] + 0.2
    end

    return f
end