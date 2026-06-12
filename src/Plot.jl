
function triPlot(points, simplices)
    filter!(sim->sim.active, simplices)
    ps = [Polygon([Point2f(real.(points[i].coord)) for i in sim.coord]) for sim in simplices]
    
    f = Figure()
    ax = Axis(f[1, 1]; aspect=1)
    poly!(ax, ps, color = rand(RGBf, length(ps)))

    vertices = stack(real.(map(p -> p.coord, points)), dims=1)
    scatter!(ax, vertices[:, 1], vertices[:, 2])

    for i in 1:length(points)
        text!(ax, vertices[i, 1], vertices[i, 2] + 0.2, text = string(i))
    end

    return f
end

function linePlot(points, simplices)
    filter!(sim->sim.active, simplices)
    lines = reduce(vcat, [[Point2f(real(points[i].coord[1]), imag(points[i].coord[1])) for i in sim.coord] for sim in simplices])
    vertices = stack(map(p -> p.coord, points), dims=1)

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=1)
    linesegments!(ax, lines, linewidth=2, color=:red)
    scatter!(ax, real.(vertices[:]), imag.(vertices[:]))
    for i in 1:length(points)
        text!(ax, real(vertices[i]), imag(vertices[i]) + 0.2, text = string(i))
    end
    limits!(ax, -5, 5, -4, 4)
    fig    
end