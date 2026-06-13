
function triPlot(thim::thimble)
    filter!(sim->sim.active, thim.simplices)
    ps = [Polygon([Point2f(real.(thim.points[i].coord)) for i in sim.coord]) for sim in thim.simplices]
    
    f = Figure()
    ax = Axis(f[1, 1]; aspect=1)
    poly!(ax, ps, color = rand(RGBf, length(ps)))

    vertices = stack(real.(map(p -> p.coord, thim.points)), dims=1)
    scatter!(ax, vertices[:, 1], vertices[:, 2])

    for i in 1:length(thim.points)
        text!(ax, vertices[i, 1], vertices[i, 2] + 0.2, text = string(i))
    end

    return f
end

function linePlot(thim::thimble)
    filter!(sim->sim.active, thim.simplices)
    lines = reduce(vcat, [[Point2f(real(thim.points[i].coord[1]), imag(thim.points[i].coord[1])) for i in sim.coord] for sim in thim.simplices])
    vertices = stack(map(p -> p.coord, thim.points), dims=1)

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=1)
    linesegments!(ax, lines, linewidth=2, color=:red)
    scatter!(ax, real.(vertices[:]), imag.(vertices[:]))
    for i in 1:length(thim.points)
        text!(ax, real(vertices[i]), imag(vertices[i]) + 0.2, text = string(i))
    end
    limits!(ax, -5, 5, -4, 4)
    fig    
end