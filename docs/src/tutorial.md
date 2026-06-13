```@meta
CurrentModule = PicardLefschetzIntegration
```

# Tutorial

## Fresnel integral
```math
I = \int_{-\infty}^\infty e^{i x^2}\mathrm{d}x
```

```@example tutorial1
using PicardLefschetzIntegration

pars = parameters(δ = 0.5, τ = -10., ϵ = 0.1, N = 20, n = 5, dim = 1)

thim = initialGrid([-4], [4], pars)
S(p) = p[1]^2
flow(S, thim, pars)
linePlot(thim)
```

```@example tutorial1
PL(S, thim, pars)
```

## Pearcey integral

```@example tutorial2
using PicardLefschetzIntegration, ProgressMeter, Base.Threads, CairoMakie

pars = parameters(δ = 0.5, τ = -10., ϵ = 0.1, N = 50, n = 5, dim = 1)

S(p, x₁, x₂, ω) = ω * (p[1]^4 + x₂ * p[1]^2 + x₁ * p[1])

aRange_L = range(-2., 2., 8)
bRange_L = range(-2., 2., 8)

thimbles = Array{thimble}(undef, length(aRange_L), length(bRange_L));
let ω = 1
    @showprogress Threads.@threads for (i, j) in collect(Iterators.product(eachindex(aRange_L),eachindex(bRange_L)))
        x₁, x₂ = aRange_L[i], bRange_L[j]
        
        thim = initialGrid([-4], [4], pars)
        flow(p -> S(p, x₁, x₂, ω), thim, pars)
        thimbles[i, j] = thim
    end
end

map(t->length(t.simplices), thimbles)
```

```@example tutorial2
aRange = range(-2., 2., 200)
bRange = range(-2., 2., 200)

let ω = 20.
    data = zeros(Complex, length(aRange), length(bRange))

    @showprogress Threads.@threads for (i, j) in collect(Iterators.product(eachindex(aRange),eachindex(bRange)))
        x₁, x₂ = aRange[i], bRange[j]
        ii, jj = find_closest(aRange_L, x₁), find_closest(bRange_L, x₂)
        
        data[i, j] = PL(p -> S(p, x₁, x₂, ω), thimbles[ii, jj], pars)
    end

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1) 
    heatmap!(ax, aRange, bRange, abs.(data).^2, interpolate=true)
    fig
end
```

## Elliptic integral
```@example tutorial3
using PicardLefschetzIntegration, ProgressMeter, Base.Threads, CairoMakie

pars = parameters(δ = 0.5, τ = -10., ϵ = 0.1, N = 20, n = 5, dim = 2)

S(p, x₁, x₂, x₃, ω) = ω * (p[1]^3 - 3. * p[1] * p[2]^2 - x₃ * (p[1]^2 + p[2]^2) - x₂ * p[2] - x₁ * p[1])

aRange_L = range(-2., 2., 8)
bRange_L = range(-2., 2., 8)

thimbles = Array{thimble}(undef, length(aRange_L), length(bRange_L));
let ω = 1, x₃ = 1
    @showprogress Threads.@threads for (i, j) in collect(Iterators.product(eachindex(aRange_L),eachindex(bRange_L)))
        x₁, x₂ = aRange_L[i], bRange_L[j]
        
        thim = initialGrid([-4, -4], [4, 4], pars)
        flow(p -> S(p, x₁, x₂, x₃, ω), thim, pars)
        thimbles[i, j] = thim
    end
end

map(t->length(t.simplices), thimbles)
```

```@example tutorial3
aRange = range(-2., 2., 100)
bRange = range(-2., 2., 100)

let ω = 20, x₃ = 1
    data = zeros(Complex, length(aRange), length(bRange))

    @showprogress Threads.@threads for (i, j) in collect(Iterators.product(eachindex(aRange),eachindex(bRange)))
        x₁, x₂ = aRange[i], bRange[j]
        ii, jj = find_closest(aRange_L, x₁), find_closest(bRange_L, x₂)
        
        data[i, j] = PL(p -> S(p, x₁, x₂, x₃, ω), thimbles[ii, jj], pars)
    end

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1) 
    heatmap!(ax, aRange, bRange, abs.(data).^2, interpolate=true)
    fig
end
```