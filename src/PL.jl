# General functions
struct parameters
    δ           #The maximal size of the simplices
    τ           #Threshold on the real part of the exponent
    ϵ           #The steps in the downward flow
    N           #Number of steps in the downward flow
    n           #Order of quadrature integration
    dim         #Dimension of the integral

    function parameters(;δ, τ, ϵ, N, n, dim)
        return new(δ, τ, ϵ, N, n, dim)
    end
end

mutable struct point
    coord::Vector{Complex}
    active
    point(coord) = new(coord, true)
end

mutable struct index
    coord::Vector{Int}
    active
    index(coord) = new(coord, true)
end

mutable struct thimble
    points
    simplices
end

function normalize(p)
    if norm(p) > 1.
        return p / norm(p)
    else
        return p
    end
end

function find_closest(A::AbstractArray{T}, b::T) where {T<:Real}
    if length(A) <= 1
        return firstindex(A)
    end

    i = searchsortedfirst(A, b)

    if i == firstindex(A)
        return i
    elseif i > lastindex(A)
        return lastindex(A)
    else
        prev_dist = b - A[i-1]
        next_dist = A[i] - b

        if prev_dist < next_dist
            return i - 1
        else
            return i
        end
    end
end

# Setup the initial grid
function initialGrid(min, max, pars::parameters)
    points = [
        point([min[1], min[2]]), 
        point([min[1], max[2]]), 
        point([max[1], max[2]]),
        point([max[1], min[2]])] 
    simplices = [index([1, 2, 3, 4])]
    @show points

    divide_rep(points, simplices, pars)
    return (points, simplices)
end

function divide(points, simplices, pars::parameters)
    for i in eachindex(simplices)
        sim = simplices[i]
        if sim.active
            v1, v2, v3, v4 = sim.coord[1], sim.coord[2], sim.coord[3], sim.coord[4]
            V1, V2, V3, V4 = points[v1].coord, points[v2].coord, points[v3].coord, points[v4].coord

            if norm(V1 - V2) > pars.δ || norm(V2 - V3) > pars.δ || norm(V3 - V4) > pars.δ || norm(V4 - V1) > pars.δ
                l = length(points)
                append!(points, [
                    point((V1 + V2)/ 2.),            # l + 1
                    point((V2 + V3)/ 2.),            # l + 2
                    point((V3 + V4)/ 2.),            # l + 3
                    point((V4 + V1)/ 2.),            # l + 4
                    point((V1 + V2 + V3 + V4)/ 4.)]) # l + 5
                simplices[i].active = false
                append!(simplices, [
                    index([v1, l + 1, l + 5, l + 4]),
                    index([l + 1, v2, l + 2, l + 5]),
                    index([l + 5, l + 2, v3, l + 3]),
                    index([l + 4, l + 5, l + 3, v4])])
            end 
        end
    end
    filter!(sim->sim.active, simplices)
end

function divide_rep(points, simplices, pars::parameters)
    n_old = length(simplices)
    n_new = n_old + 1

    while n_old != n_new
        n_old = n_new
        divide(points, simplices, pars)
        n_new = length(simplices)
    end
end

# Flow
function flow_point(S, p, pars::parameters)
    precision = 1e-10
    ϵ1, ϵ2 = [precision, 0], [0, precision]
    
    der = [S(p + ϵ1) - S(p - ϵ1), S(p + ϵ2) - S(p - ϵ2)] / (2. * precision)
    
    return p - pars.ϵ * normalize(conj(im * der))
end

function flow_points(S, points, simplices, pars::parameters)
    for p in points
        if p.active
            p.coord = flow_point(S, p.coord, pars)
        end
    end

    for i in eachindex(simplices)
        if simplices[i].active
            for v in simplices[i].coord
                if real(im * S(points[v].coord)) < pars.τ
                    simplices[i].active = false
                    points[v].active = false
                end
            end
        end
    end

    divide_rep(points, simplices, pars)
end

function flow(S, points, simplices, pars::parameters)
    for i in 1:pars.N
        flow_points(S, points, simplices, pars)
    end

    nothing
end

# Integrate
function mapping(p, p1, p2, p3, p4)
    return (p1 .* (1. - p[1]) * (1. - p[2]) + 
            p2 .* (1. + p[1]) * (1. - p[2]) + 
            p3 .* (1. + p[1]) * (1. + p[2]) + 
            p4 .* (1. - p[1]) * (1. + p[2])) / 4.
end

function jacobian(p, p1, p2, p3, p4)
    A = +(p1[1] - p3[1]) * (p2[2] - p4[2]) - (p1[2] - p3[2]) * (p2[1] - p4[1])
    B = -(p1[1] - p2[1]) * (p3[2] - p4[2]) + (p1[2] - p2[2]) * (p3[1] - p4[1])
    C = +(p2[1] - p3[1]) * (p1[2] - p4[2]) - (p2[2] - p3[2]) * (p1[1] - p4[1])
    return (A + B * p[1] + C * p[2]) / 8
end

function IntegrateQuad(integrand, quad, n)
    latice, weights = gausslegendre(n)
    return IntegrateQuad(integrand, quad, n, latice, weights)
end

function IntegrateQuad(integrand, quad, n, latice, weights)
    sum = 0
    for i=1:n, j=1:n
        sum = sum + weights[i] * weights[j] * integrand([latice[i], latice[j]], quad[1], quad[2], quad[3], quad[4])
    end
    return sum
end

function PL(S, points, simplices, n)
    latice, weights = gausslegendre(n)

    points_r = map(p->p.coord, points)
    simplices_r = map(sim->sim.coord, simplices)

    function integrand(p, p1, p2, p3, p4)
        return jacobian(p, p1, p2, p3, p4) * exp(im * S(mapping(p, p1, p2, p3, p4)))    
    end

    sum = 0
    for i in eachindex(simplices)
        sum = sum + IntegrateQuad(integrand, points_r[simplices_r[i]], n, latice, weights)
    end
    return sum
end

function PL(S, thimble, pars::parameters)
    lattice, weights = gausslegendre(pars.n)

    function integrand(p, p1, p2, p3, p4)
        return jacobian(p, p1, p2, p3, p4) * exp(im * S(mapping(p, p1, p2, p3, p4)))    
    end

    sum = 0
    for quad in thimble
        sum = sum + IntegrateQuad(integrand, quad, pars.n, lattice, weights)
    end
    return sum
end
