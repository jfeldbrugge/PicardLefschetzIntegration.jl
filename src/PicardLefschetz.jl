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

##############################
# Setup the initial grid
##############################

function hypercube_triangulation_indices(n::Int)
    perms_list = collect(permutations(1:n))
    
    simplex_indices = Vector{index}()
    for σ in perms_list
        # Build the n+1 vertex indices for simplex S_σ
        indices = Vector{Int}(undef, n + 1)
        
        # v_0 = origin (all zeros) → index 0
        indices[1] = 1
        
        # Build intermediate vertices
        for k in 1:n
            # v_k has 1s at positions σ(1), ..., σ(k)
            # Convert to index: sum of 2^(σ[i]-1) for i in 1:k
            index_k = 0
            for i in 1:k
                # Position σ[i] corresponds to bit (σ[i] - 1)
                index_k += 2^(σ[i] - 1)
            end
            indices[k + 1] = index_k + 1
        end
        
        push!(simplex_indices, index(indices))
    end
    
    return simplex_indices
end

function initialGrid(min, max, pars::parameters)
    points = Vector{point}()
    
    # List the corners of the hypercube
    bounds = [[min[i], max[i]] for i in 1:pars.dim]
    for corner in Iterators.product(bounds...)
        push!(points, point(collect(corner)))
    end
    
    # Triangulate the hypercube
    simplices = hypercube_triangulation_indices(pars.dim)

    # Orient the simplices
    for i in eachindex(simplices)
        sim = deepcopy(simplices[i])
        vertices = map(p -> p.coord, points[sim.coord])
        vol = real(volume_simplex(stack(vertices, dims=1)))
        
        if vol < 0
            simplices[i].coord[end - 1], simplices[i].coord[end] = sim.coord[end], sim.coord[end - 1]
        end
    end

    # Subdivide the triangles
    thim = thimble(points, simplices)
    divide_rep(thim, pars)
    return thim
end

function volume_simplex(vertices)
    return det(transpose(vertices[2:end, :]) .- vertices[1, :]) / factorial(size(vertices, 1) - 1)
end

##############################
# Divide
##############################

function divide(thim, pars::parameters)
    dictionary = Dict()
    for i in eachindex(thim.simplices)
        sim_coord = thim.simplices[i].coord
        if thim.simplices[i].active
            edges = filter(v -> issorted(v), collect(permutations(sim_coord, 2)))
    
            # Check whether an edge of the simplex sim if it has already been subdivided
            subdivided = false
            for edge in edges
                if haskey(dictionary, edge) && subdivided == false
                    subdivided = true
                    thim.simplices[i].active = false
    
                    index_newPoint = dictionary[edge]
    
                    newsim1 = copy(sim_coord)
                    newsim2 = copy(sim_coord)
                    replace!(newsim1, edge[1] => index_newPoint)
                    replace!(newsim2, edge[2] => index_newPoint)
                    
                    push!(thim.simplices, index(newsim1))
                    push!(thim.simplices, index(newsim2))
                end
            end
    
            # Check the length of the edges and subdivide the simpex when the longest edge exceeds pars.δ
            if subdivided == false 
                edge_lengths = [norm(thim.points[edge[1]].coord - thim.points[edge[2]].coord) for edge in edges]
            
                index_longest = argmax(edge_lengths)
                if edge_lengths[index_longest] > pars.δ
                    thim.simplices[i].active = false
    
                    newPoint = point((thim.points[edges[index_longest][1]].coord + thim.points[edges[index_longest][2]].coord) ./ 2)
            
                    push!(thim.points, newPoint)
                    index_newPoint = length(thim.points)
    
                    newsim1 = copy(sim_coord)
                    newsim2 = copy(sim_coord)
                    replace!(newsim1, edges[index_longest][1] => index_newPoint)
                    replace!(newsim2, edges[index_longest][2] => index_newPoint)
                    
                    push!(thim.simplices, index(newsim1))
                    push!(thim.simplices, index(newsim2))
    
                    dictionary[edges[index_longest]] = index_newPoint
                end
            end
        end
    end
end

function divide_rep(thim::thimble, pars::parameters)
    n_old = length(thim.simplices)
    n_new = n_old + 1

    while n_old != n_new
        n_old = n_new
        divide(thim, pars)
        n_new = length(thim.simplices)
    end
    filter!(sim->sim.active, thim.simplices)
end

function remove_inactive_points(thim::thimble)
    dictionary = Dict()
    counter = 0
    for i in eachindex(thim.points)
        if thim.points[i].active
            counter += 1
            dictionary[i] = counter
        end
    end

    for i in eachindex(thim.simplices)
        coord = copy(thim.simplices[i].coord)
        thim.simplices[i].coord = [dictionary[j] for j in coord]
    end

    filter!(p -> p.active, thim.points)

    return thim
end

##############################
# Flow
##############################

function flow_point(S, p, pars::parameters)
    precision = 1e-10

    grad = zeros(ComplexF64, pars.dim)
    
    for i in 1:pars.dim
        ϵ = zeros(pars.dim)
        ϵ[i] = precision

        grad[i] = (S(p .+ ϵ) - S(p .- ϵ)) / (2. * precision)
    end

    return p - pars.ϵ * normalize(conj(im * grad))
end

function flow_points(S, thim::thimble, pars::parameters)
    for p in thim.points
        if p.active
            p.coord = flow_point(S, p.coord, pars)
        end
    end

    for i in eachindex(thim.simplices)
        if thim.simplices[i].active
            for v in thim.simplices[i].coord
                if real(im * S(thim.points[v].coord)) < pars.τ
                    thim.simplices[i].active = false
                    thim.points[v].active = false
                end
            end
        end
    end

    divide_rep(thim, pars)
end

function flow(S, thim::thimble, pars::parameters)
    for _ in 1:pars.N
        flow_points(S, thim, pars)
    end
    remove_inactive_points(thim)
end

##############################
# Integrate
##############################

function get_affine_map(vertices)
    v₀ = vertices[1, :]
    B = transpose(vertices[2:end, :]) .- v₀
    return v₀, B
end

function mapping(p, vertices)
    v₀, B = get_affine_map(vertices)
    return v₀ + B * p
end

function jacobian(vertices)
    return det(get_affine_map(vertices)[2])
end

function integrateSimplex(f, vertices, X, W)
    integrand(p) = jacobian(vertices) * f(mapping(p, vertices))
    return sum(W[i] * integrand(X[i,:]) for i in 1:length(W))
end

function PL_integrate(S, thim::thimble, pars::parameters)
    X, W = simplexquad(pars.n, pars.dim)

    points_r = map(p->p.coord, thim.points)
    simplices_r = map(sim->sim.coord, thim.simplices)

    sum = 0
    for sim in simplices_r
        if maximum(real.(im * map(S, points_r[sim]))) > pars.τ
            sum += integrateSimplex(p -> exp(im * S(p)), stack(points_r[sim], dims=1), X, W)
        end
    end

    return sum
end

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
