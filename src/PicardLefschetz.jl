# General functions

"""
    parameters

Contains the parameters of the Picard Lefschetz integration method.

    dim         The dimension of the integral.
    δ           The maximal size of the edges of the simplices.
    τ           Threshold on the real part of the exponent.
    ϵ           The steps size in the downward flow.
    N           The number of steps in the downward flow.
    n           The order of the simplex quadrature integration scheme.
"""
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

"""
    point

A point with two boolian attributes. The first determines whether the point is active. The second determines whether the point flows or stays fixed.

    coord::Vector{Complex}  The coordinates of the point.
    active                  When the point is active, it is still part of the Picard Lefschetz flow.
    frozen                  When a point is fozen, it does not flow. This is particularly relevant for points on the boundary of the integration domain.
"""
mutable struct point
    coord::Vector{Complex}
    active
    frozen
    point(coord, fixed = false) = new(coord, true, fixed)
end

"""
    index

The indices that span a simplex. 
    coord::Vector{Int}  The indices of the simplex refering to vertices in the vector points.
    active              When the simplex is active, it is still part of the Picard Lefschetz flow.
"""
mutable struct index
    coord::Vector{Int}
    active
    index(coord) = new(coord, true)
end

"""
    thimble

The thimble consisting of a set of points and a set of simplices.
"""
mutable struct thimble
    points
    simplices
end


"""
    normalize(p, thres = 1.)

Normalize a vector when the norm exceeds the threshold.
"""
function normalize(p, thres = 1.)
    if norm(p) > thres
        return p / norm(p)
    else
        return p
    end
end

"""
    find_closest(A::AbstractArray{T}, b::T) where {T<:Real}

Find indeces of closest point in a vector.
"""
# function find_closest(A::AbstractArray{T}, b::T) where {T<:Real}
function find_closest(A, b)
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

"""
    hypercube_triangulation_indices(dim::Int)

Triangulate a dim dimensional hypercube and output the simplices.
"""
function hypercube_triangulation_indices(dim::Int)
    perms_list = collect(permutations(1:dim))
    
    simplex_indices = Vector{index}()
    for σ in perms_list
        # Build the n+1 vertex indices for simplex S_σ
        indices = Vector{Int}(undef, dim + 1)
        
        # v_0 = origin (all zeros) → index 0
        indices[1] = 1
        
        # Build intermediate vertices
        for k in 1:dim
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

"""
    initialGrid(min, max, pars::parameters, frozen = false)

Initialize a cubic grid with the lower left corner min and the upper right corner max. When forzen = true, the boundary points do not flow.
"""
function initialGrid(min, max, pars::parameters, frozen = false)
    @assert length(min) == pars.dim "The dimension of min should be the same as pars.dim"
    @assert length(max) == pars.dim "The dimension of max should be the same as pars.dim"

    points = Vector{point}()
    
    # List the corners of the hypercube
    bounds = [[min[i], max[i]] for i in 1:pars.dim]
    for corner in Iterators.product(bounds...)
        push!(points, point(collect(corner), frozen))
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

"""
    volume_simplex(vertices)

The volume of a simplex.
"""
function volume_simplex(vertices)
    return det(transpose(vertices[2:end, :]) .- vertices[1, :]) / factorial(size(vertices, 1) - 1)
end

##############################
# Divide
##############################

"""
    divide(thim, pars::parameters)

Subdivide the simplices that have an edge longer than the threshold pars.δ.
"""
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

"""
    divide_rep(thim::thimble, pars::parameters)

Repeat subdivision till all edges of the simplices are shorter than par.δ.
"""
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

"""
    remove_inactive_points(thim::thimble)

Remove inactive points.
"""
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

"""
    flow_point(S, p, pars::parameters)

Push a point one step in the downward flow.
"""
function flow_point(S, p, pars::parameters)
    @assert length(p) == pars.dim "The dimension of p should be the same as pars.dim"
    precision = 1e-10

    grad = zeros(ComplexF64, pars.dim)
    
    for i in 1:pars.dim
        ϵ = zeros(pars.dim)
        ϵ[i] = precision

        grad[i] = (S(p .+ ϵ) - S(p .- ϵ)) / (2. * precision)
    end

    return p - pars.ϵ * normalize(conj(im * grad))
end

"""
    flow_points(S, thim::thimble, pars::parameters)

Flow the thimble one step in the downward flow.
"""
function flow_points(S, thim::thimble, pars::parameters)
    for p in thim.points
        if p.active && p.frozen == false
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

"""
    flow!(S, thim::thimble, pars::parameters)

Flow the thimble pars.N steps forward in the downward flow.
"""
function flow!(S, thim::thimble, pars::parameters)
    for _ in 1:pars.N
        flow_points(S, thim, pars)
    end
    remove_inactive_points(thim)
end

##############################
# Integrate
##############################

"""
    get_affine_map(vertices)

Affine map mapping a general simplex to a standard simplex.
"""
function get_affine_map(vertices)
    v₀ = vertices[1, :]
    B = transpose(vertices[2:end, :]) .- v₀
    return v₀, B
end

"""
    mapping(p, vertices)

Map the standard simplex to a point in the simplex spanned by vertices.
"""
function mapping(p, vertices)
    @assert length(p) == size(vertices, 2) "The dimension of p should be the same as the dimension of the points in vertices"
    v₀, B = get_affine_map(vertices)
    return v₀ + B * p
end

"""
    jacobian(vertices)

Evaluate the Jacobian of the affine map
"""
function jacobian(vertices)
    return det(get_affine_map(vertices)[2])
end

"""
    integrateSimplex(f, vertices, X, W)

Integrae f over a simplex spanned by vertices with a simplex quadrature method with the points X and the weights W.
"""
function integrateSimplex(f, vertices, X, W)
    integrand(p) = jacobian(vertices) * f(mapping(p, vertices))
    return sum(W[i] * integrand(X[i,:]) for i in 1:length(W))
end

"""
    PL_integrate(S, thim::thimble, pars::parameters)

Evaluate the integral ∫exp(im * S(x))dx along the thimble.
"""
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

# function triPlot(thim::thimble)
#     filter!(sim->sim.active, thim.simplices)
#     ps = [Polygon([Point2f(real.(thim.points[i].coord)) for i in sim.coord]) for sim in thim.simplices]
    
#     f = Figure()
#     ax = Axis(f[1, 1]; aspect=1)
#     poly!(ax, ps, color = rand(RGBf, length(ps)))

#     vertices = stack(real.(map(p -> p.coord, thim.points)), dims=1)
#     scatter!(ax, vertices[:, 1], vertices[:, 2])

#     for i in 1:length(thim.points)
#         text!(ax, vertices[i, 1], vertices[i, 2] + 0.2, text = string(i))
#     end

#     return f
# end
