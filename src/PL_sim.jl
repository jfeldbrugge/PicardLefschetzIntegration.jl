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

function initialGrid_sim(min, max, pars::parameters)
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
    divide_rep_sim(points, simplices, pars)
    return points, simplices
end

function volume_simplex(vertices)
    return det(transpose(vertices[2:end, :]) .- vertices[1, :]) / factorial(size(vertices, 1) - 1)
end

##############################
# Divide
##############################

function divide_sim(points, simplices, pars::parameters)
    dictionary = Dict()
    for i in eachindex(simplices)
        sim_coord = simplices[i].coord
        if simplices[i].active
            edges = filter(v -> issorted(v), collect(permutations(sim_coord, 2)))
    
            # Check whether an edge of the simplex sim if it has already been subdivided
            subdivided = false
            for edge in edges
                if haskey(dictionary, edge) && subdivided == false
                    subdivided = true
                    simplices[i].active = false
    
                    index_newPoint = dictionary[edge]
    
                    newsim1 = copy(sim_coord)
                    newsim2 = copy(sim_coord)
                    replace!(newsim1, edge[1] => index_newPoint)
                    replace!(newsim2, edge[2] => index_newPoint)
                    
                    push!(simplices, index(newsim1))
                    push!(simplices, index(newsim2))
                end
            end
    
            # Check the length of the edges and subdivide the simpex when the longest edge exceeds pars.δ
            if subdivided == false 
                edge_lengths = [norm(points[edge[1]].coord - points[edge[2]].coord) for edge in edges]
            
                index_longest = argmax(edge_lengths)
                if edge_lengths[index_longest] > pars.δ
                    simplices[i].active = false
    
                    newPoint = point((points[edges[index_longest][1]].coord + points[edges[index_longest][2]].coord) ./ 2)
            
                    push!(points, newPoint)
                    index_newPoint = length(points)
    
                    newsim1 = copy(sim_coord)
                    newsim2 = copy(sim_coord)
                    replace!(newsim1, edges[index_longest][1] => index_newPoint)
                    replace!(newsim2, edges[index_longest][2] => index_newPoint)
                    
                    push!(simplices, index(newsim1))
                    push!(simplices, index(newsim2))
    
                    dictionary[edges[index_longest]] = index_newPoint
                end
            end
        end
    end
end

function divide_rep_sim(points, simplices, pars::parameters)
    n_old = length(simplices)
    n_new = n_old + 1

    while n_old != n_new
        n_old = n_new
        divide_sim(points, simplices, pars)
        n_new = length(simplices)
    end
    filter!(sim->sim.active, simplices)
end

function remove_inactive_points(points, simplices)
    dictionary = Dict()
    counter = 0
    for i in eachindex(points)
        if points[i].active
            counter += 1
            dictionary[i] = counter
        end
    end

    for i in eachindex(simplices)
        coord = copy(simplices[i].coord)
        simplices[i].coord = [dictionary[j] for j in coord]
    end

    filter!(p -> p.active, points)

    return points, simplices
end

##############################
# Flow
##############################

function flow_point_sim(S, p, pars::parameters)
    precision = 1e-10

    grad = zeros(ComplexF64, pars.dim)
    
    for i in 1:pars.dim
        ϵ = zeros(pars.dim)
        ϵ[i] = precision

        grad[i] = (S(p .+ ϵ) - S(p .- ϵ)) / (2. * precision)
    end

    return p - pars.ϵ * normalize(conj(im * grad))
end

function flow_points_sim(S, points, simplices, pars::parameters)
    for p in points
        if p.active
            p.coord = flow_point_sim(S, p.coord, pars)
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

    divide_rep_sim(points, simplices, pars)
end

function flow_sim(S, points, simplices, pars::parameters)
    for _ in 1:pars.N
        flow_points_sim(S, points, simplices, pars)
    end
    remove_inactive_points(points, simplices)
end

##############################
# Integrate
##############################

function get_affine_map(vertices)
    v₀ = vertices[1, :]
    B = transpose(vertices[2:end, :]) .- v₀
    return v₀, B
end

function mapping_sim(p, vertices)
    v₀, B = get_affine_map(vertices)
    return v₀ + B * p
end

function jacobian_sim(vertices)
    return det(get_affine_map(vertices)[2])
end

function integrateSimplex(f, vertices, X, W)
    integrand(p) = jacobian_sim(vertices) * f(mapping_sim(p, vertices))
    return sum(W[i] * integrand(X[i,:]) for i in 1:length(W))
end

function PL_sim(S, thim::thimble, pars::parameters)
    X, W = simplexquad(pars.n, pars.dim)

    points_r = map(p->p.coord, thim.points)
    simplices_r = map(sim->sim.coord, thim.simplices)

    sum = 0
    for sim in simplices_r
        sum += integrateSimplex(p -> exp(im * S(p)), stack(points_r[sim], dims=1), X, W)
    end
    return sum
end