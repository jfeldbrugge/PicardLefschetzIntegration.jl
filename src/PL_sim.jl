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
    
    bounds = [[min[i], max[i]] for i in 1:pars.dim]
    for corner in Iterators.product(bounds...)
        push!(points, point(collect(corner)))
    end
    
    simplices = hypercube_triangulation_indices(pars.dim)

    divide_rep_sim(points, simplices, pars)
    return points, simplices
end

##############################
# Divide
##############################

function divide_sim(points, simplices, pars::parameters)
    for i in eachindex(simplices)
        sim = simplices[i]

        edges = filter(v -> issorted(v), collect(permutations(sim.coord, 2)))
        edge_lengths = [norm(points[edge[1]].coord - points[edge[2]].coord) for edge in edges]
    
        index_longest = argmax(edge_lengths)
        if edge_lengths[index_longest] > pars.δ
            simplices[i].active = false

            newPoint = point((points[edges[index_longest][1]].coord + points[edges[index_longest][2]].coord) ./ 2)
    
            push!(points, newPoint)
            index_newPoint = length(points)
                
            newsim1 =  setdiff(sim.coord, [edges[index_longest][1]]) 
            push!(newsim1, index_newPoint)
            
            newsim2 =  setdiff(sim.coord, [edges[index_longest][2]]) 
            push!(newsim2, index_newPoint)
            
            push!(simplices, index(newsim1))
            push!(simplices, index(newsim2))
        end
    end
    filter!(sim->sim.active, simplices)
end

function divide_rep_sim(points, simplices, pars::parameters)
    n_old = length(simplices)
    n_new = n_old + 1

    while n_old != n_new
        n_old = n_new
        divide_sim(points, simplices, pars)
        n_new = length(simplices)
    end
end

##############################
# Integrate
##############################
