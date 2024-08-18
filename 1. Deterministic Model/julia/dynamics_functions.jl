module dynamics

"""
    generate_initial_conditions(N::Int, L::Int)

Generate the initial conditions of a system with `N` agents and `L` labels.

Returns:
- `L_i`: Array with the label of each agent
- `N_i`: Array with the number of agents in each label at two times
"""
function generate_initial_conditions(N::Int, L::Int)
    L_i = fill(0, N) 
    N_i = [fill(0, L), fill(0, L)]  
    p_i = fill(0, L)

    # Assigning a random label for each agent (We suppose that the previous label is also random)
    for i in 1:N
        N_i[1][rand(1:L)] += 1
        L_i[i] = rand(1:L)
        N_i[2][L_i[i]] += 1
    end

    # Calculating the momentum of each label
    for i in 1:L
        p_i[i] = N_i[2][i] - N_i[1][i]
    end

    N_i[1] .= N_i[2]

    return L_i, N_i, p_i
end

"""
    interaction(i::Int, j::Int, N_i::Vector{Vector{Int}}, L_i::Vector{Int}, p_i::Vector{Int}, p_crit::Int)

Define the interaction between agents `i` and `j` in the system.
"""
function interaction(i::Int, j::Int, N_i::Vector{Vector{Int}}, L_i::Vector{Int}, p_i::Vector{Int}, p_crit::Int, L::Int)
    if p_i[L_i[j]] > p_i[L_i[i]]
        N_i[2][L_i[i]] -= 1
        N_i[2][L_i[j]] += 1
        L_i[i] = L_i[j]
    elseif p_i[L_i[i]] < p_crit
        k = rand(1:L)
        if N_i[2][k] == 0
            N_i[2][L_i[i]] -= 1
            N_i[2][k] += 1
            L_i[i] = k
        end
    end
end

"""
    measures(N_i::Vector{Vector{Int}}, N::Int, L::Int)

Calculate various measures of the system.

Returns:
- `S`: Entropy
- `Pc/N`: Maximum proportion of agents in any label
- `Sc/N^2`: Sum of squares of proportions of agents in labels
"""
function measures(N_i::Vector{Vector{Int}}, p_i::Vector{Int}, N::Int, L::Int)
    Pc = maximum(N_i[2])
    Sc = 0.0
    S = 0.0
    for i in 1:L
        p_i[i] = N_i[2][i] - N_i[1][i]
        if N_i[2][i] > 0
            S += -(N_i[2][i] / N) * log(N_i[2][i] / N)
        end
        Sc += N_i[2][i]^2
    end
    Sc += -Pc^2
    N_i[1] .= N_i[2]
    return S, (Pc / N), (Sc / N^2)
end

"""
    generate_arrays(N, L, ensembles)

Set up the arrays that will store the data during the simulation.

"""
function generate_arrays(N, L, ensembles)
    iteration = []
    entropy = []
    percolation = []
    suscept = []
    τ = []
    
    for n in eachindex(N)
        push!(iteration, [])
        push!(entropy, [])
        push!(percolation, [])
        push!(suscept, [])
        push!(τ, [])
        for l in eachindex(L)
            push!(iteration[n], [])
            push!(entropy[n], [])
            push!(percolation[n], [])
            push!(suscept[n], [])
            push!(τ[n], [])
            for aux1 in 1:ensembles
                push!(iteration[n][l], [])
                push!(entropy[n][l], [])
                push!(percolation[n][l], [])
                push!(suscept[n][l], [])
            end
        end
    end

    return iteration, entropy, percolation, suscept, τ
end

end 
