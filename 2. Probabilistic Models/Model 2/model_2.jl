include("dynamics_functions.jl")
include("plots.jl")

# Parameters of the simulation
# N := number of agents, L := number of labels, β := temperature
# ensembles := # of ensembles for each configuration, n_greatest := # of the greatest labels to track in each configuration and ensemble
N = [10, 31, 10^2, 316, 10^3]
L = [10, 10^2, 10^3]
ensembles = 20

# Arrays with the measures
iteration, entropy, percolation, suscept, τ = dynamics.generate_arrays(N, L, ensembles)

# Passes through the values of N
for n in eachindex(N)
    # Passes through the values of L
    for l in eachindex(L)
        fixation_time = 0
        for aux1 in 1:ensembles
            L_i, N_i, p_i = dynamics.generate_initial_conditions(N[n], L[l])

            t = 0
            # The dynamics will occur until one label is fully occupied
            while maximum(N_i[2] / N[n]) != 1
                # Makes the interaction of two random agents N times
                for aux2 in 1:N[n]
                    i = rand(1:N[n])
                    j = rand(1:N[n])
                    dynamics.interaction(i, j, N_i, L_i, p_i, N[n])
                end

                # Makes some relevant measures
                S, Pc, Sc = dynamics.measures(N_i, p_i, N[n], L[l])
                push!(iteration[n][l][aux1], t)
                push!(entropy[n][l][aux1], S)
                push!(percolation[n][l][aux1], Pc)
                push!(suscept[n][l][aux1], Sc)

                t += 1
            end
            fixation_time += t
        end
        τ[n][l] = fixation_time / ensembles
    end
end

plots.plot_measures(iteration, entropy, percolation, suscept, N, L, ensembles, 10^3, 10^2)
plots.plot_τ_fixed_L(τ, N, L)
plots.plot_τ_fixed_N(τ, N, L)