include("dynamics_functions.jl")
include("plots.jl")

# Parameters of the simulation
# N := number of agents, L := number of labels, β := temperature
# ensembles := # of ensembles for each configuration, n_greatest := # of the greatest labels to track in each configuration and ensemble
N = [10, 31, 10^2, 316, 10^3]
L = [10^2]
β = [1, 1e1, 25, 5e1, 1e2]
ensembles = 20
n_greatest = 5

# Arrays with the measures
iteration, entropy, percolation, suscept, τ, labels_evolution = dynamics.generate_arrays(N, L, β, ensembles, n_greatest)

# Passes through the values of N
for n in eachindex(N)
    # Passes through the values of L
    for l in eachindex(L)
        # Passes through the values of β
        for b in eachindex(β)
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
                        dynamics.interaction(i, j, N_i, L_i, p_i, N[n], β[b])
                    end

                    # Makes some relevant measures
                    S, Pc, Sc = dynamics.measures(N_i, p_i, N[n], L[l])
                    push!(iteration[n][l][b][aux1], t)
                    push!(entropy[n][l][b][aux1], S)
                    push!(percolation[n][l][b][aux1], Pc)
                    push!(suscept[n][l][b][aux1], Sc)

                    greatest_values = dynamics.find_n_greatest(N_i[2]/N[n], n_greatest)
                    for i in 1:n_greatest
                        push!(labels_evolution[n][l][b][aux1][i], greatest_values[i])
                    end

                    t += 1
                end
                fixation_time += t
            end
            τ[n][l][b] = fixation_time / ensembles
        end
    end
end

plots.plot_fermi_function(β)
plots.plot_measures(iteration, entropy, percolation, suscept, N, L, β, ensembles, 10^3, 10^2, 1e2)
plots.plot_labels_evolution(iteration, labels_evolution, N, L, β, ensembles, n_greatest, 10^3, 10^2, 1e2)
plots.plot_τ_fixed_l(τ, N, L, β, 10^2)
plots.plot_τ_over_N(τ, N, L, β, 10^2)