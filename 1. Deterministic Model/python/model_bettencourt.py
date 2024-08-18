import numpy as np
import random
from dynamics_functions import generate_initial_conditions, interaction, measures, generate_arrays
from plots import plot_measures, plot_tau_fixed_L, plot_tau_fixed_N

# Parameters of the simulation
# N := number of agents, L := number of labels, p_crit := critical momentum threshold
# ensembles := number of ensembles for each configuration
N = [10**3, 10**4, 10**5]
L = [10, 10**2, 10**3]
p_crit = 1
ensembles = 1

# Arrays to store the measures
iteration, entropy, percolation, suscept, tau = generate_arrays(N, L, ensembles)

# Loop through the values of N
for n_index, n in enumerate(N):
    # Loop through the values of L
    for l_index, l in enumerate(L):
        fixation_time = 0
        for _ in range(ensembles):
            L_i, N_i, p_i = generate_initial_conditions(n, l)

            t = 0
            # The dynamics will occur until one label is fully occupied
            while (max(N_i[1])/n) < 1:
                # Perform interaction for two random agents N times
                for it in range(n):
                    i = random.randint(0, n - 1)
                    j = random.randint(0, n - 1)
                    interaction(i, j, N_i, L_i, p_i, p_crit, l)

                # Compute relevant measures
                S, Pc, Sc = measures(N_i, p_i, n, l)
                iteration[n_index][l_index][_].append(t)
                entropy[n_index][l_index][_].append(S)
                percolation[n_index][l_index][_].append(Pc)
                suscept[n_index][l_index][_].append(Sc)

                t += 1

            fixation_time += t
        tau[n_index][l_index] = fixation_time / ensembles

# Plotting the results
plot_measures(iteration, entropy, percolation, suscept, N, L, ensembles, 10**5, 10**3)
plot_tau_fixed_L(tau, N, L)
plot_tau_fixed_N(tau, N, L)
