import random
import imageio.v2 as imageio
import os
import networkx as nx
from dynamics_functions import generate_initial_conditions, interaction, measures, generate_arrays, generate_erdos_renyi
from plots import plot_measures, create_frame, plot_tau_fixed_L, plot_tau_fixed_N

# Initial parameters
N = [100, 316, 1000, 3162, 10000] # Number of agents
L = [10, 100, 1000, 10000] # Number of labels
p_crit = 1 # Value of p_crit
p = 0.1  # Probability for edge creation (Érdos-Renyi)
ensembles = 20 # Number of ensembles for each simulation

# Arrays to store the measures
iteration, entropy, percolation, suscept, tau = generate_arrays(N, L, ensembles)

# Get the directory where the script is located
#script_dir = os.path.dirname(os.path.abspath(__file__))
# Define a relative path for saving frames and GIF
#relative_path = 'animation_frames'
#save_dir = os.path.join(script_dir, relative_path)
#os.makedirs(save_dir, exist_ok=True)  # Create the directory if it doesn't exist
#filenames = []

# Loop through the values of N
for n_index, n in enumerate(N):
    # Loop through the values of L
    for l_index, l in enumerate(L):
        fixation_time = 0
        for _ in range(ensembles):
            t = 0
            L_i, N_i, p_i = generate_initial_conditions(n, l)
            G = generate_erdos_renyi(n, p)
            #pos = nx.spring_layout(G)  # positions for all nodes

            #filenames = create_frame(G, pos, save_dir, filenames, t, L_i, l)

            # The dynamics will occur until one label is fully occupied
            # while (max(N_i[1])/n) < 1:
            for aux in range(100):
                # Perform interaction for two random agents N times
                for it in range(n):
                    i = random.randint(0, n - 1)
                    neighbors = list(G.neighbors(i))
                    j = random.choice(neighbors)
                    interaction(i, j, N_i, L_i, p_i, p_crit, l)

                # Compute relevant measures
                S, Pc, Sc = measures(N_i, p_i, n, l)
                iteration[n_index][l_index][_].append(t)
                entropy[n_index][l_index][_].append(S)
                percolation[n_index][l_index][_].append(Pc)
                suscept[n_index][l_index][_].append(Sc)

                t += 1
                #filenames = create_frame(G, pos, save_dir, filenames, t, L_i, l)

            fixation_time += t
        tau[n_index][l_index] = fixation_time / ensembles

# Create a GIF
#gif_path = os.path.join(save_dir, 'ER_evolution.gif')
#with imageio.get_writer(gif_path, mode='I', duration=1) as writer:
#    for filename in filenames:
#        image = imageio.imread(filename)
#        writer.append_data(image)

# Clean up temporary image files
#for filename in filenames:
#    os.remove(filename)

#plot_measures(iteration, entropy, percolation, suscept, N, L, ensembles, N[0], L[0])
plot_tau_fixed_L(tau, N, L)
plot_tau_fixed_N(tau, N, L)