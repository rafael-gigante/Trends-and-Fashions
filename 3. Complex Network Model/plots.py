import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import networkx as nx
from scipy.stats import linregress

def plot_measures(iteration, entropy, percolation, suscept, N, L, ensembles, n, l):
    """
    Plot the entropy, percolation, and susceptibility given a value of 'n ∈ N' and 'l ∈ L' used in the simulation.
    """
    n_index = N.index(n)
    l_index = L.index(l)

    # Plot entropy
    plt.figure()
    for i in range(ensembles):
        plt.plot(iteration[n_index][l_index][i], entropy[n_index][l_index][i], label=None, color='black')
    plt.axhline(np.log(l), color='red', linestyle='--', linewidth=1, label=r'$\log(L)$')
    plt.xlabel("time [iterations]")
    plt.ylabel(r"$S$")
    plt.ylim(0, np.log(l) + 0.2)
    plt.title(f"N={n}, L={l}")
    plt.legend(loc='right')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()

    # Plot percolation
    plt.figure()
    for i in range(ensembles):
        plt.plot(iteration[n_index][l_index][i], percolation[n_index][l_index][i], label=None, color='black')
    plt.axhline(1, color='black', linestyle='--', linewidth=1)
    plt.xlabel("time [iterations]")
    plt.ylabel(r"$P_c$")
    plt.ylim(0, 1.02)
    plt.title(f"N={n}, L={l}")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()

    # Plot susceptibility
    plt.figure()
    for i in range(ensembles):
        plt.plot(iteration[n_index][l_index][i], suscept[n_index][l_index][i], label=None, color='black')
    plt.xlabel("time [iterations]")
    plt.ylabel(r"$S_c$")
    plt.title(f"N={n}, L={l}")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()

def fit_data(x, y):
    """
    Receive two arrays and fit the data with a linear regression.

    Returns:
    - 'data': DataFrame with original and predicted data
    - 'slope': Float with the value of the angular coefficient
    - 'intercept': Float with the value of the linear coefficient
    """
    x = np.array(x)
    y = np.array(y)
    
    # Perform linear regression
    slope, intercept, _, _, _ = linregress(x, y)
    
    # Generate predicted values
    y_pred = intercept + slope * x

    # Create a DataFrame for easy handling
    data = pd.DataFrame({
        "X": x,
        "Y": y,
        "Y_pred": y_pred
    })

    return data, round(slope, 2), round(intercept, 2)

def plot_tau_fixed_L(tau, N, L):
    """
    Make a scatter plot of log(τ) vs log(N) for each value of l ∈ L and fit the data with a linear regression.
    """
    plt.figure()
    plt.xlabel(r"$\log(N)$")
    plt.ylabel(r"$\log(\tau)$")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    for l_index, l in enumerate(L):
        temp_tau = [tau[n_index][l_index] for n_index in range(len(N))]
        data, slope, intercept = fit_data(np.log10(N), np.log10(temp_tau))

        # Scatter plot
        plt.scatter(data["X"], data["Y"], label=f"L = {l}", color=plt.cm.tab10(l_index))

        # Plot fitted line
        plt.plot(data["X"], data["Y_pred"], linestyle='--', color=plt.cm.tab10(l_index),
                 label=f"$Y = {intercept} + {slope} * X$")

    plt.legend()
    plt.show()

def plot_tau_fixed_N(tau, N, L):
    """
    Make a scatter plot of log(τ) vs log(L) for each value of n ∈ N and fit the data with a linear regression.
    """
    plt.figure()
    plt.xlabel(r"$\log(L)$")
    plt.ylabel(r"$\log(\tau)$")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    for n_index, n in enumerate(N):
        temp_tau = [tau[n_index][l_index] for l_index in range(len(L))]
        data, slope, intercept = fit_data(np.log10(L), np.log10(temp_tau))

        # Scatter plot
        plt.scatter(data["X"], data["Y"], label=f"N = {n}", color=plt.cm.tab10(n_index))

        # Plot fitted line
        plt.plot(data["X"], data["Y_pred"], linestyle='--', color=plt.cm.tab10(n_index),
                 label=f"$Y = {intercept} + {slope} * X$")

    plt.legend()
    plt.show()

def add_values_and_colors(G, L_i, L):
    """
    Add the value of the label to the node and associate with a color.
    """
    for node in G.nodes():
        value = L_i[node]
        G.nodes[node]['value'] = value
        G.nodes[node]['color'] = plt.cm.gnuplot(value / L)  # Normalize the value to get a color

def plot_graph(G, L_i, L, pos):
    """
    Add the value of the label to the node and associate with a color.
    """

    add_values_and_colors(G, L_i, L)

    # Get the node values and corresponding colors
    node_values = nx.get_node_attributes(G, 'value')
    node_colors = [G.nodes[node]['color'] for node in G.nodes()]
    
    # Draw nodes with their values as labels
    plt.figure(figsize=(8, 6))
    nx.draw(G, pos, with_labels=True, labels=node_values, node_size=400, node_color=node_colors, font_size=10, font_weight='bold', edge_color='gray')

def create_frame(G, pos, save_dir, filenames, t, L_i, L):
    """
    Create a frame that shows the graph with the labels for the gif creation.
    """

    plot_graph(G, L_i, L, pos)
    # Save each frame
    filename = os.path.join(save_dir,f"frame_{t}.png")
    filenames.append(filename)
    plt.savefig(filename)
    plt.close()
    return filenames