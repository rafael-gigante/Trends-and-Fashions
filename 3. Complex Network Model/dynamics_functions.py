import numpy as np
import random
import networkx as nx

def generate_initial_conditions(N, L):
    """
    Generate the initial conditions of a system with `N` agents and `L` labels.

    Returns:
    - L_i: List with the label of each agent
    - N_i: List with the number of agents in each label at two times
    - p_i: List with the momentum of each label
    """
    L_i = [0] * N
    N_i = [[0] * L, [0] * L]
    p_i = [0] * L

    # Assigning a random label for each agent (We suppose that the previous label is also random)
    for i in range(N):
        N_i[0][random.randint(0, L-1)] += 1
        L_i[i] = random.randint(0, L-1)
        N_i[1][L_i[i]] += 1

    # Calculating the momentum of each label
    for i in range(L):
        p_i[i] = N_i[1][i] - N_i[0][i]

    N_i[0] = N_i[1].copy()

    return L_i, N_i, p_i


def interaction(i, j, N_i, L_i, p_i, p_crit, L):
    """
    Define the interaction between agents `i` and `j` in the system.
    """
    if p_i[L_i[j]] > p_i[L_i[i]]:
        N_i[1][L_i[i]] -= 1
        N_i[1][L_i[j]] += 1
        L_i[i] = L_i[j]
    elif p_i[L_i[i]] < p_crit:
        k = random.randint(0, L-1)
        if N_i[1][k] == 0:
            N_i[1][L_i[i]] -= 1
            N_i[1][k] += 1
            L_i[i] = k


def measures(N_i, p_i, N, L):
    """
    Calculate various measures of the system.

    Returns:
    - S: Entropy
    - Pc/N: Maximum proportion of agents in any label
    - Sc/N^2: Sum of squares of proportions of agents in labels
    """
    Pc = max(N_i[1])
    Sc = 0.0
    S = 0.0
    for i in range(L):
        p_i[i] = N_i[1][i] - N_i[0][i]
        if N_i[1][i] > 0:
            S += -(N_i[1][i] / N) * np.log(N_i[1][i] / N)
        Sc += N_i[1][i]**2
    Sc -= Pc**2
    N_i[0] = N_i[1].copy()
    return S, Pc / N, Sc / (N**2)


def generate_arrays(N, L, ensembles):
    """
    Set up the arrays that will store the data during the simulation.
    """
    iteration = []
    entropy = []
    percolation = []
    suscept = []
    tau = []

    for _ in range(len(N)):
        iteration.append([])
        entropy.append([])
        percolation.append([])
        suscept.append([])
        tau.append([])
        for _ in range(len(L)):
            iteration[-1].append([[] for _ in range(ensembles)])
            entropy[-1].append([[] for _ in range(ensembles)])
            percolation[-1].append([[] for _ in range(ensembles)])
            suscept[-1].append([[] for _ in range(ensembles)])
            tau[-1].append([])

    return iteration, entropy, percolation, suscept, tau

def generate_erdos_renyi(n, p):
    """
    Generate a Erdős-Rényi graph with n nodes and probability p of edge creation.
    """
    G = nx.erdos_renyi_graph(n, p)

    # Simplifying the network by considering only the largest component
    G.remove_edges_from(nx.selfloop_edges(G))
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    G = G.subgraph(Gcc[0])
    G = nx.convert_node_labels_to_integers(G, first_label=0)

    return G

def generate_watts_strogatz(n, k, p):
    """
    Generate a Watts-Strogatz graph with n nodes, k initial neighbours and probability p of rewire.
    """

    G = nx.watts_strogatz_graph(n, k, p)

    # Simplifying the network by considering only the largest component
    G.remove_edges_from(nx.selfloop_edges(G))
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    G = G.subgraph(Gcc[0])
    G = nx.convert_node_labels_to_integers(G, first_label=0)

    return G

def generate_barabasi_albert(n, m):
    """
    Generate a Barabasi-Albert network using preferential attachment.
    """
    G = nx.barabasi_albert_graph(n, m)

    # Simplifying the network by considering only the largest component
    G.remove_edges_from(nx.selfloop_edges(G))
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    G = G.subgraph(Gcc[0])
    G = nx.convert_node_labels_to_integers(G, first_label=0)

    return G