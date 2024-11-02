import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import linregress

def plot_measures(iteration, entropy, percolation, suscept, N, L, ensembles, n, l):
    """
    Plot the entropy, percolation, and susceptibility given a value of 'n ∈ N' and 'l ∈ L' used in the simulation.
    """
    n_index = N.index(n)
    l_index = L.index(l)

    # Plot entropy
    plt.figure(figsize=(7, 3.5))
    plt.rc('font', size=10)        # Default font size
    plt.rc('axes', titlesize=12)   # Title font size
    plt.rc('axes', labelsize=10)   # X and Y label size
    plt.rc('xtick', labelsize=8)   # X tick label size
    plt.rc('ytick', labelsize=8)   # Y tick label size
    plt.rc('legend', fontsize=8)   # Legend font size
    plt.rc('lines', linewidth=2)   # Line width

    for i in range(ensembles):
        plt.plot(iteration[n_index][l_index][i], entropy[n_index][l_index][i], '--ko', markersize=5, label=None)
    plt.axhline(np.log(l), color='red', linestyle='--', label=r'$\ln(L)$')
    plt.xlabel("tempo [iterações]")
    plt.ylabel(r"$S$")
    plt.ylim(0, np.log(l) + 0.2)
    #plt.title(f"N={n}, L={l}")
    plt.legend(loc='right')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.savefig('1. Deterministic Model/python/plots/entropy.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Plot percolation
    plt.figure(figsize=(7, 3.5))
    for i in range(ensembles):
        plt.plot(iteration[n_index][l_index][i], percolation[n_index][l_index][i], '--ko', markersize=5, label=None)
    plt.xlabel("tempo [iterações]")
    plt.ylabel(r"$P_{\rm c}$")
    plt.ylim(0, 1.02)
    #plt.title(f"N={n}, L={l}")
    plt.grid(True, which='both', linestyle='--')
    plt.savefig('1. Deterministic Model/python/plots/percolation.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Plot susceptibility
    plt.figure(figsize=(7, 3.5))
    for i in range(ensembles):
        plt.plot(iteration[n_index][l_index][i], suscept[n_index][l_index][i], '--ko', markersize=5, label=None)
    plt.xlabel("tempo [iterações]")
    plt.ylabel(r"$S_{\rm c}$")
    #plt.title(f"N={n}, L={l}")
    plt.grid(True, which='both', linestyle='--')
    plt.savefig('1. Deterministic Model/python/plots/susceptibility.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_distribution(L, N_i, t):
    indices = np.arange(0, L, 1)
    plt.figure(figsize=(7, 3.5))
    plt.stem(indices, N_i[1], linefmt ='black', markerfmt=" ", basefmt=" ")
    plt.xlabel("índice da tendência")
    plt.ylabel("número de agentes")
    plt.xlim(0, L+1)
    plt.ylim(0, max(N_i[1])+ 0.1*max(N_i[1]))
    plt.grid(True, which='both', linestyle='--')
    plt.text(0.02, 0.95, f"tempo = {t}", transform=plt.gca().transAxes,
         fontsize=12, verticalalignment='top', bbox=dict(facecolor='white'))
    plt.savefig(f'1. Deterministic Model/python/plots/distribution_{t}.png', dpi=300, bbox_inches='tight')
    plt.close()


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
    plt.close()

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
    plt.close()
