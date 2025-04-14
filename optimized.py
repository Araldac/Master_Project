import numpy as np
import random
from pathlib import Path
from scipy.integrate import solve_ivp
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.optimize import minimize_scalar

#efficient code is in efficient code
def system_dynamics_efficient(t, y, C, species_contrib, w, K, d, S, M, one_minus_l):
    N = y[:S]  # Species populations
    R = y[S:S+M]  # Resource concentrations

    # Initialize derivatives
    dN_dt = np.zeros(S)
    dR_dt = np.zeros(M)

    # Update dR_dt with the consumption and production
    dR_dt= (K-R)*d #resource entrance
    dR_dt-= (N@C)*R #conusmption
    production_contrib = species_contrib*N[:, None, None]*R[None, None, :] 
    dR_dt += np.array([np.sum(production_contrib[:, resourceA, :]) for resourceA in range(M)])#production
   
    
    # Precompute resource-specific contributions
    resource_contrib = w * R  # Shape (M,)

    # Compute growth rates and species dynamics in a vectorized way
    growth_rates = np.sum(one_minus_l * C * resource_contrib, axis=1)  # Shape (S,)
    dN_dt = N * (growth_rates - d)  # Compute dN/dt for all species at once

    # Apply constraints on population sizes
    dN_dt[N < 10e-6] = 0.0
    return np.concatenate((dN_dt, dR_dt))  # Return combined derivative of populations and resources

#from stoch_matrix import stoch_matrix
def stoch_matrix_pb(M, nestedness_pb, printplot=0):
    if nestedness_pb== -1:
        matrix = np.zeros((M, M))
    else: 
        matrix = np.full((M, M), 0.01)
        Rnumb=list(range(M))
        Rank = np.digitize(Rnumb, np.linspace(0, M))  # Assign ranks based on bins
        prob=0.5+nestedness_pb/2

        # Create the matrix based on ranks
        for j in Rnumb: #column
            for i in Rnumb: #rows 
                if Rank[j] < Rank[i]:
                    matrix[i, j] = np.random.choice([0, 1], size=1, p=[1-prob, prob])[0]
                else:
                    matrix[i, j] = np.random.choice([0, 1], size=1, p=[prob, 1-prob])[0]  # Random value for same or lower rank

        np.fill_diagonal(matrix, 0)  # Set diagonal elements to 0
        matrix=matrix*np.random.uniform(0,1, size=(M,M))
        col_sums = matrix.sum(axis=0, keepdims=True)
        col_sums[col_sums == 0] = 1  # Prevent division by zero
        matrix /= col_sums

        sparcity= 1-np.sum(matrix>0)/(M*M-M)
        if printplot == 1: 
            # Step 5: Visualize the matrix
            plt.figure(figsize=(12, 12))
            heat_map = sns.heatmap(matrix, linewidth=1, annot=False, cmap="YlGnBu")
            plt.gca().set_ylabel("Resource Number Produced")
            plt.gca().set_xlabel("Resource Number Consumed")
            plt.title(f"Stoichiometric Matrix with a nested probability of {nestedness_pb} and sparcity {sparcity}")
            plt.show()
    return matrix

#from count_coex import countcoex
def countcoex(S, N_results, solution,  num_points=50, threshold=0.05):
    count4 = 0
    surviving_sp=[]
    final_abundance=[]
    for species in range(S):
        if N_results[species, -1] > 10e-6:
            # Extract the last `num_points` from the data
            time_segment = solution.t[-num_points:]  # Last 50 time points
            data_segment = N_results[species, -num_points:]  # Last 50 corresponding data points

            # Apply linear regression: np.polyfit returns [slope, intercept]
            slope, intercept = np.polyfit(time_segment, np.log(data_segment + 1e-10), 1)
            if slope > -.01:
                count4+=1
                surviving_sp.append(species)
                final_abundance.append(N_results[species, -1])
    return count4, surviving_sp, final_abundance

def sample_with_target_equitability(M, target_eq, tolerance=0.05, n_avg=50):
    """Optimally find Dirichlet alpha to match target equitability (average over several samples)."""
    def entropy_diff(alpha):
        alphas = np.ones(M) * max(alpha, 0.01)
        entropies = []
        for _ in range(n_avg):
            proportions = np.random.dirichlet(alphas)
            H = -np.sum(proportions * np.log(proportions))
            entropies.append(H / np.log(M))
        avg_eq = np.mean(entropies)
        return abs(avg_eq - target_eq)
    
    res = minimize_scalar(entropy_diff, bounds=(0.01, 20), method='bounded')
    alpha_opt = max(res.x, 0.01)
    return np.random.dirichlet(np.ones(M) * alpha_opt)


def generate_consumption_matrix(S, M, specialization, printplot=0, Crowsum=1, sd=0.05):
    if specialization == -1:
        Et = None
    else:
        Et = 1 - specialization

    consumption_matrix = np.zeros((S, M))

    for species in range(S):
        # Target equitability
        if Et is None:
            equitability = np.round(np.random.random(), 2)
        else:
            equitability = np.clip(np.random.normal(Et, sd), 0, 1)

        proportions = sample_with_target_equitability(M, equitability)
        consumption_matrix[species, :] = proportions

    if printplot == 1:
        plt.figure(figsize=(12, 12))
        heat_map = sns.heatmap(consumption_matrix, linewidth=1, annot=False, cmap="YlGnBu")
        plt.gca().set_ylabel("Species Number")
        plt.gca().set_xlabel("Resource Number Consumed")
        plt.title(f"Consumption Matrix with a specialization of {specialization}")
        plt.show()

    return consumption_matrix


def run_simulation(S, M, specialization, nestedness_pb, t_max,  d, replica_number):
    seed1 = random.randrange(2**32 - 1)
    np.random.seed(seed1)
    K = np.array([320/Rgiven]*M)
    ngiv= np.random.choice(range(M), size=M-Rgiven, replace =False)
    K[ngiv] = 0

    # Generate consumer and stoichiometry matrices
    C = generate_consumption_matrix(S, M, specialization, printplot=0) #size (S,M)
    D = stoch_matrix_pb(M, nestedness_pb, printplot=0) #Size (M,M)

    # Initialize w and l to match target sums
    w = np.random.uniform(0, 1, M)
    w *= M / np.sum(w)
    l = np.random.uniform(0, 1, (S, M))
  # Precompute one minus l
    one_minus_l = 1 - l  # Shape (S, M)

    w_ratio = np.outer(1 / w, w)  # Shape (M, M)
    weighted_l = l[:, :, None] * w_ratio  # Shape (S, M, M)
    # Compute species-specific production contributions (constant for the simulation)
    species_c = np.zeros((S, M, M))  # Shape (S, M, M)
    for species in range(S):
        species_c[species] = D * C[species,None, :]  # Shape (M, M)
    #print(species_c.shape)    
    species_contrib=np.multiply(species_c,weighted_l)

    # Set initial conditions
    N0 = np.full(S, 0.01)
    R0 = np.full(M, 0)
    y0 = np.concatenate((N0, R0))
    t_span = (0, t_max)
    #time_points = np.linspace(0, 200, 200)

    # Run simulation
    solution = solve_ivp(lambda t, y: system_dynamics_efficient(t, y, C, species_contrib, w, K, d, S, M, one_minus_l), t_span, y0, "LSODA", t_eval=np.linspace(0, t_max, t_max*4))
    N_results = solution.y[:S, :]
    coex_sp, survivors, final_ab = countcoex(S, N_results, solution)
    
    new_C=C[survivors,:]

    #Shannon index end
    H_b=[]
    for i in range(S):
        H_b.append(1-(-sum([proportion*np.log(proportion) for proportion in  C[i,:] if proportion>0 ]))/np.log(M))
    #Shannon index beggining
    H_e=[]
    for i in range(len(survivors)):
        H_e.append(1-(-sum([proportion*np.log(proportion) for proportion in  new_C[i,:] if proportion>0]))/np.log(M))

    #Network analysis
    #b_module, b_conn= network(S,M, production=(C@D.T), consumption=C, nestedness=nestedness_pb, specialization=round(np.mean(H_b),2), printplot=1)
    #e_module, e_conn= network(len(survivors),M, production= (new_C@D.T),  consumption=new_C, nestedness=nestedness_pb, specialization=round(np.mean(H_e),2),printplot=1)
    #print(f'B/Ae module {b_module,e_module} B/A connectivity {b_conn, e_conn}')
    #plt.figure(figsize=(12, 6))
    #for species in range(S):
    #    plt.plot(solution.t, N_results[species, :], label=f'Species {species + 1}')
    #plt.xlabel('Time')
    #plt.ylabel('Population Size x10^8 (cells/ml)')
    #plt.title(f'Population Dynamics Over Time: {coex_sp} coexisting sp')
    #plt.legend(loc='upper right')
    #plt.grid()
    #plt.show()

    #print(f"seed{seed1} with {coex_sp} coexisting species")
    np.set_printoptions(legacy='1.25')
    return {
        "replica_number": replica_number,
        "seed": seed1,
        "S": S,
        "M": M,
        "specialization": specialization,
        "Spebefore": np.mean(H_b),
        "Speafter": np.mean(H_e),
        "nestedness": nestedness_pb,
        "t_max": t_max,
        "K": K.tolist(),
        "Resources given": np.sum(K>0),
        "dilution": d,
        "l_sum": np.sum(l),
        "nb_coexisting_species": coex_sp,
        "coexisting_species": survivors,
        "final_abundances": final_ab,
        #"D matrix":np.ravel(D),
        #"Consumer_matrix":np.ravel(C),
        #"Cbefore": b_conn,
        #"Cafter": e_conn,
        #"Mbefore": b_module,
        #"Mafter" : e_module

    }

if __name__ == "__main__":
    # Default simulation parameters
    M = 80
    t_max = 2500
    d=1
    replicas=10

    # Check if command-line arguments are provided
    if len(sys.argv) == 5:
        S= int(sys.argv[1])
        specialization = float(sys.argv[2])
        nestedness_pb = float(sys.argv[3])
        Rgiven= int(sys.argv[4])

       # ngiv= np.random.choice(range(M), size=M-Rgiven, replace =False)
       # K[ngiv] = 0
        #print(K)
        #K = np.array([10,10,10,10,0]*4)  # Adjusted K values
        # Run simulations for each replicate and.out store results
        results = []
        for replica_number in range(replicas): 
            #print(K)
            result = run_simulation(S, M, specialization, nestedness_pb, t_max, d, replica_number)
            results.append(result)


        # Save results to DataFrame and CSV
        df = pd.DataFrame(results)
        #print(df)
        csv_filename = f'/scratch/araldaco/masterproject/coex/results/try100/{S}species_{M}M_{Rgiven}Rgiven_{specialization}spe_{nestedness_pb}nested.csv'
        df.to_csv(csv_filename, index=False)

