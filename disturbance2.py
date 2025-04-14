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
def system_dynamics_efficient(t, y, C, species_contrib, w, K, d, S, M, one_minus_l, e =0):
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
    dN_dt = N * (growth_rates - (d+e))  # Compute dN/dt for all species at once

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
def countcoex(S, N_results, solution, num_points=50, threshold=0.05):
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

def countcoex2(S, N_results, solution, num_points=50, threshold=0.05):
    count8 = 0
    surviving_sp8=[]
    final_abundance8=[]
    count16 = 0
    surviving_sp16=[]
    final_abundance16=[]
    for species in range(S):
        if N_results[species, 8] > 10e-6:
            # Extract the last `num_points` from the data
            time_segment = solution.t[0:8]  # Last 50 time points
            data_segment = N_results[species, 0:8]  # Last 50 corresponding data points

            # Apply linear regression: np.polyfit returns [slope, intercept]
            slope, intercept = np.polyfit(time_segment, np.log(data_segment + 1e-10), 1)
            if slope > -.01:
                count8+=1
                surviving_sp8.append(species)
                final_abundance8.append(N_results[species, 8])

        if N_results[species, 16] > 10e-6:
            # Extract the last `num_points` from the data
            time_segment = solution.t[8:16]  # Last 50 time points
            data_segment = N_results[species, 8:16]  # Last 50 corresponding data points

            # Apply linear regression: np.polyfit returns [slope, intercept]
            slope, intercept = np.polyfit(time_segment, np.log(data_segment + 1e-10), 1)
            if slope > -.01:
                count16+=1
                surviving_sp16.append(species)
                final_abundance16.append(N_results[species, 16])

    return count8, surviving_sp8, final_abundance8, count16, surviving_sp16, final_abundance16

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

#from tryshannon import generate_consumption_matrix
#def generate_consumption_matrix(S, M, specialization, printplot=0, Crowsum=1, sd=0.05):
#    if specialization == -1 :
#        consumption_matrix = np.zeros((S, M))
#        for species in range(S):
            # Step 1: Draw equitability for the species from a normal distribution centered on Et
#            equitability = np.round(np.random.random(),2)

            # Step 2: Generate random proportions (pi) that yield the target equitability E
            # Initialize proportions
#            proportions = np.random.dirichlet(np.ones(M))
            
            # Adjust proportions iteratively to match the equitability
#            achieved_equitability = -np.sum(proportions * np.log(proportions)) / np.log(M)
#            if M> 15: 
#                tolerance =0.1
#            else:
#                tolerance = 0.05  # Set tolerance for equitability matching

            # Use iterative adjustment to match equitability
#            while abs(achieved_equitability - equitability) > tolerance:
                # Adjust proportions based on whether we need higher or lower equitability
#                proportions = np.random.dirichlet(np.ones(M) * (1 + (equitability - achieved_equitability)))
#                achieved_equitability = -np.sum(proportions * np.log(proportions)) / np.log(M)
            
            # Assign the adjusted proportions to the consumption matrix
#            consumption_matrix[species, :] = proportions
#    else:
#        Et=1-specialization
#        consumption_matrix = np.zeros((S, M))

#        for species in range(S):
            # Step 1: Draw equitability for the species from a normal distribution centered on Et
#            equitability = np.clip(np.random.normal(Et, sd), 0, 1)

            # Step 2: Generate random proportions (pi) that yield the target equitability E
            # Initialize proportions
#            proportions = np.random.dirichlet(np.ones(M))
            
            # Adjust proportions iteratively to match the equitability
#            achieved_equitability = -np.sum(proportions * np.log(proportions)) / np.log(M)
            #tolerance = 0.01  # Set tolerance for equitability matching
#            if M> 15: 
#                tolerance =0.1
#            else:
#                tolerance = 0.05  # Set tolerance for equitability matching

            # Use iterative adjustment to match equitability
#            while abs(achieved_equitability - equitability) > tolerance:
                # Adjust proportions based on whether we need higher or lower equitability
#                proportions = np.random.dirichlet(np.ones(M) * (1 + (equitability - achieved_equitability)))
#                achieved_equitability = -np.sum(proportions * np.log(proportions)) / np.log(M)
            
#            # Assign the adjusted proportions to the consumption matrix
#            consumption_matrix[species, :] = proportions
#    if printplot ==1: 
#        # Step 5: Visualize the matrix
#        plt.figure(figsize=(12, 12))
#        heat_map = sns.heatmap(consumption_matrix, linewidth=1, annot=False, cmap="YlGnBu")
#        plt.gca().set_ylabel("Species Number")
#        plt.gca().set_xlabel("Resource Number Consumed")
#        plt.title(f"Consumption Matrix with a specialization of {specialization}")
#        plt.show()

#    return consumption_matrix


def run_simulation(S, M, specialization, nestedness_pb, t_max, seed1, d, replica_number, Rgiven):
    np.random.seed(seed1)
    K = np.array([320/Rgiven]*M)
    ngiv= np.random.choice(range(M), size=M-Rgiven, replace =False)
    K[ngiv] = 0
    print(K)
    

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
    solution = solve_ivp(lambda t, y: system_dynamics_efficient(t, y, C, species_contrib, w, K, d, S, M, one_minus_l), t_span, y0, "LSODA",t_eval=np.linspace(0, t_max, t_max))
    N_results = solution.y[:S, :]
    coex_sp, survivors, final_ab = countcoex(S, N_results, solution)
    
    new_C=C[survivors,:]

    #Shannon index end
    H_b=[]
    for i in range(S):
        H_b.append(1-(-sum([proportion*np.log(proportion) for proportion in  C[i,:] if proportion>0 ]))/np.log(M))
    #Shannon index beggining
    H_e=[]
    #for i in range(len(survivors)):
    #    H_e.append(1-(-sum([proportion*np.log(proportion) for proportion in  new_C[i,:] if proportion>0]))/np.log(M))

    plt.figure(figsize=(12, 6))
    for species in range(S):
        plt.plot(range(0, 2500), N_results[species, :], label=f'Species {species + 1}')
    plt.xlabel('Time')
    plt.ylabel('Population Size x10^8 (cells/ml)')
    plt.title(f'Population Dynamics after disturbance from {coex_sp} to coexsiting species')
    plt.legend(loc='upper right')
    plt.grid()
    plt.axvline(x=2500, color= "r", linestyle= "dashed" )
    plt.show()

    #np.set_printoptions(legacy='1.25')
    result = {
        "replica_number": None,
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
        "e": 0,
        "l_sum": np.sum(l),
        "nb_coexisting_species": coex_sp,
        "coexisting_species": survivors,
        "final_abundances": final_ab,
        "nb_coexisting_species_end": None,
        "coexisting_species_end":None ,
        "final_abundances_end":None ,
        "nb_coexisting_species8": None,
        "coexisting_species8": None,
        "final_abundances8":None ,
        "nb_coexisting_species16": None,
        "coexisting_species16": None,
        "final_abundances16": None,
        "type": "coex",
        "Value": None,
        "RS": None}

    results.append(result)
    
    ###Disturbance part
    
    y0 = np.array(solution.y[:, -1])
    t_max=2500
    t_span = (0, t_max) 

    for type_d in ["dil", "Ant", "R", "S"]: 
        if type_d == "dil": 

            for dil in range(1, 40, 2):
                d2=dil
                solution2 = solve_ivp(lambda t, y: system_dynamics_efficient(t, y, C, species_contrib, w, K, d2, S, M, one_minus_l), t_span, y0, "LSODA", t_eval=np.linspace(0, t_max, t_max))
                N_results2 = solution2.y[:S, :]
                coex_sp2, survivors2, final_ab2 = countcoex(S, N_results2, solution2)
                coex_sp8, survivors8, final_ab8, coex_sp16, survivors16, final_ab16  = countcoex2(S, N_results2, solution2)

                totalNResults=  np.hstack((N_results, N_results2))
        

                plt.figure(figsize=(12, 6))
                for species in range(S):
                    plt.plot(range(0, 5000), totalNResults[species, :], label=f'Species {species + 1}')
                plt.xlabel('Time')
                plt.ylabel('Population Size x10^8 (cells/ml)')
                plt.title(f'Population Dynamics after disturbance {type_d} from {coex_sp} to {coex_sp2} coexsiting species')
                plt.legend(loc='upper right')
                plt.grid()
                plt.axvline(x=2500, color= "r", linestyle= "dashed" )
                plt.show()

                result= {
                "replica_number": None,
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
                "dilution": d2,
                "e": 0,
                "l_sum": np.sum(l),
                "nb_coexisting_species": coex_sp,
                "coexisting_species": survivors,
                "final_abundances": final_ab,
                "nb_coexisting_species_end": coex_sp2,
                "coexisting_species_end": survivors2,
                "final_abundances_end": final_ab2,
                "nb_coexisting_species8": coex_sp8,
                "coexisting_species8": survivors8,
                "final_abundances8": final_ab8,
                "nb_coexisting_species16": coex_sp16,
                "coexisting_species16": survivors16,
                "final_abundances16": final_ab16,
                "type": type_d,
                "Value":d2,
                "RS": None}
            
                results.append(result)

        if type_d == "R":
            #Make a *copy* of K to work with
            K2 = K.copy()
            positive_indices = np.where(K2 > 0)[0]
            max_removal = len(positive_indices)

            for RR in range(0, max_removal + 1):
                if RR <= max_removal:
                    to_remove = random.sample(list(positive_indices), k=RR)
                    K2[to_remove] = 0

                solution2 = solve_ivp(lambda t, y: system_dynamics_efficient(t, y, C, species_contrib, w, K2, d, S, M, one_minus_l), t_span, y0, "LSODA",t_eval=np.linspace(0, t_max, t_max))
                N_results2 = solution2.y[:S, :]
                coex_sp2, survivors2, final_ab2 = countcoex(S, N_results2, solution2)
                coex_sp8, survivors8, final_ab8, coex_sp16, survivors16, final_ab16 = countcoex2(S, N_results2, solution2)
                totalNResults=  np.hstack((N_results, N_results2))
          
                   

                plt.figure(figsize=(12, 6))
                for species in range(S):
                    plt.plot(range(0, 5000), totalNResults[species, :], label=f'Species {species + 1}')
                plt.xlabel('Time')
                plt.ylabel('Population Size x10^8 (cells/ml)')
                plt.title(f'Population Dynamics after disturbance {type_d} from {coex_sp} to {coex_sp2} coexsiting species')
                plt.legend(loc='upper right')
                plt.grid()
                plt.axvline(x=2500, color= "r", linestyle= "dashed" )
                plt.show()

                result= {
                "replica_number": None,
                "seed": seed1,
                "S": S,
                "M": M,
                "specialization": specialization,
                "Spebefore": np.mean(H_b),
                "Speafter": np.mean(H_e),
                "nestedness": nestedness_pb,
                "t_max": t_max,
                "K": K2.tolist(),
                "Resources given": np.sum(K2>0),
                "dilution": d,
                "e": 0,
                "l_sum": np.sum(l),
                "nb_coexisting_species": coex_sp,
                "coexisting_species": survivors,
                "final_abundances": final_ab,
                "nb_coexisting_species_end": coex_sp2,
                "coexisting_species_end": survivors2,
                "final_abundances_end": final_ab2,
                "nb_coexisting_species8": coex_sp8,
                "coexisting_species8": survivors8,
                "final_abundances8": final_ab8,
                "nb_coexisting_species16": coex_sp16,
                "coexisting_species16": survivors16,
                "final_abundances16": final_ab16,
                "type": type_d,
                "Value":RR/max_removal,
                "RS":None}

                results.append(result)

        
        if type_d == "S":
            for RS in range(0, 22):
                indices_y0 = np.where(y0[:S] > 10e-6 )[0]
                y02=y0
                if len(indices_y0) >= RS:  # Randomly select one index and set it to 0
                    random_index = random.sample(list(indices_y0), k=RS)
                    y02[random_index] = 0

                solution2 = solve_ivp(lambda t, y: system_dynamics_efficient(t, y, C, species_contrib, w, K, d, S, M, one_minus_l), t_span, y02, "LSODA", t_eval=np.linspace(0, t_max, t_max))
                N_results2 = solution2.y[:S, :]
                coex_sp2, survivors2, final_ab2 = countcoex(S, N_results2, solution2)
                coex_sp8, survivors8, final_ab8, coex_sp16, survivors16, final_ab16 = countcoex2(S, N_results2, solution2)
                totalNResults=  np.hstack((N_results, N_results2))
        

                plt.figure(figsize=(12, 6))
                for species in range(S):
                    plt.plot(range(0, 5000), totalNResults[species, :], label=f'Species {species + 1}')
                plt.xlabel('Time')
                plt.ylabel('Population Size x10^8 (cells/ml)')
                plt.title(f'Population Dynamics after disturbance {type_d} from {coex_sp} to {coex_sp2} coexsiting species')
                plt.legend(loc='upper right')
                plt.grid()
                plt.axvline(x=2500, color= "r", linestyle= "dashed" )
                plt.show()

                result= {
                "replica_number": None,
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
                "e": 0,
                "l_sum": np.sum(l),
                "nb_coexisting_species": coex_sp,
                "coexisting_species": survivors,
                "final_abundances": final_ab,
                "nb_coexisting_species_end": coex_sp2,
                "coexisting_species_end": survivors2,
                "final_abundances_end": final_ab2,
                "nb_coexisting_species8": coex_sp8,
                "coexisting_species8": survivors8,
                "final_abundances8": final_ab8,
                "nb_coexisting_species16": coex_sp16,
                "coexisting_species16": survivors16,
                "final_abundances16": final_ab16,
                "type": type_d,
                "Value": RS,
                "RS":RS}

                results.append(result)

        if type_d == "Ant":
            for A in range(0, 30, 2):
                #e=np.random.uniform(1, A, size=S)
                e=np.random.normal(A, A/4, size=S)
 
                solution2 = solve_ivp(lambda t, y: system_dynamics_efficient(t, y, C, species_contrib, w, K, d, S, M, one_minus_l, e), t_span, y0, "LSODA", t_eval=np.linspace(0, t_max, t_max))
                N_results2 = solution2.y[:S, :]
                coex_sp2, survivors2, final_ab2 = countcoex(S, N_results2, solution2)
                coex_sp8, survivors8, final_ab8, coex_sp16, survivors16, final_ab16  = countcoex2(S, N_results2, solution2)
                totalNResults=  np.hstack((N_results, N_results2))
        

                plt.figure(figsize=(12, 6))
                for species in range(S):
                    plt.plot(range(0, 5000), totalNResults[species, :], label=f'Species {species + 1}')
                plt.xlabel('Time')
                plt.ylabel('Population Size x10^8 (cells/ml)')
                plt.title(f'Population Dynamics after disturbance {A} {e} from {coex_sp} to {coex_sp2} coexsiting species')
                plt.legend(loc='upper right')
                plt.grid()
                plt.axvline(x=2500, color= "r", linestyle= "dashed" )
                plt.show()

                result= {
                "replica_number": None,
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
                "e": A, 
                "l_sum": np.sum(l),
                "nb_coexisting_species": coex_sp,
                "coexisting_species": survivors,
                "final_abundances": final_ab,
                "nb_coexisting_species_end": coex_sp2,
                "coexisting_species_end": survivors2,
                "final_abundances_end": final_ab2,
                "nb_coexisting_species8": coex_sp8,
                "coexisting_species8": survivors8,
                "final_abundances8": final_ab8,
                "nb_coexisting_species16": coex_sp16,
                "coexisting_species16": survivors16,
                "final_abundances16": final_ab16,
                "type": type_d,
                "Value": A,
                "RS": None}

                results.append(result)
        #totalNResults=  np.hstack((N_results, N_results2))
    

        #plt.figure(figsize=(12, 6))
        #for species in range(S):
        #    plt.plot(range(0, 5000), totalNResults[species, :], label=f'Species {species + 1}')
        #plt.xlabel('Time')
        #plt.ylabel('Population Size x10^8 (cells/ml)')
        #plt.title(f'Population Dynamics after disturbance {disturbance_type} from {coex_sp} to {coex_sp2} coexsiting species')
        #plt.legend(loc='upper right')
        #plt.grid()
        #plt.axvline(x=2500)
        #plt.show()  



if __name__ == "__main__":
    # Default simulation parameters
    M = 80
    t_max = 2500
    d=1
    replica_number=1

    # Check if command-line arguments are provided
    if len(sys.argv) == 6:
        S= int(sys.argv[1])
        specialization = float(sys.argv[2])
        nestedness_pb = float(sys.argv[3])
        Rgiven= int(sys.argv[4])
        seed1= int(sys.argv[5])
        #disturbance_type= sys.argv[6]
        #K = np.array([int(x.strip()) for x in sys.argv[6].split(',')])
        #print(K)
        results = []

        run_simulation(S, M, specialization, nestedness_pb, t_max, seed1, d, replica_number, Rgiven)
   
        # Save results to DataFrame and CSV
        df = pd.DataFrame(results)
        #print(df)
        #csv_filename = f'/scratch/araldaco/masterproject/coex/results/disturbance2/{seed1}seed{Rgiven}Rgiven_{specialization}spe_{nestedness_pb}nested.csv'
        #csv_filename = f'trydist.csv'
        
        df.to_csv(csv_filename, index=False)
