# Implementation of the DIAMOnD method from Ghiassian et al

import numpy as np
from scipy.stats import hypergeom

# Avoids doing the costly hypergeometric calculation if the result will definitely be lower than an existing candidate
def _compare_to_existing(processed_list, seed_conns, total_conns):
    less_likely = [(a, b) for (a,b) in processed_list if a>=seed_conns and b<=total_conns]
    if len(less_likely) == 0: 
        return False
    return True

""" Returns list of scores, based on using assoc_gene_vector as the initial vector, for top number_to_rank predictions.

Inputs:
adjacency_matrix: Unnormalized adjacency matrix of shape (number_genes, number_genes)
assoc_gene_vector: Numpy array of shape (number_genes,) with 1's for seed proteins, and 0's elsewhere
alpha: The weight given to seeds
number_to_rank: The number of iterations to run DIAMOnD
"""
def diamond_scores(adjacency_matrix, assoc_gene_vector, alpha=5, number_to_rank=100):
    num_genes = assoc_gene_vector.shape[0]
    edges_per_gene = np.sum(adjacency_matrix, axis=0)
    scores = np.zeros(assoc_gene_vector.shape)
    seeds = np.copy(assoc_gene_vector)
    connections_to_seeds = np.sum(adjacency_matrix[:, np.nonzero(assoc_gene_vector)[0]], axis = 1)
    num_gene_edges = edges_per_gene + (alpha-1)*connections_to_seeds
    N = num_genes + np.sum(assoc_gene_vector)*(alpha-1)
    connections_to_seeds = connections_to_seeds * (alpha)
    num_seeds = alpha*np.sum(assoc_gene_vector)
    for index in range(1, number_to_rank+1):
        potential_cand = np.nonzero(connections_to_seeds*(1-seeds)>=1)[0]
        num_candidates = potential_cand.shape[0]
        if num_candidates == 0: break
        best_cand = -1
        best_conn = 1
        existing_calculations = []
        for i in range(num_candidates):
            cand_index = potential_cand[i]
            if _compare_to_existing(existing_calculations, connections_to_seeds[cand_index], num_gene_edges[cand_index]):
                continue
            conn = hypergeom.sf(connections_to_seeds[cand_index]-1, N, num_seeds, num_gene_edges[cand_index])
            existing_calculations.append((connections_to_seeds[cand_index], num_gene_edges[cand_index]))
            if conn < best_conn:
                best_conn = conn
                best_cand = cand_index
        connections_to_seeds += adjacency_matrix[:, best_cand]
        seeds[best_cand] = 1
        scores[best_cand] = 1.0/index
        num_seeds += 1
    return scores

