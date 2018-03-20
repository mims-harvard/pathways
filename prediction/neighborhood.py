import numpy as np

""" Returns list of scores for the Neighborhood method, based on using assoc_gene_vector as the initial vector.
Scores are reported as 0 for the initial seed genes.

Inputs:
adjacency_matrix: Adjacency matrix of shape (number_genes, number_genes)
assoc_gene_vector: Numpy array of shape (number_genes,) with 1's for seed proteins, and 0's elsewhere
"""

def neighborhood_scores(adjacency_matrix, assoc_gene_vector):

	assoc_neighbors = np.dot(adjacency_matrix, assoc_gene_vector)
	total_neighbors = np.sum(adjacency_matrix, axis = 1)
	scores = assoc_neighbors/(total_neighbors+1e-50)
	scores = scores * (1 - assoc_gene_vector) 
	return scores
