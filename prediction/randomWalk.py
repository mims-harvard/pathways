import numpy as np

""" Returns list of scores, based on using assoc_gene_vector as the initial vector.
Scores are reported as 0 for the initial seed genes.

Inputs:
adjacency_matrix: Normalized adjacency matrix of shape (number_genes, number_genes)
assoc_gene_vector: Numpy array of shape (number_genes,) with 1's for seed proteins, and 0's elsewhere
return_prob: The probability of restart in the random walk
"""
def random_walk_scores(adjacency_matrix, assoc_gene_vector, return_prob=0.75):
	ratio = return_prob 
	convergence_metric = 1
	p0 = assoc_gene_vector/np.sum(assoc_gene_vector)
	old_vector = p0
	while (convergence_metric>1e-6):
		new_vector = (1-ratio) * np.dot(adjacency_matrix, old_vector) + ratio * p0
		convergence_metric = np.linalg.norm(new_vector-old_vector)
		old_vector = np.copy(new_vector)
	scores = old_vector * (1 - assoc_gene_vector) 
	return scores
