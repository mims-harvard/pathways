from itertools import product
import numpy as np
import networkx as nx
import math


"""
Returns community properties for disease modules

Inputs:
nx_graph: Networkx graph built from data
assoc_gene_vector: Numpy array of shape (number_genes,) with 1's for seed proteins, and 0's elsewhere
"""
def getCommunityScores(nx_graph, assoc_gene_vector):
    adjacency_matrixnx.adjacency_matrix(nx_graph)
    disease_properties = {}
    disease_indices = np.nonzero(assoc_gene_vector)[0]
    num_total_nodes = adjacency_matrix.shape[0]
    num_total_edges = np.sum(adjacency_matrix)/2

    num_nodes = len(disease_indices)
    disease_properties["num_nodes"] = num_nodes

    subgraph = nx_graph.subgraph(disease_indices)
    sliced_adj_matrix = adjacency_matrix[disease_indices, :]
    num_disease_node_edges = np.sum(sliced_adj_matrix)/2
    sub_adj_matrix = sliced_adj_matrix[:, disease_indices]
    num_internal_edges = np.sum(sub_adj_matrix)/2
    external_edges = num_disease_node_edges - num_internal_edges
    
    #Internal connectivity
    disease_properties["density"] = nx.density(subgraph)
    disease_properties["average_degree"] = 2*num_internal_edges/num_nodes
    disease_properties["average_internal_clustering"] = nx.average_clustering(subgraph)

    conn_comps = nx.connected_components(subgraph)
    sorted_cc = [len(c) for c in sorted(conn_comps, key=len, reverse=True)]
    disease_properties["size_largest_connected_component"] = sorted_cc[0]
    disease_properties["percent_in_largest_connected_component"] = float(sorted_cc[0])/num_nodes
    disease_properties["number_connected_components"] = len(sorted_cc)

    #External connectivity
    disease_properties["expansion"] = external_edges/num_nodes
    disease_properties["cut_ratio"] = external_edges/(num_nodes*(num_total_nodes-num_nodes))

    #External and internal connectivity
    disease_properties["conductance"] = external_edges/(2*num_internal_edges+external_edges)
    disease_properties["normalized_cut"] = disease_props["conductance"] + external_edges/(2*(num_total_edges-num_internal_edges)+external_edges)

    return disease_properties



