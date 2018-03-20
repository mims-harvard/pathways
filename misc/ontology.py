from goatools.obo_parser import GODag
import numpy as np
from collections import defaultdict

"""
Ontology class that maps a list of diseases to a provided disease ontology_map.
Contains disease_class_map that maps a disease passed in to a map with the disease hierarchy.
Contains category_map that maps a disease class to diseases in it.

Input:
disease: List of disease names to be mapped
obo_file: The .obo file to be parsed that contains the disease ontology_map
"""
class Ontology(object):
	def __init__ (self, diseases, obo_file):
		self.diseases = diseases
		self.process_ontology(obo_file)

	def process_ontology(self, obo_file):
		obo_dag = GODag(obo_file)
		obo_dag.populate_terms()	
		diseases = self.diseases
		ontology_map = {}
		self.category_map = defaultdict(list)
		for item_id, item in obo_dag.items():
			# Considers those diseases whose names are subsets of the disease ontology names
			correlated_diseases = [disease for disease in diseases if disease in item.name]
			if len(correlated_diseases) > 0:
				d = {}
				for parent in (obo_dag.paths_to_top(item.id)[0]):
					d[parent.level] = parent.name
					for corr_disease in correlated_diseases:
						self.category_map[parent.name] += [corr_disease]
				for corr_disease in correlated_diseases:
					# Chooses the most specific disease mapping
					if corr_disease in ontology_map:
						if len(ontology_map[corr_disease]) > len(d):
							continue
					ontology_map[corr_disease] = d
		self.dag = obo_dag
		self.disease_class_map = ontology_map

