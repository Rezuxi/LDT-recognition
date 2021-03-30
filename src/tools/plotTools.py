from .GraphTools import *
import networkx as nx
import asymmetree.treeevolve as te
from asymmetree.datastructures import PhyloTree
from asymmetree.hgt import ldt_graph
import matplotlib.pyplot as plt
import numpy as np

S = te.simulate_species_tree(10, model='innovation')
TGT = te.simulate_dated_gene_tree(S, dupl_rate = 0.5, loss_rate = 0.5, hgt_rate = 0.5, prohibit_extinction = "per_family", replace_prob=0.0)
OGT = te.observable_tree(TGT)
ldt = ldt_graph(OGT, S)

G = InvestigateGraph(ldt)


'''
	For the LDT graph and the perturbed graph just compare the amount of regions and P3 density of those regions

	for the perturbed graph and the edited graph, we can look in which regions (of the perturbed graph) the edits
	occurred.


	The size of a region (amount of nodes) can be represented by the size of a circle and the P3 density
	can be represented by the color of the circle. 


	Since cograph editing seems to fix the set of triples very often, check where the edits occurr
	in order to get some insight of why this is the case. For example, say we have many perturbed graphs and 
	apply cograph editing to them. Now in the cases we dont get an ldt graph, do the edits happen outside the regions?
	And in the cases we do get ldt graphs, do the edits happen inside? Are there any connections to be made between
	P3s and P4s? How likely are we to get an LDT-graph when the edits occurr inside/outside regions?

	Are dense P3 regions an indication of the set of triples being incompatible?
'''
def visualize_regions(G):
	pass



def plot_data(G):

	"""
		label should include variables for the TGT

		Bar plot for both cograph edit (blue) and triples edit (red) on perturbed graphs that are neither
		cographs nor have consistent set of triples.
	"""
	if G._count_dG_notCograph_notConsistent > 0:
		cograph_freq = G._count_cographEdit_fixed_consistency/G._count_dG_notCograph_notConsistent
	else:
		cograph_freq = 0
	if G._count_dG_notCograph_notConsistent > 0:
		triples_freq = G._count_triplesEdit_fixed_cograph/G._count_dG_notCograph_notConsistent
	else:
		triples_freq = 0
	

	height = [cograph_freq, triples_freq]
	bar_labels = ["cographEdit", "triplesEdit"]
	y_pos = np.arange(len(bar_labels))


	plt.bar(y_pos, height, color=['blue', 'red'])
	plt.xticks(y_pos, bar_labels)

	plt.show()

def plot_success_data(G):
	pass