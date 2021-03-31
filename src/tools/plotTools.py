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
		label should include variables for the TGT(?)

	"""

	# Plots the frequency of cograph editing fixing set of triples and triples editing fixing cograph
	# frequency of each editing heuristic turns the graph into an LDT-graph (properly colored cograph with a compatible set of triples)
	#

	if G._count_dG_notCograph_notConsistent > 0:
		cograph_freq = G._count_cographEdit_fixed_consistency/G._count_dG_notCograph_notConsistent
	else:
		cograph_freq = 0
	if G._count_dG_notCograph_notConsistent > 0:
		triples_freq = G._count_triplesEdit_fixed_cograph/G._count_dG_notCograph_notConsistent
	else:
		triples_freq = 0
	
	if G._count_dG_not_cograph > 0:
		cograph_to_LDT_freq = G._count_cographEdit_to_LDT/G._count_dG_not_cograph
	else:
		cograph_to_LDT_freq = 0
	if G._count_dG_not_consistent > 0:
		triples_to_LDT_freq = G._count_triplesEdit_to_LDT/G._count_dG_not_consistent
	else:
		triples_to_LDT_freq = 0




	edits_fixing_other = [cograph_freq, triples_freq]
	edits_to_LDT = [cograph_to_LDT_freq, triples_to_LDT_freq]
	bar_labels = ["cographEdit", "triplesEdit"]
	y_pos = np.arange(len(bar_labels))

	#fig, axs = plt.subplots(1, 2)

	plt.title("Frequency of cograph editing making the set of triples compatible and \ntriples editing turning the graph into a cograph")
	plt.bar(y_pos, edits_fixing_other, color=['blue', 'red'], width=[0.5, 0.5])
	plt.xticks(y_pos, bar_labels)
	plt.ylim(0, 1.05)

	#axs[0].bar(y_pos, edits_fixing_other, color=['blue', 'red'])
	#axs[1].bar(y_pos, edits_to_LDT, color=['blue', 'red'])


	plt.show()
	plt.title("Frequency of cograph/triples editing turning the graph into an LDT-graph")
	plt.bar(y_pos, edits_to_LDT, color=['blue', 'red'], width=[0.5, 0.5])
	plt.xticks(y_pos, bar_labels)
	plt.ylim(0, 1.05)
	plt.show()

def plot_success_data(G):
	pass