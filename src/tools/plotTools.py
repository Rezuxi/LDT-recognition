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

#G = InvestigateGraph(ldt)


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




def plot_frequencies(G, n, filename="cographEdit_and_weighted_mincut_triples_edit"):
	# -1 means there were no cases.

	# frequency of cograph editing turning the graph into a cograph
	if G._count_dG_not_cograph > 0:
		cograph_f1 = G._count_cographEdit_success / G._count_dG_not_cograph
	else:
		cograph_f1 = -0.1

	# frequency of cograph editing making the set of triples consistent. that is, the set of triples goes from inconsistent to consistent as a result of cograph editing.
	if G._count_dG_notCograph_notConsistent > 0:
		cograph_f2 = G._count_cographEdit_fixed_consistency/G._count_dG_notCograph_notConsistent
	else:
		cograph_f2 = -0.1

	# frequency of the set of triples remaining consistent after cograph editing
	if G._count_dG_notCograph_consistent > 0:
		cograph_f3 = G._count_cographEdit_remained_consistent / G._count_dG_notCograph_consistent
	else:
		cograph_f3 = -0.1

	# frequency of the graph remaining properly colored after cograph editing.
	if G._count_dG_not_cograph > 0:
		cograph_f4 = G._count_cographEdit_remain_properly_colored / G._count_dG_not_cograph
	else:
		cograph_f4 = -0.1

	# frequency of cograph editing turning the graph into an LDT-graph
	if G._count_dG_not_cograph > 0:
		cograph_f5 = G._count_cographEdit_to_LDT / G._count_dG_not_cograph
	else:
		cograph_f5 = -0.1



	# frequency of triples editing making the set of triples consistent.
	if G._count_dG_not_consistent > 0:
		triples_f1 = G._count_triplesEdit_success / G._count_dG_not_consistent
	else:
		triples_f1 = -0.1

	# frequency of triples editing turning the graph into a cograph.
	if G._count_dG_notCograph_notConsistent > 0:
		triples_f2 = G._count_triplesEdit_fixed_cograph / G._count_dG_notCograph_notConsistent
	else:
		triples_f2 = -0.1

	# frequency of graph remaining a cograph after triples ediiting.
	if G._count_dG_cograph_notConsistent > 0:
		triples_f3 = G._count_triplesEdit_remained_cograph / G._count_dG_cograph_notConsistent
	else:
		triples_f3 = -0.1

	# frequency of the graph remaining properly colored after triples editing.
	if G._count_dG_not_consistent > 0:
		triples_f4 = G._count_triplesEdit_remain_properly_colored / G._count_dG_not_consistent
	else:
		triples_f4 = -0.1

	# frequency of triples editing turning the graph into an LDT-graph.
	if G._count_dG_not_consistent > 0:
		triples_f5 = G._count_triplesEdit_to_LDT / G._count_dG_not_consistent
	else:
		triples_f5 = -0.1

	# get new ID for current filename
	new_name = "graphs/" + filename + "_6.png"

	# plot cograph and triples frequencies side by side.

	X = ["G1", "G2", "G3", "G4"]
	cograph_values = [cograph_f1, cograph_f2, cograph_f4, cograph_f5]
	triples_values = [triples_f1, triples_f2, triples_f4, triples_f5]
	colors = ['blue', 'green']
	X_axis = np.arange(len(X))
	
	plt.bar(X_axis-0.2, cograph_values, 0.4, label = 'Cograph edit')
	plt.bar(X_axis+0.2, triples_values, 0.4, label = 'Triples edit')

	plt.xticks(X_axis, X)
	plt.title("Cograph vs triples editing on 100 perturbed graphs with {} nodes".format(n))
	plt.ylim(0.1, 2.0)
	plt.legend()
	#plt.savefig(new_name)
	plt.show()

def plot_success_data(G):
	pass