from tools.GraphTools import *
from tools.plotTools import *
import networkx as nx
import asymmetree.treeevolve as te
from asymmetree.datastructures import PhyloTree
from asymmetree.hgt import ldt_graph
from tools.LDT_ILP import LDTEditor
import asymmetree.tools.GraphTools as gt
import os

S = te.simulate_species_tree(20, model='innovation')
TGT = te.simulate_dated_gene_tree(S, dupl_rate = 0.5, loss_rate = 0.5, hgt_rate = 0.5, prohibit_extinction = "per_family", replace_prob=0.0)
OGT = te.observable_tree(TGT)
ldt = ldt_graph(OGT, S)

colors = gt.sort_by_colors(ldt)


#print("edges of G: \n{}".format(G._G.edges()))
#a, b, c = get_P3_data(G._G)
#print("\nThe regions of P3s: {}".format(a))
#print("\nThe amounts in the regions: {}".format(b))
#print("\nThe distance between regions: {}\n".format(c))




print("Amount of nodes: {}".format(len(ldt.nodes())))
print("Amount of colors: {}".format(len(colors)))
print("Amount of edges: {}".format(len(ldt.edges())))

def run_investigation():
	successes = 0
	tested_graphs = 0
	m = 10
	G = InvestigateGraph(ldt)
	#if len(ldt.nodes()) == n:
	for i in range(m):
		perturbed = G.perturb_graph_terminate(0.4, 0.2)
		if not perturbed:
			continue
		tested_graphs += 1
		edited_to_ldt_graph = G.LDT_editing1(d= True)
		if edited_to_ldt_graph:
			successes += 1
			print("Edges in perturbed graph: {}".format(len(G._G_perturbed.edges())))
			print("Edges in edited ldt graph: {}".format(len(edited_to_ldt_graph.edges())))

	print("Amount of nodes: {}".format(len(ldt.nodes())))
	print("Amount of colors: {}".format(len(colors)))
	print("Amount of edges in the original ldt-graph: {}".format(len(ldt.edges())))
	if tested_graphs > 0:
		print("Frequency of triples+cograph editing turning G into an LDT-graph: {}".format(successes/tested_graphs))













if __name__ == "__main__":
	run_investigation()

	pass