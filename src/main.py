from tools.GraphTools import *
import networkx as nx
import asymmetree.treeevolve as te
from asymmetree.datastructures import PhyloTree
from asymmetree.hgt import ldt_graph

S = te.simulate_species_tree(10, model='innovation')
TGT = te.simulate_dated_gene_tree(S, dupl_rate = 0.5, loss_rate = 0.5, hgt_rate = 0.5, prohibit_extinction = "per_family", replace_prob=1.0)
OGT = te.observable_tree(TGT)
ldt = ldt_graph(OGT, S)

G = InvestigateGraph(ldt)

def run_investigation():
	# do n times
		# load new species and gene tree
		# create LDT graph
		# save data for this LDT graph
		# do x times
			# perturb LDT graph
			# do cograph edit
			# do triples edit
			# save data in some way (perhaps plot)
		# compare edited graph to LDT graph data
	G.perturb_graph()

	if not G._is_cograph:
		print("Not cograph")
		edited_G = G.cograph_editing()
		edited_G_cograph = is_cograph(edited_G)
		if edited_G_cograph:
			print("Edited to cograph!")
	if not G._is_compatible:
		print("Not compatible")

if __name__ == "__main__":
	run_investigation()