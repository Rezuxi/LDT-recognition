from tools.GraphTools import *
from tools.plotTools import *
import networkx as nx
import asymmetree.treeevolve as te
from asymmetree.datastructures import PhyloTree
from asymmetree.hgt import ldt_graph

S = te.simulate_species_tree(10, model='innovation')
TGT = te.simulate_dated_gene_tree(S, dupl_rate = 0.5, loss_rate = 0.3, hgt_rate = 0.6, prohibit_extinction = "per_family", replace_prob=0.0)
OGT = te.observable_tree(TGT)
ldt = ldt_graph(OGT, S)



G = InvestigateGraph(ldt)
#print("edges of G: \n{}".format(G._G.edges()))
#a, b, c = get_P3_data(G._G)
#print("\nThe regions of P3s: {}".format(a))
#print("\nThe amounts in the regions: {}".format(b))
#print("\nThe distance between regions: {}\n".format(c))
print("amount of nodes in the ldt graph: {}".format(len(ldt.nodes())))


def run_investigation():

	for i in range(100):

		G.perturb_graph()

		'''
			COGRAPH EDIT
		'''
		if not G._is_cograph:

			G._count_dG_not_cograph += 1

			if not G._is_compatible:
				G._count_dG_notCograph_notConsistent += 1
			else:
				G._count_dG_notCograph_consistent += 1

			edited_G = G.cograph_editing()
			edited_G_is_cograph = is_cograph(edited_G)
			edited_G_is_properly_colored = is_properly_colored(edited_G, G._G)

			if edited_G_is_cograph:
				
				if edited_G_is_properly_colored:
					G._count_cographEdit_remain_properly_colored += 1
				edited_G_is_compatible = is_compatible(edited_G, G._G)
				#triples, leaves = find_all_P3(edited_G, get_triples=True, colored_G=G._G)
				#print("The set of triples for the cograph is: \n{}".format(triples))
				if edited_G_is_compatible and edited_G_is_properly_colored:
					G._count_cographEdit_to_LDT += 1
					#G.print_perturbed_G_data()
					#G.print_symmetric_diff(edited_G)
				else:
					# store some data about the non ldt graphs like density of the graph and something with P3s
					pass
				if edited_G_is_compatible and G._is_compatible:
					G._count_cographEdit_remained_consistent += 1
				elif edited_G_is_compatible and not G._is_compatible:
					G._count_cographEdit_fixed_consistency += 1
				elif not edited_G_is_compatible and G._is_compatible:
					G._count_cographEdit_broke_consistency += 1

		'''
			TRIPLES EDIT
		'''
		if not G._is_compatible:
			
			G._count_dG_not_consistent += 1

			if G._is_cograph:
				G._count_dG_cograph_notConsistent += 1

			edited_G, _ = G.triples_editing()

			if edited_G == None:
				# if the set of triples is empty, then it's considered consistent(?)
				# Build alg says it's inconsistent.
				pass
			else:

				edited_G_is_cograph = is_cograph(edited_G)

				if edited_G_is_cograph:
					G._count_triplesEdit_to_LDT += 1

				# LDT
				if edited_G_is_cograph and G._is_cograph:
					G._count_triplesEdit_remained_cograph += 1
				elif edited_G_is_cograph and not G._is_cograph:
					G._count_triplesEdit_fixed_cograph += 1
				# not LDT
				elif not edited_G_is_cograph and G._is_cograph:
					G._count_triplesEdit_broke_cograph += 1

	G.print_data()

	plot_data(G)
	'''
		To get the number of times a graph went from not a cograph and consistent to a cograph and still consistent
		we count (count_cographEdit_to_LDT - count_cographEdit_fixed_consistency)
	'''

	'''
		

	'''

	'''
		Also want to present the edit distance between perturbed graph and graph after edit.
		This is done using the symmetric difference tool.
	'''




if __name__ == "__main__":
	run_investigation()
	pass