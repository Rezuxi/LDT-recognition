from tools.GraphTools import *
from tools.plotTools import *
import networkx as nx
import asymmetree.treeevolve as te
from asymmetree.datastructures import PhyloTree
from asymmetree.hgt import ldt_graph
from tools.LDT_ILP import LDTEditor
import asymmetree.tools.GraphTools as gt
import os

S = te.simulate_species_tree(10, model='innovation')
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


def run_investigation():
	G = InvestigateGraph(ldt)
	global p3_average_triples
	for i in range(10):

		G.perturb_graph(0.5, 0.5)

		'''
			COGRAPH EDIT
		'''
		if not G._is_cograph:

			G._count_dG_not_cograph += 1

			if not G._is_compatible:
				G._count_dG_notCograph_notConsistent += 1
			else:
				G._count_dG_notCograph_consistent += 1

			# check that the edited graph is a cograph, properly colored and has a set of triples that is compatible.
			edited_G = G.cograph_editing()
			edited_G_is_cograph = is_cograph(edited_G)
			edited_G_is_properly_colored = is_properly_colored(edited_G, G._G)
			edited_G_is_compatible = is_compatible(edited_G, G._G)
			#print("Is edited cograph compatible? {}".format(edited_G_is_compatible))

			if edited_G_is_cograph:
				
				G._count_cographEdit_success += 1

				if edited_G_is_properly_colored:
					G._count_cographEdit_remain_properly_colored += 1

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
		#triples, _ = find_all_P3(edited_G, G._G)
		#triples2, _ = find_all_P3(G._G_perturbed)
		#print("Amount of P3s before cograph editing: {}".format(len(triples2)))
		#print("Amount of P3s in the edited cograph: {}".format(len(triples)))
		#print("The density of the edited graph is: {}".format(nx.density(edited_G)))
		#print("Is the edited cograph an LDT-graph? {}".format((edited_G_is_cograph and edited_G_is_properly_colored and edited_G_is_compatible)))

	
		'''
			TRIPLES EDIT
		'''
		# Always check after triples editing that the new set of triples is consistent.
		# and store the frequency of that.
		
		if not G._is_compatible:
			
			G._count_dG_not_consistent += 1

			if G._is_cograph:
				G._count_dG_cograph_notConsistent += 1

			edited_G, _ = G.triples_editing()

			# in this case, no edit was made since the set of triples was empty. the empty set is a consistent set of triples.
			if edited_G == None:
				edited_G_is_compatible = True
				edited_G_is_cograph = G._is_cograph
				edited_G_is_properly_colored = True
			else:
				edited_G_is_compatible = is_compatible(edited_G, G._G)
				edited_G_is_cograph = is_cograph(edited_G)
				edited_G_is_properly_colored = is_properly_colored(edited_G, G._G)

			if edited_G_is_properly_colored:
				G._count_triplesEdit_remain_properly_colored += 1

			if edited_G_is_compatible:
				G._count_triplesEdit_success += 1


			if edited_G_is_cograph and edited_G_is_compatible and edited_G_is_properly_colored:
				G._count_triplesEdit_to_LDT += 1

			
			if edited_G_is_cograph and G._is_cograph:
				G._count_triplesEdit_remained_cograph += 1
			elif edited_G_is_cograph and not G._is_cograph:
				G._count_triplesEdit_fixed_cograph += 1
			
			elif not edited_G_is_cograph and G._is_cograph:
				G._count_triplesEdit_broke_cograph += 1
		
		G._is_cograph = True
		G._is_compatible = True



	G.print_data()
	#plot_data(G)
	#plot_frequencies(G, len(G._G.nodes()))

def compare_P3(G, perturbed_G):
	"""
		Compare regions of P3s and their densities of LDT graphs vs non LDT-graphs
	"""
	regions, amounts, region_distances, P3_amount = get_P3_data(G)
	edited_regions, edited_amounts, edited_region_distances, edited_G_P3_amount = get_P3_data(perturbed_G)

	print("-----------------------------P3 data for LDT graph-----------------------------")
	print(regions)
	print(amounts)
	print(P3_amount)
	print(region_distances)

	print("-----------------------------P3 data for non LDT graph-----------------------------")
	print(edited_regions)
	print(edited_amounts)
	print(edited_G_P3_amount)
	print(edited_region_distances)

	

def ILP_solver(G):
	"""
		takes an LDT graph and perturbs it (ensuring it's no longer an ldt graph)
	"""
	print("Amount of nodes: {}".format(len(G.nodes())))
	IG = InvestigateGraph(G)
	IG.perturb_graph()

	solver = LDTEditor(IG._G_perturbed)
	solver.build_model()
	solver.optimize(time_limit=None)

	sol_graph, sol_distance = solver.get_solution()

	# if sol distance is 0 no edits were made and the graph is still not an LDT graph.
	# since the perturbed graph is guaranteed to not be an LDT graph

	properly_colored = is_properly_colored(sol_graph)
	cograph = is_cograph(sol_graph)
	compatible = is_compatible(sol_graph)

	edit_dist = gt.symmetric_diff(IG._G_perturbed, sol_graph)
	print("The value of the ILP: {}".format(sol_distance))
	print("The value of the symmetric distance: {}".format(edit_dist))
	print("Runtime: {}".format(solver.get_solve_time()))
	if properly_colored and cograph and compatible:
		print("Saving data...")
		solver._save_ILP_data(IG._G_perturbed, sol_graph, solver.get_solve_time(), edit_dist, only_add=False, only_delete=False)
	else:
		print("No solution found!")


def compare_edits_to_exact_results():
	path = "exact_results/"
	for file in os.listdir(path):
		filename = os.path.join(path, file)
		G, sol_graph, only_add, only_delete, min_edit_dist = LDTEditor.get_ILP_data(filename)
		# do cograph and triples editing, check if LDT-graph and compare edit distance
		IG = InvestigateGraph(sol_graph, G)
		
		cograph_edited_G = IG.cograph_editing()
		if is_cograph(cograph_edited_G) and is_properly_colored(cograph_edited_G, G) and is_compatible(cograph_edited_G, G):
			edit_dist = gt.symmetric_diff(cograph_edited_G, G)
			difference = gt.symmetric_diff(cograph_edited_G, sol_graph)
			print("\nThe graph G in the file '{}' was turned into an LDT-graph and the edit distance (of cograph editing) is : {}".format(filename, edit_dist))
			print("The edit distance between G and the exact result is: {}".format(min_edit_dist))
			print("The difference between the cograph edited graph and the exact solution is: {}".format(difference))
		else:
			print("\nThe graph G in the file '{}' was not turned into an LDT-graph as a result of cograph editing!".format(filename))


if __name__ == "__main__":
	run_investigation()
	#ILP_solver(ldt)
	#G = InvestigateGraph(ldt)
	#G.perturb_graph()
	#compare_P3(ldt, G._G_perturbed)
	#compare_edits_to_exact_results()
	pass