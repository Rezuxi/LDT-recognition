from tools.GraphTools import *
from tools.plotTools import *
import networkx as nx
import asymmetree.treeevolve as te
from asymmetree.datastructures import PhyloTree
from asymmetree.hgt import ldt_graph
from tools.LDT_ILP import LDTEditor
import asymmetree.tools.GraphTools as gt
import os
import copy
import matplotlib.pyplot as plt


# use ILP solver on 100 perturbed graphs with n nodes for different n values
# save in different folder for different n values.

def generate_solutions(n):
	done = False
	while not done:
		S = te.simulate_species_tree(10, model='innovation')
		TGT = te.simulate_dated_gene_tree(S, dupl_rate = 0.5, loss_rate = 0.5, hgt_rate = 0.5, prohibit_extinction = "per_family", replace_prob=0.0)
		OGT = te.observable_tree(TGT)
		ldt = ldt_graph(OGT, S)
		if len(ldt.nodes()) == n:
			for i in range(100):
				IG = InvestigateGraph(ldt)
				IG.perturb_graph(0.5, 0.5)

				solver = LDTEditor(IG._G_perturbed)
				solver.build_model()
				solver.optimize(time_limit=None)

				sol_graph, sol_distance = solver.get_solution()

				properly_colored = is_properly_colored(sol_graph)
				cograph = is_cograph(sol_graph)
				compatible = is_compatible(sol_graph)

				edit_dist = gt.symmetric_diff(IG._G_perturbed, sol_graph)
				print("Runtime: {}".format(solver.get_solve_time()))
				if properly_colored and cograph and compatible:
					print("Saving data...")
					solver._save_ILP_data(IG._G_perturbed, sol_graph, solver.get_solve_time(), edit_dist, only_add=False, only_delete=False, filename="{}nodes/LDTEdit_exact_solution".format(n))
				else:
					print("No solution found!")
			done = True


def cograph_editing(IG):
	if not IG._is_cograph:

		IG._count_dG_not_cograph += 1

		if not IG._is_compatible:
			# dont need to count this here since its being counted in triples edit and the same IG is used for cograph and triples editing
			#IG._count_dG_notCograph_notConsistent += 1
			pass
		else:
			IG._count_dG_notCograph_consistent += 1

		# check that the edited graph is a cograph, properly colored and has a set of triples that is compatible.
		edited_G = IG.cograph_editing()
		edited_G_is_cograph = is_cograph(edited_G)
		edited_G_is_properly_colored = is_properly_colored(edited_G, IG._G)
		edited_G_is_compatible = is_compatible(edited_G, IG._G)


		if edited_G_is_cograph:
			
			IG._count_cographEdit_success += 1

			if edited_G_is_properly_colored:
				IG._count_cographEdit_remain_properly_colored += 1

			if edited_G_is_compatible and IG._is_compatible:
				IG._count_cographEdit_remained_consistent += 1
			elif edited_G_is_compatible and not IG._is_compatible:
				IG._count_cographEdit_fixed_consistency += 1
			elif not edited_G_is_compatible and IG._is_compatible:
				IG._count_cographEdit_broke_consistency += 1

			if edited_G_is_compatible and edited_G_is_properly_colored:
				IG._count_cographEdit_to_LDT += 1
				number_edges_remaining = len(edited_G.edges())
				edit_dist = gt.symmetric_diff(IG._G_perturbed, edited_G)
				return edited_G, True, number_edges_remaining, edit_dist

	return None, False, None, None


def triples_editing(IG, a):
	if not IG._is_compatible:
		IG._count_dG_not_consistent += 1

		if IG._is_cograph:
			IG._count_dG_cograph_notConsistent += 1
		else:
			IG._count_dG_notCograph_notConsistent += 1

		if a == 1:
			edited_G, _ = IG.triples_editing()
		elif a == 2:
			edited_G, _ = IG.triples_editing(destroy_K3 = 1)	
		elif a == 3:
			edited_G, _ = IG.triples_editing2()

		if edited_G == None:
			edited_G_is_compatible = True
			edited_G_is_cograph = IG._is_cograph
			edited_G_is_properly_colored = True
		else:
			edited_G_is_compatible = is_compatible(edited_G, IG._G)
			edited_G_is_cograph = is_cograph(edited_G)
			edited_G_is_properly_colored = is_properly_colored(edited_G, IG._G)

		if edited_G_is_properly_colored:
			IG._count_triplesEdit_remain_properly_colored += 1

		if edited_G_is_compatible:
			IG._count_triplesEdit_success += 1

		if edited_G_is_cograph and IG._is_cograph:
			IG._count_triplesEdit_remained_cograph += 1
		elif edited_G_is_cograph and not IG._is_cograph:
			IG._count_triplesEdit_fixed_cograph += 1	
		elif not edited_G_is_cograph and IG._is_cograph:
			IG._count_triplesEdit_broke_cograph += 1

		if edited_G_is_cograph and edited_G_is_compatible and edited_G_is_properly_colored:
			IG._count_triplesEdit_to_LDT += 1
			number_edges_remaining = len(edited_G.edges())
			edit_dist = gt.symmetric_diff(IG._G_perturbed, edited_G)
			return edited_G, True, number_edges_remaining, edit_dist
	return None, False, None, None
		
# G1: frequency of G becoming a cograph as a result of the edit
# G2: frequency of G becoming consistent as a result of the edit
# G3: frequency of G remaining properly colored after the edit
# G4: frequency of G becoming an LDT-graph as a result of the edit
def get_freq(G):
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

	# frequency of the graph remaining properly colored after cograph editing.
	if G._count_dG_not_cograph > 0:
		cograph_f3 = G._count_cographEdit_remain_properly_colored / G._count_dG_not_cograph
	else:
		cograph_f3 = -0.1

	# frequency of cograph editing turning the graph into an LDT-graph
	if G._count_dG_not_cograph > 0:
		cograph_f4 = G._count_cographEdit_to_LDT / G._count_dG_not_cograph
	else:
		cograph_f4 = -0.1



	# frequency of triples editing turning the graph into a cograph.
	if G._count_dG_notCograph_notConsistent > 0:
		triples_f1 = G._count_triplesEdit_fixed_cograph / G._count_dG_notCograph_notConsistent
	else:
		triples_f1 = -0.1

	# frequency of triples editing making the set of triples consistent.
	if G._count_dG_not_consistent > 0:
		triples_f2 = G._count_triplesEdit_success / G._count_dG_not_consistent
	else:
		triples_f2 = -0.1

	# frequency of the graph remaining properly colored after triples editing.
	if G._count_dG_not_consistent > 0:
		triples_f3 = G._count_triplesEdit_remain_properly_colored / G._count_dG_not_consistent
	else:
		triples_f3 = -0.1

	# frequency of triples editing turning the graph into an LDT-graph.
	if G._count_dG_not_consistent > 0:
		triples_f4 = G._count_triplesEdit_to_LDT / G._count_dG_not_consistent
	else:
		triples_f4 = -0.1

	triples_values = [triples_f1, triples_f2, triples_f3, triples_f4]
	cograph_values = [cograph_f1, cograph_f2, cograph_f3, cograph_f4]

	
	return triples_values, cograph_values


def get_y_values(bool_list, values_list):
	y_values = []
	for i in range(len(bool_list)):
		if bool_list[i]:
			y_values.append(values_list[i])
		else:
			y_values.append(float('NaN'))
	return y_values

def bar_dot_plot(c_f, t1_f, t2_f, t3_f, m, n, initial_num_edges, exact_sol_num_edges, exact_sol_min_dist,
											c_edge_count, c_edit_dist, c_is_ldt,
											t1_edge_count, t1_edit_dist, t1_is_ldt,
											t2_edge_count, t2_edit_dist, t2_is_ldt,
											t3_edge_count, t3_edit_dist, t3_is_ldt
											):

	X_bar = ["G1", "G2", "G3", "G4"]
	X_axis = np.arange(len(X_bar))

	# BAR PLOT
	plt.bar(X_axis + 0.00, c_f, label = 'Cograph editing', width = 0.20)
	plt.bar(X_axis + 0.20, t1_f, label = 'Triples editing 1', width = 0.20)
	plt.bar(X_axis + 0.40, t2_f, label = 'Triples editing 2', width = 0.20)
	plt.bar(X_axis + 0.60, t3_f, label = 'Triples editing 3', width = 0.20)

	plt.xticks(X_axis, X_bar)
	plt.title("Cograph vs triples editing on {} perturbed graphs with {} nodes".format(m, n))
	plt.ylim(0.1, 2.0)
	plt.legend()
	plt.show()

	x = np.arange(1, 101, 1)
	y_c_edge_count = get_y_values(c_is_ldt, c_edge_count)
	y_t1_edge_count = get_y_values(t1_is_ldt, t1_edge_count)
	y_t2_edge_count = get_y_values(t2_is_ldt, t2_edge_count)
	y_t3_edge_count = get_y_values(t3_is_ldt, t3_edge_count)

	y_init_edge_count = []
	y_min_edit_dist = []
	y_exact_sol_edge_count = []
	for i in range(100):
		if c_is_ldt[i] or t1_is_ldt[i] or t2_is_ldt[i] or t3_is_ldt[i]:
			y_init_edge_count.append(initial_num_edges[i])
			y_min_edit_dist.append(exact_sol_min_dist[i])
			y_exact_sol_edge_count.append(exact_sol_num_edges[i])
		else:
			y_init_edge_count.append(float('NaN'))
			y_min_edit_dist.append(float('NaN'))
			y_exact_sol_edge_count.append(float('NaN'))

	# DOT PLOT (EDGE COUNT)
	plt.rcParams["figure.figsize"] = (15, 6)
	plt.title("Edge count of input graph and edited LDT-graphs    {} nodes".format(n))
	plt.plot(x, y_exact_sol_edge_count, 's', color='green', label = 'exact solution graph')
	plt.plot(x, y_init_edge_count, 'o', color='brown', label = 'input graph')
	plt.plot(x, y_c_edge_count, 'o', color='black', label = 'cograph edited LDT-graph')
	plt.plot(x, y_t1_edge_count, 'o', color='purple', label = 'triples edited (1) LDT-graph')
	plt.plot(x, y_t2_edge_count, 'o', color='red', label = 'triples edited (2) LDT-graph')
	plt.plot(x, y_t3_edge_count, 'o', color='orange', label = 'triples edited (3) LDT-graph')


	plt.xticks(np.arange(min(x), max(x)+1, 2.0))
	plt.ylabel("edge count")
	plt.xlabel("graph No.")
	plt.legend(bbox_to_anchor=(1.10, 1), loc='upper right', fontsize='small')
	plt.tight_layout()
	plt.show()

	y_c_edit_dist = get_y_values(c_is_ldt, c_edit_dist)
	y_t1_edit_dist = get_y_values(t1_is_ldt, t1_edit_dist)
	y_t2_edit_dist = get_y_values(t2_is_ldt, t2_edit_dist)
	y_t3_edit_dist = get_y_values(t3_is_ldt, t3_edit_dist)

	# DOT PLOT (EDIT DISTANCE)
	plt.plot(x, y_min_edit_dist, 's', color='green', label = 'exact solution graph')
	plt.plot(x, y_c_edit_dist, 'o', color='black', label = 'cograph edited LDT-graph')
	plt.plot(x, y_t1_edit_dist, 'o', color='purple', label = 'triples edited (1) LDT-graph')
	plt.plot(x, y_t2_edit_dist, 'o', color='red', label = 'triples edited (2) LDT-graph')
	plt.plot(x, y_t3_edit_dist, 'o', color='orange', label = 'triples edited (3) LDT-graph')

	plt.xticks(np.arange(min(x), max(x)+1, 2.0))
	plt.ylabel("edit distance (symmetric difference)")
	plt.xlabel("graph No.")
	plt.legend(bbox_to_anchor=(1.10, 1), loc='upper right', fontsize='small')
	plt.tight_layout()
	plt.show()

def benchmark(n):

	c1_graphs, c1_edge_count, c1_is_ldt, c1_edit_dist = ([] for i in range(4))
	t1_graphs, t1_edge_count, t1_is_ldt, t1_edit_dist = ([] for i in range(4))
	t2_graphs, t2_edge_count, t2_is_ldt, t2_edit_dist = ([] for i in range(4))
	t3_graphs, t3_edge_count, t3_is_ldt, t3_edit_dist = ([] for i in range(4))

	init_edge_amounts = []
	exact_sol_edge_amounts = []

	exact_sol_min_dist = []

	IG1 = None
	files = os.listdir("exact_results/{}nodes".format(n))
	for f in files:

		G, edited_G, only_add, only_delete, min_edit_dist = LDTEditor.get_ILP_data("exact_results/{}nodes/".format(n) + f)	
		if not IG1:
			IG1 = InvestigateGraph(edited_G, disturbed_G = G)
			IG2 = copy.deepcopy(IG1)
			IG3 = copy.deepcopy(IG1)
		else:
			IG1.set_perturbed_graph(G)
			IG1.set_G(edited_G)

			IG2.set_perturbed_graph(G)
			IG2.set_G(edited_G)

			IG3.set_perturbed_graph(G)
			IG3.set_G(edited_G)

		initial_num_edges = len(G.edges())
		exact_sol_num_edges = len(edited_G.edges())

		cograph_edited_G, is_c1_ldt, c1_num_edges, c1_ldt_edit_dist = cograph_editing(IG1)
		triples1_edited_G, is_t1_ldt, t1_num_edges, t1_ldt_edit_dist = triples_editing(IG1, 1)
		triples2_edited_G, is_t2_ldt, t2_num_edges, t2_ldt_edit_dist = triples_editing(IG2, 2)
		triples3_edited_G, is_t3_ldt, t3_num_edges, t3_ldt_edit_dist = triples_editing(IG3, 3)

		exact_sol_min_dist.append(min_edit_dist)

		init_edge_amounts.append(initial_num_edges)
		exact_sol_edge_amounts.append(exact_sol_num_edges)

		c1_graphs.append(cograph_edited_G)
		t1_graphs.append(triples1_edited_G)
		t2_graphs.append(triples2_edited_G)
		t3_graphs.append(triples3_edited_G)

		c1_is_ldt.append(is_c1_ldt)
		t1_is_ldt.append(is_t1_ldt)
		t2_is_ldt.append(is_t2_ldt)
		t3_is_ldt.append(is_t3_ldt)

		c1_edge_count.append(c1_num_edges)
		t1_edge_count.append(t1_num_edges)
		t2_edge_count.append(t2_num_edges)
		t3_edge_count.append(t3_num_edges)

		c1_edit_dist.append(c1_ldt_edit_dist)
		t1_edit_dist.append(t1_ldt_edit_dist)
		t2_edit_dist.append(t2_ldt_edit_dist)
		t3_edit_dist.append(t3_ldt_edit_dist)

	triples1_freq, cograph_freq = get_freq(IG1)
	triples2_freq, _ 			= get_freq(IG2)
	triples3_freq, _ 			= get_freq(IG3)

	
	# create bar plot with scatter plot to the right with edit dist and edge count for each graph that became LDT graph.
	bar_dot_plot(cograph_freq, triples1_freq, triples2_freq, triples3_freq, 
				100, n, init_edge_amounts, exact_sol_edge_amounts, exact_sol_min_dist,
				c1_edge_count, c1_edit_dist, c1_is_ldt,
				t1_edge_count, t1_edit_dist, t1_is_ldt,
				t2_edge_count, t2_edit_dist, t2_is_ldt,
				t3_edge_count, t3_edit_dist, t3_is_ldt
				)

		

#benchmark(15)
generate_solutions(20)