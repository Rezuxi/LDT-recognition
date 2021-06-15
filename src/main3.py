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





# 'exact_solutions/trees/{}trees'
def generate_solutions_fromTrees(n, filename):
	# load species+gene trees
	name = filename + '/{}trees'.format(n)
	tree_files = []

	probabilities = [0.15, 0.30, 0.50]
	#nodes = [10, 14, 18]
	restrictions = ['', 'insertion', 'deletion']

	for _, _, files in os.walk(name):
		for file in files:
			tree_files.append(file)
	species_trees = []
	gene_trees 	  = []

	for f in tree_files:
		if f.startswith('species'):
			species_trees.append(f)
		else:
			gene_trees.append(f)
	ID = 0
	for i in range(len(species_trees)):
		S = PhyloTree.load(name + '/{}'.format(species_trees[i]))
		TGT = PhyloTree.load(name + '/{}'.format(gene_trees[i]))
		OGT = te.observable_tree(TGT)
		ldt = ldt_graph(OGT, S)
		# perturb using p = 0.15, 0.3, 0.5
		# for each p, solve using ILP with deletion, insertion and both
		for p1 in probabilities:
			for p2 in probabilities:
				p_i = str(p1).replace('.', '')
				if len(p_i) < 3:
					p_i = p_i + '0'
				p_d = str(p2).replace('.', '')
				if len(p_d) < 3:
					p_d = p_d + '0'
				IG = InvestigateGraph(ldt)
				perturbed = IG.perturb_graph_terminate(p1, p2)
				if not perturbed:
					print("failed")
				else:
					# solve 3 times using deletion, insertion and both
					solver1 = LDTEditor(IG._G_perturbed)
					solver1.build_model()
					solver1.optimize(time_limit=None)

					solver2 = LDTEditor(IG._G_perturbed, only_delete = True)
					solver2.build_model()
					solver2.optimize(time_limit=None)
					
					solver3 = LDTEditor(IG._G_perturbed, only_add = True)
					solver3.build_model()
					solver3.optimize(time_limit=None)

					sol_graph1, sol_distance1 = solver1.get_solution()
					sol_graph2, sol_distance2 = solver2.get_solution()
					sol_graph3, sol_distance3 = solver3.get_solution()

					properly_colored1 = is_properly_colored(sol_graph1)
					cograph1 = is_cograph(sol_graph1)
					compatible1 = is_compatible(sol_graph1)

					properly_colored2 = is_properly_colored(sol_graph2)
					cograph2 = is_cograph(sol_graph2)
					compatible2 = is_compatible(sol_graph2)

					properly_colored3 = is_properly_colored(sol_graph3)
					cograph3 = is_cograph(sol_graph3)
					compatible3 = is_compatible(sol_graph3)

					folderName = 'exact_solutions/{}_{}_{}nodes{}/'
					saveFolder1 = folderName.format(p1, p2, n, '')

					if properly_colored1 and cograph1 and compatible1:
						print("Saving data...")
						solver1._save_ILP_data(IG._G_perturbed, sol_graph1, solver1.get_solve_time(), sol_distance1, i_p = p1, d_p = p2, only_add=False, only_delete=False, saveFolder = folderName.format(p_i, p_d, n, ''), ID = ID)
					else:
						print("No solution found!")

					if properly_colored2 and cograph2 and compatible2:
						print("Saving data (deletion)...")
						solver2._save_ILP_data(IG._G_perturbed, sol_graph2, solver2.get_solve_time(), sol_distance2, i_p = p1, d_p = p2, only_add=False, only_delete=True, saveFolder = folderName.format(p_i, p_d, n, '_deletion'), ID = ID)
					else:
						print("No solution found for deletion only!")

					if properly_colored3 and cograph3 and compatible3:
						print("Saving data (insertion)...")
						solver3._save_ILP_data(IG._G_perturbed, sol_graph3, solver3.get_solve_time(), sol_distance3, i_p = p1, d_p = p2, only_add=True, only_delete=False, saveFolder = folderName.format(p_i, p_d, n, '_insertion'), ID = ID)
					else:
						print("No solution found for insertion only!")
					ID += 1


def generate_folders():
	ps = ['015', '030', '050']
	n = [10, 14, 18]
	for p in ps:
		for p2 in ps:	
			for i in n:
				#print("{}, {}, {}".format(p, p2, i))
				path1 = 'exact_solutions/{}_{}_{}nodes'.format(p, p2, i)
				path2 = 'exact_solutions/{}_{}_{}nodes_deletion'.format(p, p2, i)
				path3 = 'exact_solutions/{}_{}_{}nodes_insertion'.format(p, p2, i)

				os.mkdir(path1)
				os.mkdir(path2)
				os.mkdir(path3)



def count_files(n, filename):
	name = filename + '_{}nodes'.format(n)
	all_files = []
	for _, _, files in os.walk(name):
		for file in files:
			all_files.append(file)
	return len(all_files)


def triples_editing(IG, n = 1, deletion = False, insertion = False):
	if not IG._is_compatible:
		IG._count_dG_not_consistent += 1

		if IG._is_cograph:
			IG._count_dG_cograph_notConsistent += 1
		else:
			IG._count_dG_notCograph_notConsistent += 1


		edited_G, triples = IG.triples_editing(n = n, deletion = deletion, insertion = insertion)
		

		if triples == None:
			edited_G_is_compatible = True
			edited_G_is_cograph = is_cograph(edited_G)
			edited_G_is_properly_colored = is_properly_colored(edited_G, IG._G)
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


def cograph_editing(IG):
	if not IG._is_cograph:

		IG._count_dG_not_cograph += 1

		if not IG._is_compatible:
			IG._count_dG_notCograph_notConsistent += 1
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


def LDT_editing(IG, n = 1, deletion = False, insertion = False):
	IG._count_not_ldt += 1
	triples_edited_G, _ = IG.triples_editing(n = n, deletion = deletion, insertion = insertion)
	is_properly_colored = True
	if (deletion ^ insertion):
		# only insert or delete so we are sure to make G consistent.
		isConsistent = True
	else:
		isConsistent = is_compatible(triples_edited_G)
	isCograph = is_cograph(triples_edited_G)

	if not isCograph:
		cograph_edited_G = IG.cograph_editing(G = triples_edited_G)
	else:
		IG._count_ldtEdit_success += 1
		edges_remaining = len(triples_edited_G.edges())
		edit_dist = gt.symmetric_diff(IG._G_perturbed, triples_edited_G)
		return triples_edited_G, True, edges_remaining, edit_dist
		

	color_graph(IG._G, cograph_edited_G)
	properClrd_cograph = make_properly_colored(cograph_edited_G)
	isCograph = is_cograph(properClrd_cograph)
	isConsistent = is_compatible(properClrd_cograph)
	if isConsistent and isCograph:
		IG._count_ldtEdit_success += 1
		edges_remaining = len(properClrd_cograph.edges())
		edit_dist = gt.symmetric_diff(IG._G_perturbed, properClrd_cograph)
		return properClrd_cograph, True, edges_remaining, edit_dist
	return None, False, None, None

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

	# for ldt editing we're only interested in how often we get an ldt graph.
	if G._count_not_ldt > 0:
		ldt_f = G._count_ldtEdit_success / G._count_not_ldt
	else:
		ldt_f = -0.1


	triples_values = [triples_f1, triples_f2, triples_f3, triples_f4]
	cograph_values = [cograph_f1, cograph_f2, cograph_f3, cograph_f4]
	ldt_values 	   = [-0.1, -0.1, -0.1, ldt_f]
	#ldt_values	   = []
	
	return triples_values, cograph_values, ldt_values

def get_ratios(min_dist, edit_dist):
	ratios = []
	for i in range(len(min_dist)):
		if edit_dist[i]:
			ratios.append(edit_dist[i]/min_dist[i])
	return ratios

def box_plots(n, p1, p2, ldt_edit_distances, cograph_edit_distances, t1_edit_distances, t2_edit_distances, t3_edit_distances):

	data = [ldt_edit_distances, cograph_edit_distances, t1_edit_distances, t2_edit_distances, t3_edit_distances]
	labels = ['LDT editing', 'cograph editing', 'triples editing', 'triples editing (deletion)', 'triples editing (insertion)']

	#fig = plt.figure(figsize =(10, 7))
	
	# Creating plot
	#plt.boxplot(ldt_edit_distances)
	#plt.boxplot(cograph_edit_distances)
	#plt.boxplot(t1_edit_distances)
	#plt.boxplot(t2_edit_distances)
	#plt.boxplot(t3_edit_distances)
	plt.title("Ratio between edit distance and minimum edit distance for each heuristic.")
	plt.ylabel("ratio between edit distance and minimum edit distance")
	plt.xlabel("Different heuristics")
	plt.boxplot(data, positions = np.array(range(len(labels)))*2.0, sym='', widths=0.6)
	plt.xticks(range(0, len(labels) * 2, 2), labels, rotation=30)
	# show plot
	plt.tight_layout()
	plt.show()

def bar_plots(n, p1, p2, cograph_f, t1_f, t2_f, t3_f, ldt_f):
	S1 = "$G^*$ became a cograph"
	S2 = "$R_{G^*}$ became consistent"
	S3 = "$G^*$ remained properly colored"
	S4 = "$G^*$ became an LDT-graph"

	X_bar = [S1, S2, S3, S4]
	X_axis = np.arange(len(X_bar))
	f_size = 18
	#plt.figure(figsize =(7, 5))
	# Y axis : Frequency
	# X axis : G* -> cograph, R_G* -> consistent, G* properly colored, G* -> LDT

	# COLORS FOR ILP: 'green', INPUT GRAPH: 'brown'
	# COLORS FOR EDITS: COGRAPH: 'blue', TE1: 'purple', TE2: '#1a9c79', TE3: 'red', TE4: 'orange'
	# BAR PLOT
	
	plt.bar(X_axis + 0.00, cograph_f, label = 'cograph editing', width = 0.15, color = 'blue')
	plt.bar(X_axis + 0.15, t1_f, label = 'Triples editing', width = 0.15, color = '#1a9c79')
	plt.bar(X_axis + 0.30, t2_f, label = 'Triples editing (deletion)', width = 0.15, color = '#6a12b3')
	plt.bar(X_axis + 0.45, t3_f, label = 'Triples editing (insertion)', width = 0.15, color = 'red')
	#plt.bar(X_axis + 0.60, ldt_f, label = 'LDT editing', width = 0.15, color = 'orange')
	
	
	
	plt.ylabel("frequency", fontsize=f_size)
	plt.xticks(X_axis, X_bar, rotation = 30, fontsize=f_size)
	plt.title("different heuristics applied to 100 perturbed lDT-graphs with {} vertices.\nThese graphs were perturbed with $p=({}, {})$.".format(n, p1, p2), fontsize = f_size)
	plt.ylim(0.0, 1.0)
	plt.legend(bbox_to_anchor=(1.00, 1.00), loc="upper left", borderaxespad=0, prop={'size': 12}, fontsize=f_size)
	#plt.legend()
	plt.tight_layout()
	plt.show()


# get ldt graphs and perturb with different probabilities and apply triples editing with insertion and deletion, then compare edge count

def benchmarkLDTEdits(n, P1, P2):
	p1 = str(P1).replace('.', '')
	p1 = p1 + '0' if len(p1) == 2 else p1
	
	p2 = str(P2).replace('.', '')
	p2 = p2 + '0' if len(p2) == 2 else p2

	t1_graphs, t1_edge_count, t1_is_ldt, t1_edit_dist = ([] for i in range(4))	# ldt editing	
	t2_graphs, t2_edge_count, t2_is_ldt, t2_edit_dist = ([] for i in range(4))	# ldt editing (triples edit deletion)
	t3_graphs, t3_edge_count, t3_is_ldt, t3_edit_dist = ([] for i in range(4))	# ldt editing (triples edit insertion)		


	files = os.listdir("exact_solutions/{}_{}_{}nodes".format(p1, p2, n))
	files_del = os.listdir("exact_solutions/{}_{}_{}nodes_deletion".format(p1, p2, n))
	files_ins = os.listdir("exact_solutions/{}_{}_{}nodes_insertion".format(p1, p2, n))

	exact_sol_min_dist = []
	exact_sol_min_dist_del = []
	exact_sol_min_dist_ins = []


	IG1 = None

	for i in range(100):
		# G1 (perturbed graph is the same for these 3 lines G1=G2=G3)
		G1, edited_G1, only_add1, only_delete1, min_edit_dist1 = LDTEditor.get_ILP_data("exact_solutions/{}_{}_{}nodes/".format(p1, p2, n) + files[i])
		G2, edited_G2, only_add2, only_delete2, min_edit_dist2 = LDTEditor.get_ILP_data("exact_solutions/{}_{}_{}nodes_deletion/".format(p1, p2, n) + files_del[i])	
		G3, edited_G3, only_add3, only_delete3, min_edit_dist3 = LDTEditor.get_ILP_data("exact_solutions/{}_{}_{}nodes_insertion/".format(p1, p2, n) + files_ins[i])	
		
		if not IG1:
			IG1 = InvestigateGraph(edited_G1, disturbed_G = G1)
			IG2 = copy.deepcopy(IG1)
			IG3 = copy.deepcopy(IG1)
		else:
			IG1.set_perturbed_graph(G1)
			IG1.set_G(edited_G1)

			IG2.set_perturbed_graph(G1)
			IG2.set_G(edited_G1)

			IG3.set_perturbed_graph(G1)
			IG3.set_G(edited_G1)

		t1_edited_G, is_t1_ldt, _, t1_ldt_edit_dist = LDT_editing(IG1, n=100)	# ldt editing with triples editing allowing both deletions and insertions for n = 100.
		t2_edited_G, is_t2_ldt, _, t2_ldt_edit_dist = LDT_editing(IG2, deletion=True)	# ldt editing with triples editing allowing only deletions.
		t3_edited_G, is_t3_ldt, _, t3_ldt_edit_dist = LDT_editing(IG3, insertion=True)	# ldt editing with triples editing allowing only insertions.
		
		t1_graphs.append(t1_edited_G)
		t2_graphs.append(t2_edited_G)
		t3_graphs.append(t3_edited_G)

		t1_is_ldt.append(is_t1_ldt)
		t2_is_ldt.append(is_t2_ldt)
		t3_is_ldt.append(is_t3_ldt)

		t1_edit_dist.append(t1_ldt_edit_dist)
		t2_edit_dist.append(t2_ldt_edit_dist)
		t3_edit_dist.append(t3_ldt_edit_dist)

		exact_sol_min_dist.append(min_edit_dist1)
		exact_sol_min_dist_del.append(min_edit_dist2)
		exact_sol_min_dist_ins.append(min_edit_dist3)

	# get frequencies
	_, _, ldt1_freq				= get_freq(IG1)
	_, _, ldt2_freq				= get_freq(IG2)
	_, _, ldt3_freq				= get_freq(IG3)

	# get ratios.
	# we compare with min dist with no restrictions for all three methods since we only restrict triples editing to deletion/insertion
	ldt1_edit_ratios = get_ratios(exact_sol_min_dist, t1_edit_dist)
	ldt2_edit_ratios = get_ratios(exact_sol_min_dist, t2_edit_dist)	
	ldt3_edit_ratios = get_ratios(exact_sol_min_dist, t3_edit_dist) 


	
	font_size = 18
	bar_title = "Frequency of LDT editing (with different restrictions for triples editing)\nbeing successful on 100 perturbed LDT-graphs with {} vertices.".format(n) 
	bar_title =	bar_title + "These graphs were perturbed with $p=({}, {})$.".format(P1, P2)
	box_title = "Ratios $C_i/X_i$ for $i=1,...,100$ where $X_i$ is the (min) edit distance between $G_i$ and $G^{'}_i$\nand $C_i$ is the edit distance between $G^*_i$ and $G_i$."
	box_title = box_title + "These graphs have {} vertices and were obtained by perturbing LDT-graphs with $p=({}, {})$.".format(n, P1, P2)
	
	y_label_bar = "Frequency"
	y_label_box = "Ratio"

	x_rotation = 30
	x_ticks = ["LDT editing (i)", "LDT editing (ii)", "LDT editing (iii)"]

	X_axis = np.arange(len(x_ticks))

	# bar plot for frequency of ldt edits being successful
	bar_data = [ldt1_freq[3], ldt2_freq[3], ldt3_freq[3]]
	plt.title(bar_title, fontsize=font_size)
	plt.ylabel(y_label_bar, fontsize=font_size)
	plt.ylim(0.0, 1.0)
	plt.xticks(X_axis, x_ticks, rotation = 30, fontsize=font_size)
	plt.bar(x_ticks, bar_data, width=0.20)
	#plt.bar(X_axis + 0.00, ldt1_freq[3], width = 0.15)
	#plt.bar(X_axis + 0.15, ldt2_freq[3], width = 0.15)
	#plt.bar(X_axis + 0.30, ldt3_freq[3], width = 0.15)
	plt.tight_layout()
	plt.show()

	# box plot
	box_data = [ldt1_edit_ratios, ldt2_edit_ratios, ldt3_edit_ratios]
	plt.title(box_title, fontsize=font_size)
	plt.ylabel(y_label_box, fontsize=font_size)
	plt.boxplot(box_data, positions = np.array(range(len(x_ticks)))*2.0, sym='', widths=0.6)
	plt.xticks(range(0, len(x_ticks) * 2, 2), x_ticks, rotation=30, fontsize=font_size)
	plt.tight_layout()
	plt.show()
	

	

def benchmark(n, p1, p2):
	c1_graphs, c1_edge_count, c1_is_ldt, c1_edit_dist = ([] for i in range(4))	# cograph editing
	t1_graphs, t1_edge_count, t1_is_ldt, t1_edit_dist = ([] for i in range(4))	# triples editing with both insertion/deletion
	t2_graphs, t2_edge_count, t2_is_ldt, t2_edit_dist = ([] for i in range(4))	# triples editing with deletion only
	t3_graphs, t3_edge_count, t3_is_ldt, t3_edit_dist = ([] for i in range(4))	# triples editing with insertion only
	t4_graphs, t4_edge_count, t4_is_ldt, t4_edit_dist = ([] for i in range(4))	# ldt editing

	init_edge_amounts = []
	exact_sol_edge_amounts = []
	exact_sol_min_dist = []

	exact_sol_edge_amounts_del = []
	exact_sol_min_dist_del = []

	exact_sol_edge_amounts_ins = []
	exact_sol_min_dist_ins = []

	IG1 = None
	#files1 = os.listdir("exact_solutions/050_050_{}nodes".format(n))
	#files1_del = os.listdir("exact_solutions/050_050_{}nodes_deletion".format(n))
	#files1_ins = os.listdir("exact_solutions/050_050_{}nodes_insertion".format(n))

	files1 = os.listdir("exact_solutions/{}_{}_{}nodes".format(p1, p2, n))
	files1_del = os.listdir("exact_solutions/{}_{}_{}nodes_deletion".format(p1, p2, n))
	files1_ins = os.listdir("exact_solutions/{}_{}_{}nodes_insertion".format(p1, p2, n))

	#files3 = os.listdir("exact_solutions/015_015_{}nodes".format(n))
	#files3_del = os.listdir("exact_solutions/015_015_{}nodes_deletion".format(n))
	#files3_ins = os.listdir("exact_solutions/015_015_{}nodes_insertion".format(n))

	#len1 = len(files1)
	#len2 = len(files2)
	#len3 = len(files3)
	#print(len1)
	#print(len2)
	#print(len3)

	for i in range(100):
		# G1 (perturbed graph is the same for these 3 lines G1=G2=G3)
		G1, edited_G1, only_add1, only_delete1, min_edit_dist1 = LDTEditor.get_ILP_data("exact_solutions/{}_{}_{}nodes/".format(p1, p2, n) + files1[i])
		G2, edited_G2, only_add2, only_delete2, min_edit_dist2 = LDTEditor.get_ILP_data("exact_solutions/{}_{}_{}nodes_deletion/".format(p1, p2, n) + files1_del[i])	
		G3, edited_G3, only_add3, only_delete3, min_edit_dist3 = LDTEditor.get_ILP_data("exact_solutions/{}_{}_{}nodes_insertion/".format(p1, p2, n) + files1_ins[i])	

		# use all editing methods on G1 
		if not IG1:
			IG1 = InvestigateGraph(edited_G1, disturbed_G = G1)
			IG2 = copy.deepcopy(IG1)
			IG3 = copy.deepcopy(IG1)
			IG4 = copy.deepcopy(IG1)
			IG5 = copy.deepcopy(IG1)
		else:
			IG1.set_perturbed_graph(G1)
			IG1.set_G(edited_G1)

			IG2.set_perturbed_graph(G1)
			IG2.set_G(edited_G1)

			IG3.set_perturbed_graph(G1)
			IG3.set_G(edited_G1)

			IG4.set_perturbed_graph(G1)
			IG4.set_G(edited_G1)

			IG5.set_perturbed_graph(G1)
			IG5.set_G(edited_G1)

		initial_num_edges = len(G1.edges())

		exact_sol_num_edges = len(edited_G1.edges())
		exact_sol_num_edges_deletion = len(edited_G2.edges())
		exact_sol_num_edges_insertion = len(edited_G3.edges())

		cograph_edited_G, is_c1_ldt, c1_num_edges, c1_ldt_edit_dist = cograph_editing(IG1)
		triples1_edited_G, is_t1_ldt, t1_num_edges, t1_ldt_edit_dist = triples_editing(IG2, n = 100)
		triples2_edited_G, is_t2_ldt, t2_num_edges, t2_ldt_edit_dist = triples_editing(IG3, deletion = True)
		triples3_edited_G, is_t3_ldt, t3_num_edges, t3_ldt_edit_dist = triples_editing(IG4, insertion = True)
		ldt_edited_G, is_t4_ldt, t4_num_edges, t4_ldt_edit_dist = LDT_editing(IG5, deletion = True)

		init_edge_amounts.append(initial_num_edges)

		exact_sol_min_dist.append(min_edit_dist1)
		exact_sol_min_dist_del.append(min_edit_dist2)
		exact_sol_min_dist_ins.append(min_edit_dist3)

		
		exact_sol_edge_amounts.append(exact_sol_num_edges)
		exact_sol_edge_amounts_del.append(exact_sol_num_edges_deletion)
		exact_sol_edge_amounts_ins.append(exact_sol_num_edges_insertion)

		c1_graphs.append(cograph_edited_G)
		t1_graphs.append(triples1_edited_G)
		t2_graphs.append(triples2_edited_G)
		t3_graphs.append(triples3_edited_G)
		t4_graphs.append(ldt_edited_G)

		c1_is_ldt.append(is_c1_ldt)
		t1_is_ldt.append(is_t1_ldt)
		t2_is_ldt.append(is_t2_ldt)
		t3_is_ldt.append(is_t3_ldt)
		t4_is_ldt.append(is_t4_ldt)

		c1_edge_count.append(c1_num_edges)
		t1_edge_count.append(t1_num_edges)
		t2_edge_count.append(t2_num_edges)
		t3_edge_count.append(t3_num_edges)
		t4_edge_count.append(t4_num_edges)

		c1_edit_dist.append(c1_ldt_edit_dist)
		t1_edit_dist.append(t1_ldt_edit_dist)
		t2_edit_dist.append(t2_ldt_edit_dist)
		t3_edit_dist.append(t3_ldt_edit_dist)
		t4_edit_dist.append(t4_ldt_edit_dist)

	_,  cograph_freq, _			= get_freq(IG1)
	triples1_freq, _, _			= get_freq(IG2)
	triples2_freq, _, _			= get_freq(IG3)
	triples3_freq, _, _			= get_freq(IG4)
	_, 		_, 		ldt_freq 	= get_freq(IG5)

	ldt_edit_ratios = get_ratios(exact_sol_min_dist, t4_edit_dist)
	cograph_edit_ratios = get_ratios(exact_sol_min_dist, c1_edit_dist)
	t1_edit_ratios = get_ratios(exact_sol_min_dist, t1_edit_dist)
	t2_edit_ratios = get_ratios(exact_sol_min_dist_del, t2_edit_dist)
	t3_edit_ratios = get_ratios(exact_sol_min_dist_ins, t3_edit_dist)

	#bar_plots(n, 0.15, 0.15, cograph_freq, triples1_freq, triples2_freq, triples3_freq, ldt_freq)	# bar plots of frequencies of different heuristics.
	box_plots(10, 0.15, 0.15, ldt_edit_ratios, cograph_edit_ratios, t1_edit_ratios, t2_edit_ratios, t3_edit_ratios)	# box plots of ratios between edit distance and minimum edit distance
	#box_plots()	# edit distance between edited G* (G* being ldt) and original ldt graph (pre perturbed)

	# bar plot using frequencies
	# box plots with either ratio between edit_dist / min_edit_dist OR just plain edit distances. In both cases we can include horizontal lines for min_dist for ins/del and both.
	# with plain edit dist we cant rly compare to min edit dist unless we draw a line for each graph edited and even then we cant see which points correspond to which line.
#generate_solutions_fromTrees(10, 'exact_solutions/trees')

#benchmark(18, '015', '015')
#benchmark(14)
#benchmark(18)
benchmarkLDTEdits(18, 0.50, 0.15)