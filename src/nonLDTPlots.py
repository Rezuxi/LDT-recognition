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
def benchmark_fromTrees(n, p1, p2, filename='exact_solutions/trees'):
	# load species+gene trees
	name = filename + '/{}trees'.format(n)
	tree_files = []

	#probabilities = [0.15, 0.30, 0.50]
	#probs = [(0.15, 0.15), (0.3, 0.3), (0.5, 0.5), (0.15, 0.5), (0.5, 0.15)]
	#nodes = [10, 14, 18]


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

	c1_graphs, c1_edge_count, c1_is_ldt, c1_edit_dist = ([] for i in range(4))	# cograph editing
	t1_graphs, t1_edge_count, t1_is_ldt, t1_edit_dist = ([] for i in range(4))	# triples editing with both insertion/deletion
	t2_graphs, t2_edge_count, t2_is_ldt, t2_edit_dist = ([] for i in range(4))	# triples editing with deletion only
	t3_graphs, t3_edge_count, t3_is_ldt, t3_edit_dist = ([] for i in range(4))	# triples editing with insertion only
	t4_graphs, t4_edge_count, t4_is_ldt, t4_edit_dist = ([] for i in range(4))	# ldt editing

	IG1 = None
	for i in range(len(species_trees)):
		print("Tree pair {}".format(i))
		S = PhyloTree.load(name + '/{}'.format(species_trees[i]))
		TGT = PhyloTree.load(name + '/{}'.format(gene_trees[i]))
		OGT = te.observable_tree(TGT)
		ldt = ldt_graph(OGT, S)


		if not IG1:
			IG1 = InvestigateGraph(ldt)
			IG2 = copy.deepcopy(IG1)
			IG3 = copy.deepcopy(IG1)
			IG4 = copy.deepcopy(IG1)
			#IG5 = copy.deepcopy(IG1)

		IG1.perturb_graph_terminate(p1, p2)
		IG2.perturb_graph_terminate(p1, p2)
		IG3.perturb_graph_terminate(p1, p2)
		IG4.perturb_graph_terminate(p1, p2)
		#IG5.perturb_graph_terminate(p1, p2)
		


		cograph_edited_G, is_c1_ldt, c1_num_edges, c1_ldt_edit_dist = cograph_editing(IG1)
		triples1_edited_G, is_t1_ldt, t1_num_edges, t1_ldt_edit_dist = triples_editing(IG2, n = 100)
		triples2_edited_G, is_t2_ldt, t2_num_edges, t2_ldt_edit_dist = triples_editing(IG3, deletion = True)
		triples3_edited_G, is_t3_ldt, t3_num_edges, t3_ldt_edit_dist = triples_editing(IG4, insertion = True)
		#ldt_edited_G, is_t4_ldt, t4_num_edges, t4_ldt_edit_dist = LDT_editing(IG5, deletion = True)

		c1_graphs.append(cograph_edited_G)
		t1_graphs.append(triples1_edited_G)
		t2_graphs.append(triples2_edited_G)
		t3_graphs.append(triples3_edited_G)
		#t4_graphs.append(ldt_edited_G)

		c1_is_ldt.append(is_c1_ldt)
		t1_is_ldt.append(is_t1_ldt)
		t2_is_ldt.append(is_t2_ldt)
		t3_is_ldt.append(is_t3_ldt)
		#t4_is_ldt.append(is_t4_ldt)

		c1_edge_count.append(c1_num_edges)
		t1_edge_count.append(t1_num_edges)
		t2_edge_count.append(t2_num_edges)
		t3_edge_count.append(t3_num_edges)
		#t4_edge_count.append(t4_num_edges)

		c1_edit_dist.append(c1_ldt_edit_dist)
		t1_edit_dist.append(t1_ldt_edit_dist)
		t2_edit_dist.append(t2_ldt_edit_dist)
		t3_edit_dist.append(t3_ldt_edit_dist)
		#t4_edit_dist.append(t4_ldt_edit_dist)

	_,  cograph_freq, _			= get_freq(IG1) 
	triples1_freq, _, _			= get_freq(IG2)
	triples2_freq, _, _			= get_freq(IG3)
	triples3_freq, _, _			= get_freq(IG4)
	#_, 		_, 		ldt_freq 	= get_freq(IG5)

	frequencies = [cograph_freq, triples1_freq, triples2_freq, triples3_freq]	
	
	return frequencies












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


	

def benchmark(n, P1, P2):
	p1 = str(P1).replace('.', '')
	p1 = p1 + '0' if len(p1) == 2 else p1
	
	p2 = str(P2).replace('.', '')
	p2 = p2 + '0' if len(p2) == 2 else p2
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

	_,  cograph_freq, _			= get_freq(IG1) # list of 4 values
	triples1_freq, _, _			= get_freq(IG2)
	triples2_freq, _, _			= get_freq(IG3)
	triples3_freq, _, _			= get_freq(IG4)
	_, 		_, 		ldt_freq 	= get_freq(IG5)

	ldt_edit_ratios = get_ratios(exact_sol_min_dist, t4_edit_dist)
	cograph_edit_ratios = get_ratios(exact_sol_min_dist, c1_edit_dist)
	t1_edit_ratios = get_ratios(exact_sol_min_dist, t1_edit_dist)
	t2_edit_ratios = get_ratios(exact_sol_min_dist_del, t2_edit_dist)
	t3_edit_ratios = get_ratios(exact_sol_min_dist_ins, t3_edit_dist)

	frequencies = [cograph_freq, triples1_freq, triples2_freq, triples3_freq]	# 4 lists of 4 values each -> 16 values
	ratios = [cograph_edit_ratios, t1_edit_ratios, t2_edit_ratios, t3_edit_ratios] # 4 lists of ratios
	# returns the frequencies for each edit on each tick for n, p.
	return frequencies, ratios

def non_LDT_plots(ls, n_values = [10, 14, 18], probs = [(0.15, 0.15), (0.3, 0.3), (0.5, 0.5)]):
	S1 = "$G^*$ became a cograph"
	S2 = "$R_{G^*}$ became consistent"
	S3 = "$G^*$ remained properly colored"
	S4 = "$G^*$ became an LDT-graph"
	lbls = ['cograph editing', 'Triples editing ($k=100$)', 'Triples editing (deletion)', 'Triples editing (insertion)']
	X_bar = [S1, S2, S3, S4]
	X_axis = np.arange(0, 8, 2)
	f_size = 16
	w = 0.3
	rot = 30
	a = len(n_values) # rows
	b = len(probs) # cols
	#colors = ['#a83832', '#6ba832', '#32a885', '#3273a8', '#9b32a8']
	colors = ['red', 'green', '#5da2f0', '#e09710', '#cc00ff']
	#probs = [(0.15, 0.15), (0.30, 0.30), (0.50, 0.50)]

	fig, axs = plt.subplots(a, b, sharey=True, sharex=True, figsize=(24, 16))
	#handles = []
	for i in range(a): # each n
		n_data = ls[i]
		for j in range(b): # each p
			data = n_data[j]
			
			axs[i, j].bar(X_axis + w*0, data[0], label = 'cograph editing', width = w, color = 'blue')
			axs[i, j].bar(X_axis + w*1, data[1], label = 'Triples editing ($k=100$)', width = w, color = '#1a9c79')
			axs[i, j].bar(X_axis + w*2, data[2], label = 'Triples editing (deletion)', width = w, color = '#6a12b3')
			axs[i, j].bar(X_axis + w*3, data[3], label = 'Triples editing (insertion)', width = w, color = 'red')
		
			axs[i, j].set_ylim(0, 1)
			axs[i, j].set_title('$n={}$, $p={}$'.format(n_values[i], probs[j]), fontsize=18)
			if j == 0:
				axs[i, j].set_ylabel('Frequency', fontsize=16, weight='bold', labelpad=10, fontfamily='cursive')
		
			axs[i, j].set_xticks(X_axis)
			axs[i, j].set_xticklabels(X_bar, fontsize=14, weight='heavy', fontfamily='cursive', rotation = rot)
			#if len(handles) < 4:
			#	handles.append(p1)
		
	

	fig.suptitle("different heuristics applied to 100 perturbed LDT-graphs with $n$ vertices and perturbation probability $p=(p_{ins}, p_{del})$", fontsize=18)
	handles, labels = axs[0, 0].get_legend_handles_labels()
	#fig.legend([handles[0]['bars'][0], handles[1]['bars'][0], handles[2]['bars'][0], handles[3]['bars'][0]], probs, loc='upper right', prop={'size': 12}, fontsize=14,  bbox_to_anchor=(1.00, 0.885))
	fig.legend(handles, lbls, loc='upper right', prop={'size': 12}, fontsize=14,  bbox_to_anchor=(1.00, 1.00))
	#plt.tight_layout()
	fig.savefig('nonLDTplots/bar1.png', dpi=150)



def non_LDT_boxplots(ls):
	box_title = r"Ratios between $\dfrac{A_i}{X_i}$ for resulting LDT-graphs of cograph and triples editing."
	box_title = box_title + "\nThe graphs these edits are applied to have $n$ vertices and were obtained by perturbing LDT-graphs with $p=(p_{ins}, p_{del})$."
	X_bar = ['cograph editing', 'triples editing (k=100)', 'triples editing (deletion)', 'triples editing (insertion)']
	X_axis = np.arange(len(X_bar))

	w = 0.2
	#pos = [0.1, 1.3, 2.5]
	pos1 = [0, 1, 2, 3]
	pos2 = [0.2, 1.2, 2.2, 3.2]
	pos3 = [0.4, 1.4, 2.4, 3.4]


	pos = [pos1, pos2, pos3]
	colors = ['C0', 'C5', 'C2']
	probs = ['p=(0.15, 0.15)', 'p=(0.3, 0.3)', 'p=(0.5, 0.5)']
	fig, axs = plt.subplots(2, 1, sharex=True, figsize=(20, 18))
	handles = []

	for i in range(2):
		# subplots and n values
		n_data = ls[i]
		axs[i].set_title('$n={}$'.format(10+i*8), fontsize=18)

		for j in range(3):

			 
			pos_j = pos[j]
			# ratios for all edits (p_j)
			ratios = n_data[j]
			p1 = axs[i].boxplot(ratios, positions=pos_j, widths=w, sym='', patch_artist=True, boxprops=dict(facecolor=colors[j]))

			axs[i].set_xticks(X_axis)
			axs[i].set_xticklabels(X_bar, fontsize=16, weight='heavy', fontfamily='cursive', rotation=30)
			axs[i].set_ylabel('Ratio', fontsize=16, weight='bold', fontfamily='cursive')
			axs[i].tick_params(axis='y', labelsize=16)
			
			if len(handles) < 3:
				handles.append(p1)

	#print(handles)
	fig.suptitle(box_title, fontsize=20)
	fig.legend([handles[0]['boxes'][0], handles[1]['boxes'][0], handles[2]['boxes'][0]], probs, loc='upper right', prop={'size': 14}, fontsize=14)
	#plt.tight_layout()
	#plt.show()
	
	fig.savefig('nonLDTplots/box.png', dpi=100)	

def main():
	f1, r1 = benchmark(10, 0.15, 0.15)

	f2, r2 = benchmark(10, 0.3, 0.3)
	f3, r3 = benchmark(10, 0.5, 0.5)
	
	f4, r4 = benchmark(18, 0.15, 0.15)
	f5, r5 = benchmark(18, 0.3, 0.3)
	f6, r6 = benchmark(18, 0.5, 0.5)
	

	n1_freq = [f1, f2, f3]
	n2_freq = [f4, f5, f6]

	n1_ratios = [r1, r2, r3] # list of 3 lists with 4 lists of ratios each. 12 plots
	n2_ratios = [r4, r5, r6]

	all_freq = [n1_freq, n2_freq]
	all_ratios = [n1_ratios, n2_ratios]		# ls[2][3][4] -> n, p, edits

	non_LDT_plots(all_freq)
	non_LDT_boxplots(all_ratios)


# plots for larger graphs (from trees)
def main2():
	f1 = benchmark_fromTrees(40, 0.15, 0.15)
	f2 = benchmark_fromTrees(40, 0.3, 0.3)
	f3 = benchmark_fromTrees(40, 0.5, 0.5)

	#f4 = benchmark_fromTrees(50, 0.15, 0.15)
	#f5 = benchmark_fromTrees(50, 0.3, 0.3)
	#f6 = benchmark_fromTrees(50, 0.5, 0.5)

	f7 = benchmark_fromTrees(50, 0.15, 0.15)
	f8 = benchmark_fromTrees(50, 0.3, 0.3)
	f9 = benchmark_fromTrees(50, 0.5, 0.5)

	n1_freq = [f1, f2, f3]
	#n2_freq = [f4, f5, f6]
	n3_freq = [f7, f8, f9]

	all_freq = [n1_freq, n3_freq]
	non_LDT_plots(all_freq, n_values=[40, 50])

#main()
main2()