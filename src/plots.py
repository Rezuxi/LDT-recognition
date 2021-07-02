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
import numpy as np
from matplotlib.legend import _get_legend_handles_labels
np.random.seed(1)

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
	
	return triples_values, cograph_values, ldt_f #ldt_values

def get_ratios(min_dist, edit_dist):
	ratios = []
	for i in range(len(min_dist)):
		if edit_dist[i]:
			ratios.append(edit_dist[i]/min_dist[i])
	return ratios



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

	t1_graphs, t1_edge_count, t1_is_ldt, t1_edit_dist = ([] for i in range(4))	# ldt editing	
	t2_graphs, t2_edge_count, t2_is_ldt, t2_edit_dist = ([] for i in range(4))	# ldt editing (triples edit deletion)
	t3_graphs, t3_edge_count, t3_is_ldt, t3_edit_dist = ([] for i in range(4))	# ldt editing (triples edit insertion)		

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


		IG1.perturb_graph_terminate(p1, p2)
		IG2.perturb_graph_terminate(p1, p2)
		IG3.perturb_graph_terminate(p1, p2)

		

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


	_, _, ldt1_freq				= get_freq(IG1)
	_, _, ldt2_freq				= get_freq(IG2)
	_, _, ldt3_freq				= get_freq(IG3)

	frequencies = [ldt1_freq, ldt2_freq, ldt3_freq]	
	
	return frequencies





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

	freqs = [ldt1_freq, ldt2_freq, ldt3_freq]	# freqs for LDT edit i-iii
	ratios = [ldt1_edit_ratios, ldt2_edit_ratios, ldt3_edit_ratios]	# ratios for LDT edits i-iii

	return freqs, ratios
	


def barPlots(ls):
	ticks = ['LDT editing (i)', 'LDT editing (ii)', 'LDT editing (iii)']
	#x_pos = np.arange(len(ticks))	# set step to increase distance between x ticks
	x_pos = [0, 0.6, 1.2]
	w = 0.1
	#colors = ['#a83832', '#6ba832', '#32a885', '#3273a8', '#9b32a8']
	colors = ['red', 'green', '#5da2f0', '#e09710', '#cc00ff']
	probs = [(0.15, 0.15), (0.30, 0.30), (0.50, 0.50), (0.15, 0.50), (0.50, 0.15)]

	fig, axs = plt.subplots(3, 1, sharey=True, sharex=True, figsize=(16, 10))
	for i in range(3):
		m = max(ls)
		data = ls[i] # all data for n vertices
		#axs[i].title.set_text('$n={}$'.format(10+4*i))
		axs[i].set_title('$n={}$'.format(10+4*i), fontsize=18)#, horizontalalignment='left', x=0.15)
		axs[i].set_ylim(0, 1)
		axs[i].set_ylabel('Frequency', fontsize=16, weight='bold', labelpad=10, fontfamily='cursive')
		#axs[i].tick_params(axis='y', which='major', pad=10)
		for k in range(len(probs)):	
			h1 = data[k][0] # ldt edit 1
			h2 = data[k][1]
			h3 = data[k][2]
			axs[i].bar(x_pos[0] + w*k, h1, color=colors[k], width=w, label='p={}'.format(probs[k]), alpha=1.0)
			axs[i].bar(x_pos[1] + w*k, h2, color=colors[k], width=w)
			axs[i].bar(x_pos[2] + w*k, h3, color=colors[k], width=w)

			axs[i].set_xticks(x_pos)
			axs[i].set_xticklabels(ticks, fontsize=18, weight='heavy', fontfamily='cursive')

			# text above and inside bar
			# only include text inside (i.e frequency) if its above 0.1
			'''
			axs[i].text(x_pos[0] + w*k, h1+0.1, '$p={}$'.format(probs[k]), fontsize=12, color='black', horizontalalignment='center')
			axs[i].text(x_pos[1] + w*k, h1+0.1, '$p={}$'.format(probs[k]), fontsize=12, color='black', horizontalalignment='center')
			axs[i].text(x_pos[2] + w*k, h1+0.1, '$p={}$'.format(probs[k]), fontsize=12, color='black', horizontalalignment='center')
			
			if h1 > 0.1:
				axs[i].text(x_pos[0] + w*k, 0.1, '${}\%$'.format(h1*100), fontsize=12, color='black', horizontalalignment='center')
			if h2 > 0.1:
				axs[i].text(x_pos[1] + w*k, 0.1, '${}\%$'.format(h2*100), fontsize=12, color='black', horizontalalignment='center')
			if h3 > 0.1:
				axs[i].text(x_pos[2] + w*k, 0.1, '${}\%$'.format(h3*100), fontsize=12, color='black', horizontalalignment='center')
			'''
			if h1 > 0.90:
				t1 = h1 - 0.07	# inside
			else:
				t1 = h1 + 0.02	# above
			if h2 > 0.90:
				t2 = h2 - 0.07	# inside
			else:
				t2 = h2 + 0.02	# above
			if h3 > 0.90:
				t3 = h3 - 0.07	# inside
			else:
				t3 = h3 + 0.02	# above
			axs[i].text(x_pos[0] + w*k, t1, '{}%'.format(int(h1*100)), fontsize=12, color='black', horizontalalignment='center', weight='bold')
			axs[i].text(x_pos[1] + w*k, t2, '{}%'.format(int(h2*100)), fontsize=12, color='black', horizontalalignment='center', weight='bold')
			axs[i].text(x_pos[2] + w*k, t3, '{}%'.format(int(h3*100)), fontsize=12, color='black', horizontalalignment='center', weight='bold')

	fig.suptitle("Success rates of LDT edits on 100 perturbed LDT-graphs with $n$ vertices and perturbation probability $p=(p_{ins}, p_{del})$", fontsize=16)
	#plt.legend(probs, loc='upper left', bbox_to_anchor=(1.00, 1.00))
	#fig.legend(*_get_legend_handles_labels(fig.axes))
	axs[0].legend(loc='upper left', bbox_to_anchor=(1.00, 1.025), prop={'size': 12}, fontsize=14)
	plt.tight_layout()
	fig.savefig('testPlots/bar.png', dpi=100)
	#plt.show()

def boxPlots(ls, ylim = False, filename='testPlots/box{}.png', ID=0):
	box_title = r"Ratios between $\dfrac{C_i}{X_i}$ for LDT editing with different restrictions for triples editing."
	box_title = box_title + "\nThe graphs these edits are applied to have $n$ vertices and were obtained by perturbing LDT-graphs with $p=(p_{ins}, p_{del})$."
	ticks = ['LDT editing (i)', 'LDT editing (ii)', 'LDT editing (iii)']
	w = 0.2
	pos = [0.1, 1.3, 2.5]
	pos1 = [0.1, 0.3, 0.5, 0.7, 0.9]
	pos2 = [1.3, 1.5, 1.7, 1.9, 2.1]
	pos3 = [2.5, 2.7, 2.9, 3.1, 3.3]
	#pos4 = []
	#pos5 = []
	positsns = [pos1, pos2, pos3]
	colors = ['C0', 'C5', 'C2', 'C3', 'C4']
	probs = ['p=(0.15, 0.15)', 'p=(0.3, 0.3)', 'p=(0.5, 0.5)', 'p=(0.15, 0.5)', 'p=(0.5, 0.15)']
	fig, axs = plt.subplots(3, 1, sharex=True, figsize=(20, 12))
	handles = []
	for i in range(3):
		# subplots and n values
		n_data = ls[i]
		axs[i].set_title('$n={}$'.format(10+i*4), fontsize=18)

		for j in range(3):
			# triples edits (i-iii)
			edits_j = n_data[j] 
			pos_j = positsns[j]
			for k in range(5):
				# ratios for p1, p2, p3, p4, p5
				ratios = edits_j[k]
				p1 = axs[i].boxplot(ratios, positions=[pos_j[k]], widths=w, sym='', patch_artist=True, boxprops=dict(facecolor=colors[k]))
				axs[i].set_xticks(pos)
				axs[i].set_xticklabels(ticks, fontsize=18, weight='heavy', fontfamily='cursive')
				axs[i].set_ylabel('Ratio', fontsize=16, weight='bold', fontfamily='cursive')
				axs[i].tick_params(axis='y', labelsize=16)
				if len(handles) < 5:
					handles.append(p1)
		if ylim:
			axs[i].set_ylim(0.9, 5)
		else:
			axs[i].set_ylim(0.9, None)

	fig.suptitle(box_title, fontsize=16)
	fig.legend([handles[0]['boxes'][0], handles[1]['boxes'][0], handles[2]['boxes'][0], handles[3]['boxes'][0], handles[4]['boxes'][0]], probs, loc='upper right', prop={'size': 12}, fontsize=14,  bbox_to_anchor=(1.00, 0.885))
	#plt.tight_layout()
	#plt.show()
	file = filename.format(ID)
	fig.savefig(file, dpi=100)
def main():
	
	f1, r1 = benchmarkLDTEdits(10, 0.15, 0.15)
	f2, r2 = benchmarkLDTEdits(10, 0.3, 0.3)
	f3, r3 = benchmarkLDTEdits(10, 0.5, 0.5)
	f4, r4 = benchmarkLDTEdits(10, 0.15, 0.5)
	f5, r5 = benchmarkLDTEdits(10, 0.5, 0.15)
	
	f6, r6 = benchmarkLDTEdits(14, 0.15, 0.15)
	f7, r7 = benchmarkLDTEdits(14, 0.3, 0.3)
	f8, r8 = benchmarkLDTEdits(14, 0.5, 0.5)
	f9, r9 = benchmarkLDTEdits(14, 0.15, 0.5)
	f10, r10 = benchmarkLDTEdits(14, 0.5, 0.15)
	
	f11, r11 = benchmarkLDTEdits(18, 0.15, 0.15)
	f12, r12 = benchmarkLDTEdits(18, 0.3, 0.3)
	f13, r13 = benchmarkLDTEdits(18, 0.5, 0.5)
	f14, r14 = benchmarkLDTEdits(18, 0.15, 0.5)
	f15, r15 = benchmarkLDTEdits(18, 0.5, 0.15)
	
	
	ratios1 = [r1[0], r2[0], r3[0], r4[0], r5[0]]	# p1-5 of ldt1 10nodes
	ratios2 = [r1[1], r2[1], r3[1], r4[1], r5[1]]	# p1-5 of ldt2 10nodes
	ratios3 = [r1[2], r2[2], r3[2], r4[2], r5[2]]	# p1-5 of ldt3 10nodes

	ratios4 = [r6[0], r7[0], r8[0], r9[0], r10[0]]	# p1-5 of ldt1 14nodes
	ratios5 = [r6[1], r7[1], r8[1], r9[1], r10[1]]	# p1-5 of ldt2 14nodes
	ratios6 = [r6[2], r7[2], r8[2], r9[2], r10[2]]	# p1-5 of ldt3 14nodes

	ratios7 = [r11[0], r12[0], r13[0], r14[0], r15[0]]	# p1-5 of ldt1 18nodes
	ratios8 = [r11[1], r12[1], r13[1], r14[1], r15[1]]	# p1-5 of ldt2 18nodes
	ratios9 = [r11[2], r12[2], r13[2], r14[2], r15[2]]	# p1-5 of ldt3 18nodes

	n1_ratios = [ratios1, ratios2, ratios3]	# LDT 1-3 10n
	n2_ratios = [ratios4, ratios5, ratios6]	# LDT 1-3 14n
	n3_ratios = [ratios7, ratios8, ratios9]	# LDT 1-3 18n

	all_ratios = [n1_ratios, n2_ratios, n3_ratios]
	'''
	freqs1 = [f1[0], f2[0], f3[0], f4[0], f5[0]]	# frequencies for LDT 1 n10
	freqs2 = [f1[1], f2[1], f3[1], f4[1], f5[1]]	# frequencies for LDT 2 n10
	freqs3 = [f1[2], f2[2], f3[2], f4[2], f5[2]]	# frequencies for LDT 3 n10

	freqs4 = [f6[0], f7[0], f8[0], f9[0], f10[0]]	# frequencies for LDT 1 n14
	freqs5 = [f6[1], f7[1], f8[1], f9[1], f10[1]]	# frequencies for LDT 2 n14
	freqs6 = [f6[2], f7[2], f8[2], f9[2], f10[2]]	# frequencies for LDT 3 n14

	freqs7 = [f11[0], f12[0], f13[0], f14[0], f15[0]]	# frequencies for LDT 1 n18
	freqs8 = [f11[1], f12[1], f13[1], f14[1], f15[1]]	# frequencies for LDT 2 n18
	freqs9 = [f11[2], f12[2], f13[2], f14[2], f15[2]]	# frequencies for LDT 3 n18

	n1_freqs = [freqs1, freqs2, freqs3]	# LDT 1, 2, 3 n = 10
	n2_freqs = [freqs4, freqs5, freqs6]	# LDT 1, 2, 3 n = 14
	n3_freqs = [freqs7, freqs8, freqs9]	# LDT 1, 2, 3 n = 18

	all_freqs = [n1_freqs, n2_freqs, n3_freqs] # n = 10, 14, 18
	'''
	# f1 has 3 values. 1 for each edit
	freqs1 = [f1, f2, f3, f4, f5]
	freqs2 = [f5, f6, f7, f9, f10]
	freqs3 = [f11, f12, f13, f14, f15]
	all_freqs = [freqs1, freqs2, freqs3]
	'''
		3*5*3 = 45 bars
		for i in range(3):
			data = ls[i] 

			for k in range(5)
				h1[k][0]
				plot 3 bars 
			


	'''
	# DUMMY DATA
	#r1 = np.random.randn(40) # ratios for ldt i p = 1
	#r2 = np.random.randn(40) # ratios for ldt i p = 2
	#r3 = np.random.randn(40) # ratios for ldt i p = 3
	#r4 = np.random.randn(40) # ratios for ldt i p = 4
	#r5 = np.random.randn(40) # ratios for ldt i p = 5. use the same dummy data for ldt ii and iii
	
	#f1 = 0.1
	#f2 = 0.5
	#f3 = 0.7
	#f4 = 0.3
	#f5 = 1.0
	# should be a list of 3 lists L_i where L_i has 5 lists P_j where P_j has
	#ldt1_ratios = [p1, p2, p3, p4, p5]
	#ldt_ratios10 = [ldt1_ratios, ldt2_ratios, ldt3_ratios]
	#ldt_ratios14 = [ldt1_ratios, ldt2_ratios, ldt3_ratios]
	#ldt_ratios18 = [ldt1_ratios, ldt2_ratios, ldt3_ratios]

	#ratios = [r1, r2, r3, r4, r5]  		  # p1, p2, p3, p4, p5
	#ldt_ratios = [ratios, ratios, ratios] # ldt editing 1, 2, 3
	#all_ratios = [ldt_ratios, ldt_ratios, ldt_ratios] # n 10, 14, 18

	#freqs = [f1, f2, f3, f4, f5]


	barPlots(all_freqs)
	boxPlots(all_ratios)
	boxPlots(all_ratios, ylim=True, ID = 1)




def main2():
	f1 = benchmark_fromTrees(10, 0.15, 0.15)
	f2 = benchmark_fromTrees(10, 0.3, 0.3)
	f3 = benchmark_fromTrees(10, 0.5, 0.5)
	f4 = benchmark_fromTrees(10, 0.15, 0.5)
	f5 = benchmark_fromTrees(10, 0.5, 0.15)

	f6 = benchmark_fromTrees(10, 0.15, 0.15)
	f7 = benchmark_fromTrees(10, 0.3, 0.3)
	f8 = benchmark_fromTrees(10, 0.5, 0.5)
	f9 = benchmark_fromTrees(10, 0.15, 0.5)
	f10 = benchmark_fromTrees(10, 0.5, 0.15)

	f11 = benchmark_fromTrees(10, 0.15, 0.15)
	f12 = benchmark_fromTrees(10, 0.3, 0.3)
	f13 = benchmark_fromTrees(10, 0.5, 0.5)
	f14 = benchmark_fromTrees(10, 0.15, 0.5)
	f15 = benchmark_fromTrees(10, 0.5, 0.15)

	freqs1 = [f1, f2, f3, f4, f5]
	freqs2 = [f5, f6, f7, f9, f10]
	freqs3 = [f11, f12, f13, f14, f15]
	all_freqs = [freqs1, freqs2, freqs3]
	
	barPlots(all_freqs)

#main()