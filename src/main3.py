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

k = count_files(10, 'exact_solutions/015_030')
print(k)

#generate_solutions_fromTrees(10, 'exact_solutions/trees')

