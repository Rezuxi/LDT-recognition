import networkx as nx
import asymmetree.tools as tools
import asymmetree.cograph as cg
from asymmetree.tools.GraphTools import disturb_graph, symmetric_diff
import random
import copy
import itertools


#####################################################################################
#																					#
#						Induced paths on 3 vertices (P3)							#
#																					#
#####################################################################################

# TODO: Try breaking all triangles by cutting 2 of their edges so that they dont form new P3s.
#		also try only breaking triangles that have an edge from bad triples in them.

# total runtime O(|E|*2(|V|-1)) -> 2(|E||V| - |E|) -> O(|E||V|)
def get_species_triples(G, colored_G = None):
	'''
		Returns the set of species triples and the leaves of the species tree (colors in the graph).
	'''
	species_leaves = set()
	# list of all species triples
	species_triples = []
	# list of induced complete subgraphs on 3 vertices with pairwise distinct colors
	K3s = set()
	triangles = []
	# assign each gene triple (p3 with distinct colors) to a species triple so that we know which gene triples we need to remove in order to remove
	# a species triple.
	triples_dict = {}

	if colored_G:
		nodes = colored_G.nodes()
	else:
		nodes = G.nodes()

	for u, v in G.edges():	# O(|E|)
		u_adj = G.adj[u]
		v_adj = G.adj[v]
		A = nodes[v]['color']
		C = nodes[u]['color']

		for w in u_adj:				# O(|V|-1) (worst case scenario)
			B = nodes[w]['color']
			if A == B or A == C or B == C or w == v:
				# the 3 colors are not distinct -> not a species triple
				# if w = v then we're only looking at the vertices u, v.
				continue
			if w in v_adj :
				# check if this forms a triangle (u, v), (v, w), (w, u) with distinct colors (in case we later on want to remove one of these edges resulting in a new triples)
				uvw = sorted([u, v, w])
				# NOTE: THERE ARE DUPLICATES HERE
				triangles.append((*uvw,))
				K3s.add((*uvw,))	# this does not include duplicates
				continue
			else:
				# check if the species triple exists
				# sort the colors a,b in the triple ab|c
				AB = sorted([A, B])
				triple = (*AB, C)
				if triple in species_triples:
					# if it's in species_triples, it's been added to triples_dict aswell
					# add the gene triple to the triples dict. (v, w, u)
					vw = sorted([v, w])
					triples_dict[(*AB, C)].append((*vw, u))
				else:
					species_triples.append(triple)
					triples_dict[(*AB, C)] = []
					species_leaves.add(A)
					species_leaves.add(B)
					species_leaves.add(C)
			
		for w in v_adj:				# O(|V|-1) (worst case scenario)
			B = nodes[w]['color']
			if A == B or A == C or B == C or w == u:
				continue
			if w in u_adj:
				uvw = sorted([u, v, w])
				triangles.append((*uvw,))
				K3s.add((*uvw,))
			else:
				BC = sorted([B, C])
				triple = (*BC, A)
				if triple in species_triples:
					uw = sorted([u, w])
					triples_dict[(*BC, A)].append((*uw, v))
				else:
					species_triples.append(triple)
					triples_dict[(*BC, A)] = []
					species_leaves.add(A)
					species_leaves.add(B)
					species_leaves.add(C)
	return species_triples, species_leaves, K3s, triples_dict








# species triples introduced by adding an edge
# O(2(|V|-1)) -> O(|V|)
def get_triples_of_P3(G, u, v):
	species_triples = []
	u_adj = G.adj[u]
	v_adj = G.adj[v]
	A = G.nodes[v]['color']
	C = G.nodes[u]['color']

	for w in u_adj:
		B = G.nodes[w]['color']
		if A == B or A == C or B == C or w == v:
			continue
		if w in v_adj:
			continue
		AB = sorted([A, B])
		triple = (*AB, C)
		if triple in species_triples:
			continue
		species_triples.append(triple)
	for w in v_adj:
		B = G.nodes[w]['color']
		if A == B or A == C or B == C or w == u:
			continue
		if w in u_adj:
			continue
		BC = sorted([B, C])
		triple = (*BC, A)
		if triple in species_triples:
			continue
		species_triples.append(triple)
	return species_triples

# O(|E| + |R| * (|A| * (|B| + |K|))) = O(|E| + |R||A||B| + |R||A||K|) 
def weight_edges(G, triples_dict, triangles, cut_list):
	insert_edges = {}	# weight of insert edges
	for e in G.edges():			# O(|E|)
		# init attributes
		G.edges[e]['weight'] = 0
		G.edges[e]['unknowns'] = 0
		G.edges[e]['bad'] = 0

	# go through all species triples (A, B, C) and check if (A, B) is in cut_edge, then we want to remove this triple (A, B, C)
	# This is done by destroying all P3s forming this species triple.
	for A, B, C in triples_dict:						# O(|R|), where R is the set of species triples
		# if this is true, we need to break the species triple (A, B, C)
		if (A, B) in cut_list or (B, A) in cut_list:
			AB = sorted([A, B])
			AC = sorted([A, C])
			BC = sorted([B, C])
			species_triple = (*AB, C)								# AB|C

			# COUNT OCCURRENCES OF EDGES IN P3s THAT FORM SPECIES TRIPLES IN CUT_LIST
			for a, b, c in triples_dict[species_triple]:			# ab|c. this triple is sorted so that a < b.	O(|A|), |A| <= total amount of colored P3s of G
				ac = sorted([a, c]) 	# delete edge
				bc = sorted([b, c])		# delete edge
				ab = sorted([a, b])		# insert edge
				if not (*ab,) in insert_edges:
					insert_edges[(*ab,)] = {}
					insert_edges[(*ab,)]['weight'] = 1
					insert_edges[(*ab,)]['unknowns'] = 0
					insert_edges[(*ab,)]['bad'] = 0
				else:
					insert_edges[(*ab,)]['weight'] += 1
				#if G.has_edge(*ac,) and G.has_edge(*bc,)
				G.edges[(*ac,)]['weight'] += 1
				G.edges[(*bc,)]['weight'] += 1
				

				a_adj = G.adj[a]
				b_adj = G.adj[b]

				###################################################################
				#						WEIGHT INSERT EDGE 						  #
				###################################################################
				# list of species triples being formed by adding the edge (a, b)
				potential_species_triples = get_triples_of_P3(G, a, b)

				for X, Y, Z in potential_species_triples:		# O(|B|), where B is the set of new P3s introduced as a resutl of inserting the edge (a, b)
					if (X, Y) in cut_list or (Y, X) in cut_list:
						insert_edges[(*ab,)]['bad'] += 1
					if (X, Y, Z) not in triples_dict:
						insert_edges[(*ab,)]['unknowns'] += 1
					# dont think checking that the colors are different is necessary since it's a triple (3 distinct colors).
				###################################################################
				#						WEIGHT DELETE EDGES 					  #
				###################################################################
				for triangle in triangles:			# O(|K3s|) = O(|K|)
					if c in triangle:
						# edge (a, c)
						if a in triangle:
							w = (set(triangle) ^ set(ac)).pop()
	
							W = G.nodes[w]['color']
							triple = (*AC, W)
							if (A, C) in cut_list or (C, A) in cut_list:
								G.edges[ac]['bad'] += 1
							if triple not in triples_dict:
								G.edges[ac]['unknowns'] += 1
						# edge (b, c)
						elif b in triangle:
							w = (set(triangle) ^ set(bc)).pop()
							W = G.nodes[w]['color']
							triple = (*BC, W)
							if (B, C) in cut_list or (C, B) in cut_list:
								G.edges[bc]['bad'] += 1
							if triple not in triples_dict:
								G.edges[bc]['unknowns'] += 1
	return insert_edges


# O(|E||A|(|B| + |K|))
def weight_edges2(G, triples_dict, K3s, cut_list):
	K3_edge_weights = weight_K3s(K3s)

	insert_edges = {}	# weight of insert edges
	for e in G.edges():		# O(|E|)
		# init attributes
		G.edges[e]['weight'] = 0
		G.edges[e]['unknowns'] = 0
		G.edges[e]['bad'] = 0


	for A, B, C in triples_dict:	# O(|R|)

		if (A, B) in cut_list or (B, A) in cut_list:
			AB = sorted([A, B])
			AC = sorted([A, C])
			BC = sorted([B, C])
			species_triple = (*AB, C)								# AB|C


			for a, b, c in triples_dict[species_triple]:			# ab|c. this triple is sorted so that a < b.	# O(|A|)
				ac = sorted([a, c]) 	# delete edge
				bc = sorted([b, c])		# delete edge
				ab = sorted([a, b])		# insert edge
				if not (*ab,) in insert_edges:
					insert_edges[(*ab,)] = {}
					insert_edges[(*ab,)]['weight'] = 1
					insert_edges[(*ab,)]['unknowns'] = 0
					insert_edges[(*ab,)]['bad'] = 0
				else:
					insert_edges[(*ab,)]['weight'] += 1
				G.edges[(*ac,)]['weight'] += 1
				G.edges[(*bc,)]['weight'] += 1
				

				a_adj = G.adj[a]
				b_adj = G.adj[b]

				###################################################################
				#						WEIGHT INSERT EDGE 						  #
				###################################################################
				# list of species triples being formed by adding the edge (a, b)
				potential_species_triples = get_triples_of_P3(G, a, b)

				for X, Y, Z in potential_species_triples:		# O(|B|)
					if (X, Y) in cut_list or (Y, X) in cut_list:
						insert_edges[(*ab,)]['bad'] += 1
					if (X, Y, Z) not in triples_dict:
						insert_edges[(*ab,)]['unknowns'] += 1
					# dont think checking that the colors are different is necessary since it's a triple (3 distinct colors).
				###################################################################
				#						WEIGHT DELETE EDGES 					  #
				###################################################################
				
				for triangle in triangles:		# O(|K|)
					if c in triangle:
						# edge (a, c)
						if a in triangle:
							w = (set(triangle) ^ set(ac)).pop()
	
							W = G.nodes[w]['color']
							triple = (*AC, W)
							if (A, C) in cut_list or (C, A) in cut_list:
								G.edges[ac]['bad'] += 1
							if triple not in triples_dict:
								G.edges[ac]['unknowns'] += 1
						# edge (b, c)
						elif b in triangle:
							w = (set(triangle) ^ set(bc)).pop()
							W = G.nodes[w]['color']
							triple = (*BC, W)
							if (B, C) in cut_list or (C, B) in cut_list:
								G.edges[bc]['bad'] += 1
							if triple not in triples_dict:
								G.edges[bc]['unknowns'] += 1
	return insert_edges




#####################################################################################
#																					#
#									Graph Editing Tools 							#
#																					#
#####################################################################################

# O(|V||E| + BUILD)
def is_compatible(G, colored_G = None):
	if colored_G:
		triples, species_leaves, _, _ = get_species_triples(G, colored_G=colored_G)
	else:
		triples, species_leaves, _, _ = get_species_triples(G)

	# if triples , leaves are empty then it's consistent
	if triples == []:
		#print("The set of species triples was empty!")
		return True
	
	B = tools.Build(triples, species_leaves)
	tree_triples = B.build_tree()

	return True if tree_triples else False
	
#  time complexity of cograph alg
def is_cograph(G):
	cotree = cg.Cotree.cotree(G)
	if cotree:
		return True
	return False

# O(|E|)
def is_properly_colored(G, colored_G = None):
	'''
		Checks if a graph is properly colored. (no adjacent vertices have the same color)
	'''
	if colored_G:
		nodes = colored_G.nodes()
	else:
		nodes = G.nodes()

	for e in G.edges():
		if nodes[e[0]]['color'] == nodes[e[1]]['color']:
			return False
	return True

# O(|V|)
def color_graph(colored_G, G):
	for n in colored_G.nodes():
		G.nodes[n]['color'] = colored_G.nodes[n]['color']

def make_properly_colored(G):
	nodes = G.nodes()
	copy_G = G.copy()
	for e in G.edges():
		if nodes[e[0]]['color'] == nodes[e[1]]['color']:
			copy_G.remove_edge(e[0], e[1])
	return copy_G

class InvestigateGraph:

	def __init__(self, G, disturbed_G = None):
		'''
			G is an LDT graph
		'''

		self._G = G
		self._G_perturbed = disturbed_G

		self._is_cograph = True
		self._is_compatible = True

		if disturbed_G:
			self._is_cograph = is_cograph(disturbed_G)
			self._is_compatible = is_compatible(disturbed_G)

		self._triplesEdit_to_LDT = False
		self._cographEdit_to_LDT = False

		'''
			Count the P3s and regions of the graphs that are/become LDT-graphs and compare them to
			non LDT-graphs. See if there is anything standing out.
		'''

		# count how many times the perturbed graph remains properly colored after cograph editing and also for triples editing when adding edges
		self._count_cographEdit_remain_properly_colored = 0
		self._count_triplesEdit_remain_properly_colored = 0

		# count how many times any of the edits results in an LDT-graph.
		self._count_cographEdit_to_LDT = 0
		self._count_triplesEdit_to_LDT = 0
		
		# count how many times a graph went from being consistent to inconsistent by doing cograph editing
		# and the same for triple editing		
		self._count_triplesEdit_broke_cograph = 0
		self._count_cographEdit_broke_consistency = 0

		# count how many times a graph went from broken consistency to fixed by doing cograph editing
		# and the same for triple editing
		self._count_triplesEdit_fixed_cograph = 0
		self._count_cographEdit_fixed_consistency = 0

		# count how many times a graph went from broken consistency to fixed by doing cograph editing
		# and the same for triple editing
		self._count_triplesEdit_remained_cograph = 0
		self._count_cographEdit_remained_consistent = 0

		self._count_not_ldt = 0
		self._count_ldtEdit_success = 0

		# count how many times the disturbed graph is not a cograph and not consistent
		# so that we can compare the success rate of each heuristic.
		self._count_dG_not_cograph = 0
		self._count_dG_not_consistent = 0
		
		# count how many times the disturbed graph remains a cograph or consistent so that 
		# we can compare the frequency at which each heuristic breaks the other property of 
		# LDT-graphs. i.e. cograph editing breaking consistency and vice versa.
		self._count_dG_cograph = 0
		self._count_dG_consistent = 0

		self._count_dG_notCograph_notConsistent = 0
		self._count_dG_cograph_notConsistent = 0
		self._count_dG_notCograph_consistent = 0

		self._count_cographEdit_success = 0
		self._count_triplesEdit_success = 0

		self._amount_of_P3 = 0



	
	def perturb_graph(self, i_rate = None, d_rate = None):
		'''
			perturbs a graph until it is not an LDT-graph
			by default random values for deletion/insertion rate
		'''

		# if still ldt graph after n loops, perhaps print some data about the graph
		while self._is_cograph and self._is_compatible:
			# randomize probabilities again
			if i_rate == None:
				i_rate = round(random.random(), 1)
			if d_rate == None:
				d_rate = round(random.random(), 1)
			if i_rate == 0.0 and d_rate == 0.0:
				i_rate = round(random.uniform(0.1, 1.0), 1)
				d_rate = round(random.uniform(0.1, 1.0), 1)
			self._G_perturbed = disturb_graph(self._G, insertion_prob=i_rate, deletion_prob=d_rate)
			if is_cograph(self._G_perturbed):
				self._is_cograph = True
			else:
				self._is_cograph = False
			if is_compatible(self._G_perturbed):
				self._is_compatible = True
			else:
				self._is_compatible = False
			#print("adding noise again!")


	def perturb_graph_terminate(self, i_rate = None, d_rate = None, max_attempts = 100):
		"""
			Returns true if the graph was perturbed to the point where it is no longer an LDT-graph
			after max_attempts it returns false if it is still LDT
		"""
		for i in range(max_attempts):
			if i_rate == None:
				i_rate = round(random.random(), 1)
			if d_rate == None:
				d_rate = round(random.random(), 1)
			if i_rate == 0.0 and d_rate == 0.0:
				i_rate = round(random.uniform(0.1, 1.0), 1)
				d_rate = round(random.uniform(0.1, 1.0), 1)
			self._G_perturbed = disturb_graph(self._G, insertion_prob=i_rate, deletion_prob=d_rate)
			if is_cograph(self._G_perturbed):
				self._is_cograph = True
			else:
				self._is_cograph = False
			if is_compatible(self._G_perturbed):
				self._is_compatible = True
			else:
				self._is_compatible = False

			if not (self._is_compatible and self._is_cograph):
				return True
		return False	


	def set_perturbed_graph(self, G):
		isCograph = is_cograph(G)
		isCompatible = is_compatible(G)
		self._G_perturbed = G

		self._is_cograph = isCograph
		self._is_compatible = isCompatible

	def set_G(self, G):
		self._G = G

	#########################################################################################################
	#																										#
	#										EDITING HEURISTICS												#
	#																										#
	#########################################################################################################


	def triples_editing(self, mincut=True, weighted_mincut=True, n = 1, deletion = False, insertion = False):
		'''
			do triples editing up to n times or until the edited graphs set of triples becomes consistent
			If restricted to insertion or deletion, then we repeat the editing until consistency, because it will terminate.
			otherwise it's limited to n iterations
		'''
		copy_G = self._G_perturbed.copy()
		v = len(copy_G.nodes())
		if insertion ^ deletion:
			n = int((v * (v-1)) / 2) 			# number of edges in a complete graph. R_G will become consistent before this.
		for i in range(n):
			triples, species_leaves, triangles, triples_dict = get_species_triples(copy_G, self._G)		

			if len(triples) == 0:
				return copy_G, None

			B = tools.Build(triples, species_leaves, mincut=mincut, weighted_mincut=weighted_mincut)	
			tree_triples = B.build_tree()															
			
			if len(B.cut_list) > 0:
				insert_edges = weight_edges(copy_G, triples_dict, triangles, B.cut_list)		
				for X, Y, Z in triples_dict:										
					if (X, Y) in B.cut_list or (Y, X) in B.cut_list:					
						for a, b, c in triples_dict[(X, Y, Z)]:								

							if not (copy_G.has_edge(a, c) and copy_G.has_edge(b, c)) or copy_G.has_edge(a, b):
								continue
							
							if insertion:
								ins_edge = (a, b)	# insert edge
								copy_G.add_edge(*ins_edge)
								continue
							
								
							best_edge, to_delete, best_del_edge = self.edge_to_edit(a, b, c, insert_edges, copy_G)
							
							if deletion:
								copy_G.remove_edge(*best_del_edge)
								continue


							if to_delete:
								copy_G.remove_edge(*best_edge)
								#broken = True
							else:
								copy_G.add_edge(*best_edge)
								#broken = True
				if is_compatible(copy_G):
					return copy_G, tree_triples
				else:
					copy_G = copy_G.copy()
			else:
				break
		return copy_G, tree_triples


	def LDT_editing1(self, n = 1, d = False, i = False):
		'''	
			1. Triples editing -> cograph editing -> color editing (if not properly colored, makes it so)
			2. check if consistent and cograph. If yes -> return LDT-graph else return None
		'''
		triples_edited_G, _ = self.triples_editing(n = n, deletion = d, insertion = i)
		is_properly_colored = True
		if (d ^ i):
			# only insert or delete so we are sure to make G consistent.
			isConsistent = True
		else:
			isConsistent = is_compatible(triples_edited_G)
		isCograph = is_cograph(triples_edited_G)

		if not isCograph:
			cograph_edited_G = self.cograph_editing(G = triples_edited_G)
		else:
			#print("Triples editing -> LDT")
			return triples_edited_G

		color_graph(self._G, cograph_edited_G)
		properClrd_cograph = make_properly_colored(cograph_edited_G)
		isCograph = is_cograph(properClrd_cograph)
		isConsistent = is_compatible(properClrd_cograph)
		if isConsistent and isCograph:
			return properClrd_cograph
		return None


	# choose between 3 edges per P3 based on prio noted above
	def edge_to_edit(self, a, b, c, ie, G):
		e1 = (a, c)	# delete edge
		e2 = (b, c)	# delete edge
		e3 = (a, b)	# insert edge

		best_edge = None
		to_delete = True

		# find best delete edge
		if G.edges[e1]['bad'] < G.edges[e2]['bad']:
			best_edge = e1
		elif G.edges[e1]['bad'] == G.edges[e2]['bad']:
			if G.edges[e1]['unknowns'] < G.edges[e2]['unknowns']:
				best_edge = e1
			elif G.edges[e1]['unknowns'] == G.edges[e2]['unknowns']:
				if G.edges[e1]['weight'] >= G.edges[e2]['weight']:
					best_edge = e1
				else:
					best_edge = e2
			else:
				best_edge = e2
		else:
			best_edge = e2
		best_del_edge = best_edge
		# compare best delete edge to insert edge
		if ie[e3]['bad'] < G.edges[best_edge]['bad']:
			best_edge = e3
			to_delete = False
		elif ie[e3]['bad'] == G.edges[best_edge]['bad']:
			if ie[e3]['unknowns'] < G.edges[best_edge]['unknowns']:
				best_edge = e3
				to_delete = False
			elif ie[e3]['unknowns'] == G.edges[best_edge]['unknowns']:	
				if ie[e3]['weight'] >= G.edges[best_edge]['weight']:
					best_edge = e3
					to_delete = False
		return best_edge, to_delete, best_del_edge


	def cograph_editing(self, G=None):
		return cg.edit_to_cograph(G) if G else cg.edit_to_cograph(self._G_perturbed)
		

	def color_editing(self, G=None):
		'''
			removes edges between any x, y such that x and y have the same color.
		'''
		if G:
			nodes = G.nodes()
		else:
			nodes = self._G.nodes()
		copy_G = self._G.copy()
		for e in G.edges():
			if nodes[e[0]]['color'] == nodes[e[1]]['color']:
				copy_G.remove_edge(e)
		return copy_G



	#########################################################################################################
	#																										#
	#											PRINT DATA													#
	#																										#
	#########################################################################################################


	def print_symmetric_diff(self, G):
		# edit distance between edited graph and perturbed graph
		edit_dist = symmetric_diff(self._G_perturbed, G)
		print("The edit distance is: {}".format(edit_dist))
		

	def print_perturbed_G_data(self):
		print("----------------------------------------------------------------NON LDT-GRAPH----------------------------------------------------------------")
		print("The amount of nodes and edges in the perturbed graph is: {}, {}".format(len(self._G_perturbed.nodes()), len(self._G_perturbed.edges())))
		print("The density of the perturbed graph is: {}".format(nx.density(self._G_perturbed)))

		print("\t\t-------------------------------P3 related-------------------------------")
		print("Amount of P3s relative to the amount of nodes per region:")
		print("Amount of P3s relative to the amount of edges per region:")

	def print_data(self):

		# check that the variable being divided by is > 0
		'''
			Need to count the number of times the perturbed graph is neither a cograph nor "consistent" to properly count the frequency at which they are both "fixed" by each edit heuristic
			also how many times the perturbed graph is a cograph but not "consistent" and vice versa.
		'''
				
		if self._count_dG_not_cograph > 0:
			cograph_f1 = self._count_cographEdit_success / self._count_dG_not_cograph						# frequency of cograph editing turning the graph into a cograph
			cograph_f4 = self._count_cographEdit_remain_properly_colored / self._count_dG_not_cograph		# frequency of the graph remaining properly colored after cograph editing
			cograph_f5 = self._count_cographEdit_to_LDT / self._count_dG_not_cograph						# frequency of cograph editing turning the graph into an LDT-graph
		else:
			cograph_f1 = -1
			cograph_f4 = -1
			cograph_f5 = -1

		# frequency of cograph editing making the set of triples consistent. that is, the set of triples goes from inconsistent to consistent as a result of cograph editing.
		if self._count_dG_notCograph_notConsistent > 0:
			cograph_f2 = self._count_cographEdit_fixed_consistency / self._count_dG_notCograph_notConsistent
		else:
			cograph_f2 = -1

		# frequency of the set of triples remaining consistent after cograph editing
		if self._count_dG_notCograph_consistent > 0:
			cograph_f3 = self._count_cographEdit_remained_consistent / self._count_dG_notCograph_consistent
		else:
			cograph_f3 = -1



		if self._count_dG_not_consistent > 0:
			triples_f1 = self._count_triplesEdit_success / self._count_dG_not_consistent 						# frequency of triples editing making the set of triples consistent.
			triples_f4 = self._count_triplesEdit_remain_properly_colored / self._count_dG_not_consistent 		# frequency of the graph remaining properly colored after triples editing
			triples_f5 = self._count_triplesEdit_to_LDT / self._count_dG_not_consistent  						# frequency of triples editing turning the graph into an LDT-graph.
		else:
			triples_f1 = -1
			triples_f4 = -1
			triples_f5 = -1

		# frequency of triples editing turning the graph into a cograph.
		if self._count_dG_notCograph_notConsistent > 0:
			triples_f2 = self._count_triplesEdit_fixed_cograph / self._count_dG_notCograph_notConsistent
		else:
			triples_f2 = -1

		# frequency of graph remaining a cograph after triples ediiting.
		if self._count_dG_cograph_notConsistent > 0:
			triples_f3 = self._count_triplesEdit_remained_cograph / self._count_dG_cograph_notConsistent
		else:
			triples_f3 = -1



		print("\n\t\t------------------------------------Cograph editing data------------------------------------")

		print("\nFrequency of cograph editing turning the graph into a cograph: {}".format(cograph_f1))

		print("\nFrequency of cograph editing making the set of triples compatible: {}".format(cograph_f2))

		#print("\nFrequency of the set of triples remaining consistent after cograph editing: {}".format(cograph_f3))

		print("\nFrequency of the graph remaining properly colored after cograph editing: {}".format(cograph_f4))

		print("\nFrequency of cograph editing turning the graph into an LDT-graph: {}".format(cograph_f5))


		print("\n\t\t------------------------------------Triples editing data------------------------------------")

		print("\nFrequency of triples editing making the set of triples compatible: {}".format(triples_f1))

		print("\nFrequency of triples editing turning the graph into a cograph: {}".format(triples_f2))

		#print("\nFrequency of the graph remaining a cograph after triples editing: {}".format(triples_f3))

		print("\nFrequency of the graph remaining properly colored after triples editing: {}".format(triples_f4))

		print("\nFrequency of triples editing turning the graph into an LDT-graph: {}.".format(triples_f5))
		

		#print("\nOf 100 perturbed graphs, the amount of times it wasn't a cograph is: {}".format(self._count_dG_not_cograph))
		#print("Of 100 perturbed graphs, the amount of times its set of triples wasn't compatible is: {}".format(self._count_dG_not_consistent))


if __name__ == "__main__":
	pass