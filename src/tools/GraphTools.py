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

# K3s are given as sorted (a, b, c) where there is an edge between all vertices
# destroy each K3 by removing 2 edges. This method removes arbitrary edges from each K3

# O(|K|)
def destroy_K3s(G, K3s):
	#print("K3s list: {}".format(K3s))
	for a, b, c in K3s:
		#print("Destroying K3: {}".format((a, b, c)))
		# destroy the edges (a, b) and (b, c)
		if G.has_edge(a, b):
			G.remove_edge(a, b)
		if G.has_edge(b, c):
			G.remove_edge(b, c)





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
# O(|K|)
def weight_K3s(K3s):
	# weight each edge of the K3s based on the number of occurrences in K3s.
	edge_weights = {}
	for a, b, c in K3s:
		e1 = (a, b)
		e2 = (a, c)
		e3 = (b, c)
		if not e1 in edge_weights:
			edge_weights[e1] = 1	
		else:
			edge_weights[e1] += 1
		if not e2 in edge_weights:
			edge_weights[e2] = 1	
		else:
			edge_weights[e2] += 1
		if not e3 in edge_weights:
			edge_weights[e3] = 1	
		else:
			edge_weights[e3] += 1
		
	return edge_weights

# O(|K| * (|K| - 1)) = O(|K|^2)
# part of triplesEditing3
def destroy_K3s_overlapping_edge(G, e, K3s, K3_edge_weights):
	# O(|K|)
	for a, b, c in K3s:
		e1 = None 	# first edge to be removed is e (if its in a K3)
		e2 = None 	# second edge to be removed is the one with least weight -> recursively cut edges until the second edge has weight 0 or K3s is empty.
		if e == (a, b):
			e1 = e
			if K3_edge_weights[(a, c)] <= K3_edge_weights[(b, c)]:
				e2 = (a, c)
			else:
				e2 = (b, c)
		elif e == (a, c):
			e1 = e
			if K3_edge_weights[(a, b)] <= K3_edge_weights[(b, c)]:
				e2 = (a, b)
			else:
				e2 = (b, c)
		elif e == (b, c):
			e1 = e
			if K3_edge_weights[(a, b)] <= K3_edge_weights[(a, c)]:
				e2 = (a, b)
			else:
				e2 = (a, c)
		# if e1 == None, then this edge did not overlap a K3
		if e1:
			if G.has_edge(*e1):
				G.remove_edge(*e1)

			if G.has_edge(*e2):
				if K3_edge_weights[e2] > 1:
					destroy_K3s_overlapping_edge(G, e2, K3s, K3_edge_weights)		# O(|K| - 1) (worst case) ?. This runs at most |K| - 1 times since we can at most remove all K3s excluding this.
		# TODO: update list of K3s as we remove edges



# returns gene triples
def find_all_P3(G, get_triples=False, colored_G=None):
	'''
		Finds all connected triples that don't form a triangle
		if triples = True (the graph needs to be colored), we also check that the colors are distinct for each node in a triple (a, b, c)

		returns the P3s as a list of lists
		returns the triples as a list of tuples
	'''
	leaves = set()
	triples = []

	if colored_G:
		nodes = colored_G.nodes()	
	else:
		nodes = G.nodes()
	
	for e in G.edges():
		u_ID = e[0]
		v_ID = e[1]
		
		u_adj = G.adj[u_ID]
		v_adj = G.adj[v_ID]


		for w_ID in u_adj:
			# if they're the same node, skip
			if w_ID == v_ID:
				#print("{} and {} is the same node".format(w_ID, v_ID))
				continue
			# if they form a triangle, skip
			if w_ID in v_adj:
				#print("Triangle detected for the nodes: {}".format((v_ID, u_ID, w_ID)))
				continue
			# if we want a set of triples
			if get_triples == True:
				# if the nodes' colors are not pairwise distinct, skip
				# TODO: take into account that the graph might not be properly colored, so check that all colors of u, v, w are distinct. (after cograph editing, the graph might not be properly colored)
				if nodes[w_ID]['color'] == nodes[v_ID]['color'] or nodes[u_ID]['color'] == nodes[v_ID]['color'] or nodes[u_ID]['color'] == nodes[w_ID]['color']:
					#print("Not pairwise distinct colors for: {}".format((v_ID, w_ID, u_ID)))
					continue
				else:
					# if the triples has been counted already
					if (w_ID, v_ID, u_ID) in triples or (v_ID, w_ID, u_ID) in triples:
						#print("This triple is already included: {}".format((v_ID, w_ID, u_ID)))
						continue
					else:
						#print("Adding the triple: {}".format((w_ID, v_ID, u_ID)))
						triples.append((w_ID, v_ID, u_ID))
						leaves.add(u_ID)
						leaves.add(v_ID)
						leaves.add(w_ID)
			# if we want all P3
			# all P3 will be of the form (a, b, c) or (c, b, a). i.e., the edges are (a,b) and (b, c)
			else:
				if [w_ID, u_ID, v_ID] in triples or [v_ID, u_ID, w_ID] in triples:
					#print("This P3 is already included: {}".format((v_ID, u_ID, w_ID)))
					continue
				else:
					#print("Adding the P3: {}".format((w_ID, u_ID, v_ID)))
					triples.append([w_ID, u_ID, v_ID])


		for w_ID in v_adj:
			if w_ID == u_ID:
				#print("{} and {} is the same node".format(w_ID, u_ID))
				continue
			if w_ID in u_adj:
				#print("Triangle detected for the nodes: {}".format((v_ID, u_ID, w_ID)))
				continue
			if get_triples == True:
				if nodes[w_ID]['color'] == nodes[u_ID]['color'] or nodes[u_ID]['color'] == nodes[v_ID]['color'] or nodes[v_ID]['color'] == nodes[w_ID]['color']:
					#print("This P3 is already included: {}".format((u_ID, w_ID, v_ID)))
					continue
				else:
					if (w_ID, u_ID, v_ID) in triples or (u_ID, w_ID, v_ID) in triples:
						#print("This triple is already included: {}".format((u_ID, w_ID, v_ID)))
						continue
					else:
						#print("Adding the triple: {}".format((w_ID, u_ID, v_ID)))
						triples.append((w_ID, u_ID, v_ID))
						leaves.add(u_ID)
						leaves.add(v_ID)
						leaves.add(w_ID)
			else:
				if [w_ID, v_ID, u_ID] in triples or [u_ID, v_ID, w_ID] in triples:
					#print("This P3 is already included: {}".format((u_ID, v_ID, w_ID)))
					continue
				else:
					#print("Adding the P3: {}".format((w_ID, v_ID, u_ID)))
					triples.append([w_ID, v_ID, u_ID])
	return triples, leaves


def P3_regions(l, a = 1):
	'''
		l is a list of P3s (lists)
		a is the minimum number of elements regions need to have in common for them to be
		considered the same.

		returns a list of sets (regions) and a list with amount of P3s per region.
		the idx of the amount list matches the idx of the list of regions.
	'''

	regions = []
	amounts = []
	#curr_key = 0
	curr_value = 1
	while len(l) > 0:
		head, *tail = l
		head = set(head)
		n = -1
		while len(head) > n:
			n = len(head)
			rest = []
			for t in tail:
				t = set(t)
				if len(head.intersection(t)) >= a:
					curr_value += 1
					head |= t
				else:
					rest.append(t)
			tail = rest
		regions.append(head)
		amounts.append(curr_value)
		#curr_key += 1
		curr_value = 1
		l = tail
	return regions, amounts


def P3_distance(lengths, p1, p2):
	'''
		returns the min distance between the path p1 and
		assumes the shortest path between all pairs has been calculated
		no need to run nx.all_pairs_shortest_path_length() here
	'''
	min_dist = float('inf')
	for a in p1:
		d1 = lengths[a][p2[0]]
		d2 = lengths[a][p2[1]]
		d3 = lengths[a][p2[2]]
		min_value = min((d1, d2, d3))
		if min_value < min_dist:
			min_dist = min_value
	return min_dist



def unique_combinations(*l):
	'''
		get all combinations from the lists in l of length n.
	'''
	for c in itertools.combinations(l, 2):
		for pair in itertools.product(*c):
			yield pair

def regions_distance(G, regions, lengths=None):
	'''
		regions is a list of sets containing nodes. each set is its own region.
		G is the graph containing the nodes in the regions.

		return
			dictionary with key mapping to the index of the region in the regions list, the value of which is
			a dictionary with keys mapping to the index of other regions the value of which is the min distance 
			between the regions.
	'''

	if not lengths:
		lengths = dict(nx.all_pairs_shortest_path_length(G))
	region_distances = {}
	k = len(regions)
	for i in range(k):
		region_distances[i] = dict()
		for j in range(i+1, k):
			min_dist = float('inf')

			combinations = unique_combinations(regions[i], regions[j])		
			for c in combinations:
				# check if keys exist
				if c[0] in lengths:
					if c[1] in lengths[c[0]]:
						if lengths[c[0]][c[1]] < min_dist:
							min_dist = lengths[c[0]][c[1]]
			region_distances[i][j] = min_dist
	return region_distances 


def get_P3_data(G, colored_G=None):
	if colored_G:
		P3s, _ = find_all_P3(G, get_triples=True, colored_G=colored_G)
	else:
		P3s, _ = find_all_P3(G, get_triples=True)

	#print("P3s: \n{}".format(P3s))
	#print("length of P3s: {}".format(len(P3s)))
	regions, amounts = P3_regions(P3s)
	regions_distances = regions_distance(G, regions)

	return regions, amounts, regions_distances, len(P3s)


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
	# O(|E||V| + BUILD + |E| + |R||A||B| + |R||A||K| + |A||B|)  -> 	This is triples edit 1
	# O(2|E||V| + |K| + BUILD + |E| + |R||A||B| + |R||A||K| + |A||B|)  -> 	This is triples edit 2
	# TODO: destroy K3s in get_species_triples so that we dont need to rerun it and update the lists
	#		if all K3s are destroyed then we dont need to weight 'bad' and 'unknowns' for delete edges.
	def triples_editing(self, delete_edges=True, add_edges=False, mincut=True, weighted_mincut=True, destroy_K3 = None, G = None):
		"""
		TODO:
			If add_edges is set to true the editing of the original graph will delete aswell as add edges
		"""
		if G:
			copy_G = G.copy()
			color_graph(self._G, copy_G)
		else:
			copy_G = self._G_perturbed.copy()
		#copy_G = copy.deepcopy(self.disturbed_G)
		# copy_G is colored (?) so dont need to pass in a different colored graph here
		triples, species_leaves, triangles, triples_dict = get_species_triples(copy_G, self._G)		# O(|E||V|)

		
		if len(triples) == 0:
			#print("The set of triples is empty!\n")
			return None, None
		# NOTE: If G has less than two nodes an error will occur in the BUILD alg, specifically in stoer_wagner alg.

		if destroy_K3:
			#print("Amount of K3s: {}".format(len(triangles)))
			#print("Destroying K3s...")
			if destroy_K3 == 1:
				destroy_K3s(copy_G, triangles)	# O(|K|)
			elif destroy_K3 == 2:
				# do weighted K3 cuts and make sure to remove exactly 2 (max) edges per K3 (some K3s might share edges)
				pass
			# all K3s destroyed so the list is now empty
			triples, species_leaves, triangles, triples_dict = get_species_triples(copy_G, self._G)		# O(|E||V|)
			if len(triples) == 0:
				return None, None
			#print("Amount of K3s: {}".format(len(triangles)))
			


		B = tools.Build(triples, species_leaves, mincut=mincut, weighted_mincut=weighted_mincut)	
		tree_triples = B.build_tree()															# runtime for weighted mincut BUILD
		
		# NOTE: cut_list is a list of species triples to remove. for example (A, B) could be in cut_list, in which case all triples (A, B, X) for any X need to be removed.
		# weigh the edges
		#print("cut_list: {}".format(B.cut_list))
		if len(B.cut_list) > 0:
			insert_edges = weight_edges(copy_G, triples_dict, triangles, B.cut_list)	#  O(|E| + |R||A||B| + |R||A||K|) 	

			# decide which edge(s) to remove
			for X, Y, Z in triples_dict:							# O(|A|)
				#print("X, Y, Z: {}".format((X, Y, Z)))	
				if (X, Y) in B.cut_list or (Y, X) in B.cut_list:
					#print("{} is in cut_list".format(((X, Y))))
					for a, b, c in triples_dict[(X, Y, Z)]:				# O(|B|)
						

						# check if P3 is still intact.
						if not (copy_G.has_edge(a, c) and copy_G.has_edge(b, c)) or copy_G.has_edge(a, b):
							continue
						# check which edge of (a, c), (b, c) and (a, b) has the highest weight

						#print("{}: {}\t{}: {}\t{}: {}".format((a, c), (copy_G.edges[(a, c)]['bad'], copy_G.edges[(a, c)]['unknowns'], copy_G.edges[(a, c)]['weight']), 
						#											  (b, c), (copy_G.edges[(b, c)]['bad'], copy_G.edges[(b, c)]['unknowns'], copy_G.edges[(b, c)]['weight']), 
						#											  (a, b), (insert_edges[(a, b)]['bad'], insert_edges[(a, b)]['unknowns'], insert_edges[(a, b)]['weight'])))
						best_edge, to_delete, _ = self.edge_to_edit(a, b, c, insert_edges, copy_G)
						
						if to_delete:
							copy_G.remove_edge(*best_edge)
							#print("removed edge {}".format(best_edge))
							broken = True
						else:
							copy_G.add_edge(*best_edge)
							#print("added edge {}".format(best_edge))
							broken = True

		return copy_G, tree_triples

	# O(|V| + |E||V| + |K| + BUILD + |E| + |R||A||B| + |R||A||K| + |A||B||K|^2) (this last term is most likely upper bound since we can remove at most |K| K3s and we remove 1 per call)
	# This isnt fixing the set of triples 100% of the time.
	# It does fix cograph most of the time (because there are very few edges left -> less dense graph)
	def triples_editing2(self, mincut=True, weighted_mincut=True, G = None):
		"""
			1. Get triples and run weighted mincut and get cut_list
			2. The edges of those K3s that overlap any P3s forming the species triples in cut_list, are weighted based on their occurrences.
				Once those K3s are destroyed (by removing 2 of their edges), we can remove any edge from a P3 without risking creating a new P3
			3. Weight the edges ('bad', 'unknowns', 'weight')
			4. if a P3 can't be removed by introducing new bad P3s then we destroy all those K3s that overlap this edge.
		"""
		if G:
			copy_G = G.copy()
			color_graph(self._G, copy_G)		# O(|V|)
		else:
			copy_G = self._G_perturbed.copy()
		triples, species_leaves, K3s, triples_dict = get_species_triples(copy_G, self._G) # O(|E||V|)
		K3_edge_weights = weight_K3s(K3s)												  # O(|K|)
		if len(triples) == 0:
			#print("The set of triples is empty!\n")
			return None, None

		#print("Amount of K3s: {}".format(len(K3s)))
		#print("Destroying K3s...")

		B = tools.Build(triples, species_leaves, mincut=mincut, weighted_mincut=weighted_mincut)	
		tree_triples = B.build_tree()															# weighted mincut BUILD
		#print("weights of K3s: ")
		#print(K3_edge_weights)
		if len(B.cut_list) > 0:
			insert_edges = weight_edges(copy_G, triples_dict, K3s, B.cut_list) 		# O(|E| + |R||A||B| + |R||A||K|)

			for X, Y, Z in triples_dict:											# O(|A|)
				if (X, Y) in B.cut_list or (Y, X) in B.cut_list:
					for a, b, c in triples_dict[(X, Y, Z)]:							# O(|B|)
						# P3 has been broken, skip
						if not (copy_G.has_edge(a, c) and copy_G.has_edge(b, c)) or copy_G.has_edge(a, b):
							continue

						#print("{}: {}\t{}: {}\t{}: {}".format((a, c), (copy_G.edges[(a, c)]['bad'], copy_G.edges[(a, c)]['unknowns'], copy_G.edges[(a, c)]['weight']), 
						#											  (b, c), (copy_G.edges[(b, c)]['bad'], copy_G.edges[(b, c)]['unknowns'], copy_G.edges[(b, c)]['weight']), 
						#											  (a, b), (insert_edges[(a, b)]['bad'], insert_edges[(a, b)]['unknowns'], insert_edges[(a, b)]['weight'])))
						best_edge, to_delete, best_del_edge = self.edge_to_edit(a, b, c, insert_edges, copy_G)
						# only insert if the edge doesn't introduce new triples
						if to_delete:
							copy_G.remove_edge(*best_edge)
							#print("removed edge {}".format(best_edge))
							destroy_K3s_overlapping_edge(copy_G, best_edge, K3s, K3_edge_weights)			# O(|K|^2)
						else:
							# check if 'bad' and 'unknown' is 0, otherwise choose best delete edge
							if insert_edges[best_edge]['bad'] == 0 and insert_edges[best_edge]['unknowns'] == 0:
								copy_G.add_edge(*best_edge)
								#print("added edge {}".format(best_edge))
							else:
								copy_G.remove_edge(*best_del_edge)
								#print("removed edge {}".format(best_edge))
								destroy_K3s_overlapping_edge(copy_G, best_del_edge, K3s, K3_edge_weights)
		triples, species_leaves, K3s, triples_dict = get_species_triples(copy_G)
		#print("Amount of K3s: {}".format(len(K3s)))
		#print("The K3s are: {}".format(K3s))
		return copy_G, tree_triples



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