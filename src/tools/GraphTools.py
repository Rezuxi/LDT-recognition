import networkx as nx
import asymmetree.tools as tools
import asymmetree.cograph as cg
from asymmetree.tools.GraphTools import disturb_graph
import random
import copy


#####################################################################################
#																					#
#						Induced paths on 3 vertices (P3)							#
#																					#
#####################################################################################


def find_all_P3(G, get_triples=False, colored_G=None):
	'''
		Finds all connected triples that don't form a triangle
		if triples = True (the graph needs to be colored), we also check that the colors are distinct for each node in a triple (a, b, c)

		returns the P3s as a list of lists
		returns the triples as a list of tuples
	'''
	leaves = set()
	triples = []
	#nodes = None

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
				if nodes[w_ID]['color'] == nodes[v_ID]['color']:
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
				if nodes[w_ID]['color'] == nodes[u_ID]['color']:
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

		returns a list of sets (regions) and a dictionary where the key is the idx of the region and
		value is the number of P3s in that region.
	'''
	regions = []
	amounts = {}
	curr_key = 0
	curr_value = 0
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
		regions.append(first)
		amounts[curr_key] = curr_value
		curr_key += 1
		curr_value = 0
		l = tail
	return regions, amounts

'''
def P3_regions_lite(l):
	"""
		convert to graph and count components(regions)
		this method is faster but only applies to the case where 
		sharing at least 1 node merges the regions.

		returns the "regions" of the P3s in l.
	"""
	G = nx.Graph()
	def to_edges(ls):
		"""
			for the P3 (a, b, c) it adds the edges (a, b) and (b, c)
		"""
		it = iter(ls)
		last = next(it)
		for curr in it:
			yield last, curr
			last = curr
	for n in l:
		G.add_edges_from(to_edges(n))

	return nx.connected_components(G)

'''


def overlapping_P3_amount(G):
	'''
		returns the number of overlapping P3s in G
	'''
	P3, _ = find_all_P3(G)
	k = len(P3)
	if k == 0:
		return 0

	count = 0

	for i in range(k):
		for j in range(i+1, k):
			if any(item in P3[i] for item in P3[j]):
				count += 1
	return count


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

			

def regions_distance(G, regions):
	pass



def P3_distance_matrix(G):
	'''
		each row is a P3 at index i
		each entry is the distance from P3 at index i to the P3 at index j (the column) 
	'''
	# get all P3
	P3, _ = find_all_P3(G)
	k = len(P3)

	lengths = nx.all_pairs_shortest_path_length(G)
	# create k*k matrix all entries set to 0
	m = [[0 for x in range(k)] for y in range(k)]
	for i in range(k):
		for j in range(i+1, k):
			m[i][j] = P3_distance(lengths, P3[i], P3[j])
	return m


#####################################################################################
#																					#
#						   		Graph Editing tools									#
#																					#
#####################################################################################


def is_compatible(G, colored_G = None):
	if colored_G:
		triples, leaves = find_all_P3(G, get_triples=True, colored_G=colored_G)
	else:
		triples, leaves = find_all_P3(G, get_triples=True)

	B = tools.Build(triples, leaves)
	tree_triples = B.build_tree()

	if tree_triples:
		return True
	else:
		return False

def is_cograph(G):
	cotree = cg.Cotree.cotree(G)
	if cotree:
		return True
	return False





class InvestigateGraph:

	def __init__(self, G):
		'''
			G is an LDT graph
		'''
		self._G = G
		self._G_perturbed = None
		self._is_cograph = True
		self._is_compatible = True

		self._cographEdit_to_LDT = 0
		self._triplesEdit_to_LDT = 0
		
		self._triplesEdit_broke_cograph = 0
		self._cographEdit_broke_compatibility = 0
		self._triplesEdit_fixed_cograph = 0
		self._cographEdit_fixed_triples = 0

	def perturb_graph(self, i_rate = None, d_rate = None):
		# by default have random values for deletion/insertion rate
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
		# make sure the disturbed graph is not an LDT
		if self._is_cograph and self._is_compatible:
			print("adding noise again!")
			self.perturb_graph()


	# only deletes edges with highest weight for every triple to be removed
	def triples_editing(self):
		copy_G = self._G_perturbed.copy()
		#copy_G = copy.deepcopy(self.disturbed_G)
		triples, leaves = find_all_P3(copy_G, get_triples=True, colored=self._G)
		
		if len(triples) == 0:
			print("The set of triples is empty!\n")
			return None, None
		# NOTE: If G has less than two nodes an error will occur in the BUILD alg, specifically in stoer_wagner alg.

		B = tools.Build(triples, leaves, mincut=True)
		tree_triples = B.build_tree()
		#print("Cut list: {}".format(B.cut_list))
		# set weights to 0
		for a, b, c in triples:
			copy_G[a][c]['weight'] = 0
			copy_G[b][c]['weight'] = 0

		# set weights to the edges based on how often they appear in the set of triples
		if len(B.cut_list) > 0:
			for a, b, c in triples:
				copy_G[a][c]['weight'] += 1
				copy_G[b][c]['weight'] += 1

			# decide which edge(s) to remove
			for a, b, c in triples:
				if (a, b) in B.cut_list or (b, a) in B.cut_list:
					'''
						check if both edges exist. If not then we dont need to do anything since
						1 of the edges has already been cut.
					'''
					if not (copy_G.has_edge(a, c) and copy_G.has_edge(b, c)):
						continue
					# check which edge of (a, c) and (b, c) has the highest weight
					if copy_G[a][c]['weight'] >= copy_G[b][c]['weight']:
						# remove edge (a, c)
						copy_G.remove_edge(a, c)
						#print("Removing edge ({}, {})".format(a, c))
					else:
						copy_G.remove_edge(b, c)
						#print("Removing edge ({}, {})".format(b, c))
		return copy_G, tree_triples


	def cograph_editing(self):
		edited_G = cg.edit_to_cograph(self._G_perturbed)
		return edited_G



if __name__ == "__main__":
	pass