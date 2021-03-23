import unittest
import networkx as nx
import sys, os
sys.path.insert(0, os.path.abspath('..'))

from src.tools.GraphTools import *

G1 = nx.Graph()
G2 = nx.Graph()
G3 = nx.Graph()

G1.add_edges_from([(0, 1), (0, 2), (1, 4), (2, 3), (2, 4)])
G1_P3 = [[2, 0, 1], [2, 4, 1], [3, 2, 4], [3, 2, 0], [0, 2, 4], [0, 1, 4]]
P3, _ = find_all_P3(G1)

def assertListofList(l1, l2):
	'''
		Checks if the lists in the given lists l1, l2
		share the same elements (order not important)
	'''
	a = {(frozenset(item)) for item in l1}
	b = {(frozenset(item)) for item in l2}
	return True if a == b else False

class TestPathTools(unittest.TestCase):


	def test_find_all_P3(self):
		P3, _ = find_all_P3(G1)
		self.assertTrue(assertListofList(P3, G1_P3))
	

	def test_merge_regions(self):
		pass

	def test_merge_regions_lite(self):
		pass
		

	def test_overlapping_P3_amount(self):
		pass


	def test_P3_distance(self):
		pass

	def test_P3_distance_matrix(self):
		pass

#print("Current working dir: {}".format(os.getcwd()))
if __name__ == "__main__":
	unittest.main()
	