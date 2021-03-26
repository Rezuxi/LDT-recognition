import unittest
import networkx as nx
import sys, os
sys.path.insert(0, os.path.abspath('..'))

import tests.expectedResults as er
from src.tools.GraphTools import *





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
		#P3, _ = find_all_P3(er.G1)
		#self.assertTrue(assertListofList(P3, er.G1_P3))
		pass

	def test_P3_regions(self):
		P3s, _ = find_all_P3(er.G2, get_triples=True)
		regions, amounts = P3_regions(P3s)

		self.assertCountEqual(regions, er.G2_P3_regions, msg="The regions are incorrect!")
		self.assertCountEqual(amounts, er.G2_P3_region_amounts, msg="The count per regions is incorrect!")

	def test_P3_distance(self):
		pass
		

	def test_regions_distance(self):
		P3s, _ = find_all_P3(er.G2, get_triples=True)
		regions, amounts = P3_regions(P3s)

		distances = regions_distance(er.G2, regions)
		self.assertCountEqual(distances, er.G2_region_distances, msg="The distance between each region is incorrect!")


	def test_P3_distance_matrix(self):
		pass

#print("Current working dir: {}".format(os.getcwd()))
if __name__ == "__main__":
	unittest.main()
	