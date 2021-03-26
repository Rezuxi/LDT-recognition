import networkx as nx


G1 = nx.Graph()
G2 = nx.Graph()
G3 = nx.Graph()

G1.add_edges_from([(0, 1), (0, 2), (1, 4), (2, 3), (2, 4)])
G1_P3 = [[2, 0, 1], [2, 4, 1], [3, 2, 4], [3, 2, 0], [0, 2, 4], [0, 1, 4]]


# colors: 0 is blue, 1 is green, 2 is red
G2_nodes = [(0, {"color": 0}), (1, {"color": 1}), (2, {"color": 2}), (3, {"color": 1}),
			(4, {"color": 0}), (5, {"color": 2}), (6, {"color": 2}), (7, {"color": 1}),
			(8, {"color": 0}), (9, {"color": 2}), (10, {"color": 2}), (11, {"color": 1}),
			(12, {"color": 2}), (13, {"color": 0}), (14, {"color": 1}), (15, {"color": 0}),
			(16, {"color": 2})
			]
G2_edges = [(0, 1), (0, 2), (1, 4), (2, 3), (2, 4), (4, 5), (4, 6), (4, 7), (5, 7), (6, 7),
			(7, 8), (7, 9), (3, 10), (10, 11), (11, 12), (11, 13), (13, 14), (14, 15), (15, 16)
			]

# Expected results for G2
G2_colored_P3 = [(1, 0, 2), (1, 4, 2), (0, 2, 3), (4, 2, 3), (1, 4, 6), (1, 4, 5), (2, 4, 7),
				 (4, 7, 9), (5, 7, 8), (6, 7, 8), (8, 7, 9),
				 (10, 11, 13), (13, 11, 12),
				 (14, 15, 16)
				]

G2_P3_regions = [{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {10, 11, 12, 13}, {14, 15, 16}]
G2_P3_region_amounts = [11, 2, 1]
G2_region_distances = {0: {1: 1, 2: 4}, 1 : {2 : 1}}

G2.add_nodes_from(G2_nodes)
G2.add_edges_from(G2_edges)
