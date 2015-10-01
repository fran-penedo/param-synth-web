import networkx as nx
import matplotlib.pyplot as plt

def draw_graph(sparse, f):
    g = nx.DiGraph(sparse)
    nx.draw(g)

def draw_tree_path(leaf):
    t = leaf
    gs = []

    while True:
        g = nx.DiGraph(t.node.ts.m)
        for a, b in t.node.path[len(t.parent.node.path):]:
            g[a][b]['color'] = 'red'
        gs.add(g)

        if t.parent is None:
            break
