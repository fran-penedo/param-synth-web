import networkx as nx
import matplotlib.pyplot as plt

def draw_graph(sparse, f):
    g = nx.DiGraph(sparse)
    nx.draw(g)

def path_graphs(leaf):
    gs = [nx.DiGraph(leaf.node.ts.ts)]
    t = leaf

    while t.parent is not None:
        g = nx.DiGraph(t.parent.node.ts.ts)
        #for a, b in t.node.path[len(t.parent.node.path):]:
        #    g[t.node.ts.states.index(a)][t.node.ts.states.index(b)]['color'] = 'red'
        gs.append(g)
        t = t.parent

    return gs

def draw_path_graphs(leaf, prefix=""):
    gs = path_graphs(leaf)

    for i, g in enumerate(gs):
        plt.figure()
        nx.draw(g)
        plt.savefig(prefix + "{}.png".format(i), bbox_inches='tight')
