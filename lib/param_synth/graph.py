from pygraphviz import AGraph

def build_graph(ts):
    g = AGraph(directed=True)
    g.add_edges_from(ts.ts.todok().keys())
    g.graph_attr['overlap'] = 'scalexy'
    for n, s in zip(g.nodes(), ts._pwa.states):
        n.attr['label'] = s
    return g

def path_graphs(leaf):
    print leaf.node.ts.toNUSMV()
    gs = [build_graph(leaf.node.ts)]
    t = leaf

    while t.parent is not None:
        g = build_graph(t.parent.node.ts)
        for a, b in t.node.path[len(t.parent.node.path):]:
            g.get_edge(t.node.ts.states.index(a),
                       t.node.ts.states.index(b)).attr['color'] = 'red'
        gs.append(g)
        t = t.parent


    return gs

def draw_path_graphs(leaf, prefix=""):
    gs = path_graphs(leaf)

    for i, g in enumerate(gs):
        g.draw(prefix + "{}.png".format(i), prog='neato')
