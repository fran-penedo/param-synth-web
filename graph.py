from pygraphviz import AGraph

def draw_graph(sparse, f):
    g = AGraph(directed=True, strict=True, concentrate=True)
    g.add_edges_from(sparse.todok().keys())
    g.draw(f, prog='neato')
