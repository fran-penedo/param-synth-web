from flask import Flask, render_template, request, jsonify, send_file
import lib.param_synth.synth as synth
from lib.param_synth.synth import CDDMatrix, CDDMatrixUnion, vrep_pts
import numpy as np
import lib.param_synth.graph as graph
import math

app = Flask(__name__)
app.debug = True

IMG_DIR = "temp/"

@app.route("/")
def index():
    return render_template('index.html')


@app.route("/partition", methods=['POST'])
def partition():
    data = request.get_json(True)
    poly = data["poly"]
    constrs = data["constrs"]

    return jsonify(partition=build_partition(poly, constrs))


@app.route("/synthesize", methods=['POST'])
def synthesize():
    data = request.get_json(True)
    constrs = data["constrs"]
    poly = data["poly"]
    psets = data["pars"]
    init = data["init"]
    spec = data["spec"]

    pwa = synth.PWASystem(CDDMatrix([[1] + v for v in poly], False),
                    to_hrep(constrs),
                    [CDDMatrixUnion(
                        CDDMatrix([[1] +
                                   [float(i) for i in pset["A"].split(" ")] +
                                   v for v in pset["b_space"]], False))
                    for pset in psets])
    ts = synth.PWATS(pwa, init=[init], ltl=spec)

    tree = synth.synthesize(ts)
    leaf = next(l for l in synth.leaves(tree) if l.node.feas)

    prefix = "foo"
    imgs = graph.draw_path_graphs(leaf, prefix, directory=IMG_DIR)

    cur = leaf
    psets = [psets_vertices(cur.node.ts._pwa)]
    while cur.parent is not None:
        psets.append(psets_vertices(cur.parent.node.ts._pwa))
        cur = cur.parent

    result = {"psets":list(reversed(psets)), "imgs": imgs}
    return jsonify(result=result)


def psets_vertices(pwa):
    vs = []
    for l in pwa.states:
        if l == synth.PWASystem.OUT: continue
        pset = pwa.eqs[l].pset
        cents = [centroid(p)[4:] for p in pset.components()]
        verts = [vrep_pts(p)[:,4:].tolist() for p in pset.components()]
        verts = [
            sorted(v,
                   key=lambda p: math.atan2(p[1]-cent[1],p[0]-cent[0]))
             for v, cent in zip(verts, cents)]
        vs.append(verts)

    return vs


@app.route("/synthesize2", methods=['POST'])
def synthesize2():
    vs1 = [[2,0], [0, 2], [-2, 0]]
    vs2 = [[3,0], [0, 4], [-2, 0]]

    vs = [[vs1, vs2] for i in range(4)]
    return jsonify(psets=vs)


@app.route("/get_image", methods=['GET'])
def get_image():
    f = request.args['file']
    return send_file(IMG_DIR + f)


def to_hrep(constrs):
    constrs_vect = [[b[1] - a[1], a[0] - b[0]] for a, b in constrs]
    constrs_vect = [[- cons[0][0] * v[0] - cons[0][1] * v[1]] + v
                    for cons, v in zip(constrs, constrs_vect)]
    return constrs_vect


def build_partition(poly_vs, constrs):
    poly = CDDMatrix([[1] + v for v in poly_vs], False)

    constrs_vect = to_hrep(constrs)

    p = sorted(synth.partition(poly, constrs_vect).items())

    p_simple = [{"name": k,
                 "constrs": [c for c in m],
                 "centroid": centroid(m)}
                for k, m in p]

    return p_simple


def centroid(m):
    return np.average(vrep_pts(m), axis=0).tolist()


if __name__ == '__main__':
    app.run()
