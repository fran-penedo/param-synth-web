from flask import Flask, render_template, request, jsonify
import lib.param_synth.synth as synth
from lib.param_synth.synth import CDDMatrix, vrep_pts
import numpy as np

app = Flask(__name__)
app.debug = True

@app.route("/")
def index():
    return render_template('index.html')

@app.route("/partition", methods=['POST'])
def partition():
    data = request.get_json(True)
    poly = data["poly"]
    constrs = data["constrs"]
    print constrs

    return jsonify(partition=build_partition(poly, constrs))


def build_partition(poly_vs, constrs):
    poly = CDDMatrix([[1] + v for v in poly_vs], False)

    constrs_vect = [[b[1] - a[1], a[0] - b[0]] for a, b in constrs]
    constrs_vect = [[- cons[0][0] * v[0] - cons[0][1] * v[1]] + v
                    for cons, v in zip(constrs, constrs_vect)]

    p = synth.partition(poly, constrs_vect)
    print p

    p_simple = [{"name": k,
                 "constrs": [c for c in m],
                 "centroid": centroid(m)}
                for k, m in p.items()]

    return p_simple


def centroid(m):
    return np.average(vrep_pts(m), axis=0).tolist()


if __name__ == '__main__':
    app.run()
