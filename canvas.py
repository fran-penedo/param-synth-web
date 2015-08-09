from flask import Flask, render_template
import lib.param_synth.synth as synth
from synth import CDDMatrix

app = Flask(__name__)
app.debug = True

@app.route("/")
def index():
    return render_template('index.html')

def partition(poly_vs, constrs):
    poly = CDDMatrix([[1] + v for v in poly_vs], False)

    constrs_vect = [[b[1] - a[1], a[0] - b[0]] for a, b in constrs]
    constrs_vect = [[- a[0][0] * v[0] - a[0][1] * v[1]] + v
                    for cons, v in constrs, constrs_vect]

    p = synth.partition(poly, constrs_vect)
    p_simple = {}
    for k, m in p:
        p_simple[k] = [c for c in m]

    return p_simple


def center(m):
    vs = vrep_pts(m)



if __name__ == '__main__':
    app.run()
