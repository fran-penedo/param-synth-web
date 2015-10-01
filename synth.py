import cdd
import numpy as np
import scipy
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO
from scipy.integrate import odeint
from scipy.sparse import lil_matrix
from scipy.linalg import det
#from scipy.spatial import Delaunay
import operator
from math import factorial
from subprocess import Popen, PIPE
import copy
import re
import tempfile
import os


class _CDDMatrix(object):

    def __init__(self, m):
        self._m = m

    def __len__(self):
        return len(self._m)

    def __getitem__(self, key):
        return self._m.__getitem__(key)

    def __setitem__(self, key, value):
        return self._m.__setitem__(key, value)

    def __delitem__(self, key):
        return self._m.__delitem__(key)

    def canonicalize(self):
        return self._m.canonicalize()

    @property
    def col_size(self):
        return self._m.col_size

    @property
    def row_size(self):
        return self._m.row_size

    @property
    def obj_type(self):
        return self._m.obj_type

    @obj_type.setter
    def obj_type(self, value):
        self._m.obj_type = value

    @property
    def rep_type(self):
        return self._m.rep_type

    @property
    def lin_set(self):
        return self._m.lin_set

    @rep_type.setter
    def rep_type(self, value):
        self._m.rep_type = value

    def extend(self, b):
        lin = []
        notlin = []
        for i in range(b.row_size):
            if i in b.lin_set:
                lin.append(b[i])
            else:
                notlin.append(b[i])

        self._extend(lin, True)
        self._extend(notlin, False)

    def _extend(self, rows, linear):
        if len(rows) > 0:
            self._m.extend(rows, linear)

    def copy(self):
        return _CDDMatrix(self._m.copy())

    def __deepcopy__(self, memo):
        return self.copy()

    def __str__(self):
        return self._m.__str__()


def CDDMatrix(rows, ineq=True):
    m = cdd.Matrix(rows)
    if not ineq:
        m.rep_type = cdd.RepType.GENERATOR
        m = cdd.Polyhedron(m).get_inequalities()
    m.rep_type = cdd.RepType.INEQUALITY
    return _CDDMatrix(m)


def vrep(m):
    return _CDDMatrix(cdd.Polyhedron(m._m).get_generators())


def vrep_pts(m):
    return [v[1:] for v in vrep(m)]


def extend(a, b):
    x = a.copy()
    x.extend(b)
    return x


class CDDMatrixUnion(object):

    def __init__(self, ms):
        if type(ms) is not list:
            ms = [ms]
        self._components = ms

    def components(self):
        return self._components

    def copy(self):
        return CDDMatrixUnion([x.copy() for x in self.components()])

    def extend(self, x):
        for y in self.components():
            y.extend(x)

        self.prune()

    def prune(self):
        self._components = [m for m in self.components() if not pempty(m)]


def smul(s, v):
    return [s * i for i in v]


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)))


def pequal(m1, m2):
    return set(vrep(m1)) == set(vrep(m2))


def pempty(m):
    if len(m.copy().canonicalize()[0]) > 0:
        return True

    m.obj_type = cdd.LPObjType.MAX
    m.obj_func = [0 for i in range(m.col_size)]
    lp = cdd.LinProg(m._m)
    lp.solve()
    return lp.status == cdd.LPStatusType.INCONSISTENT


def pdiff(a, b):
    if isinstance(b, CDDMatrixUnion):
        raise NotImplementedError()
    return pinters(a, CDDMatrixUnion([CDDMatrix(smul(-1, v)) for v in b]))

def pinters(a, b):
    if isinstance(b, CDDMatrixUnion):
        if isinstance(a, CDDMatrixUnion):
            return CDDMatrixUnion([comp for m in b.components()
                                   for comp in extend(a, m).components()])
        else:
            return extend(b, a)
    else:
        return extend(a, b)


def bin(s):
    return [s] if s <= 1 else bin(s >> 1) + [s & 1]


def binp(s, p):
    b = bin(s)
    return [0 for i in range(p - len(b))] + b


def contobin(c):
    return ''.join(['1' if i == 1 else '0' for i in c])


def bintocon(b):
    return [1 if c == '1' else -1 for c in b]


def constr_it(n):
    i = 0
    m = 2**n
    while i < m:
        b = binp(i, m.bit_length() - 1)
        yield map(lambda x: x if x == 1 else -1, b)
        i += 1


def partition(poly, constrs):
    part = {}
    for b in constr_it(len(constrs)):
        mat = CDDMatrix([smul(m, row) for row, m in zip(constrs, b)])
        p = pinters(poly, mat)
        if not pempty(p):
            part[contobin(b)] = p

    return part


def refine(poly, part, constrs):
    newpart = partition(poly, constrs)

    prod = ((extend(a, b), a, i1, i2)
            for i1, a in part.items() for i2, b in newpart.items())

    return {i1 if pequal(c, a) else i1 + i2: c
            for c, a, i1, i2 in prod if not pempty(c)}


def simp_vol(s):
    return abs(det(np.array([v - s[0] for v in s[1:]])) / factorial(len(s) - 1))


def volume(m):
    if isinstance(m, CDDMatrixUnion):
        # FIXME THIS IS WRONG!!
        return sum(volume(x) for x in m.components())
    else:
        #pts = vrep_pts(m)
        #dt = Delaunay(pts)
        # return sum(simp_vol(s) for s in dt.points[dt.simplices])
        pts = vrep_pts(m)
        return reduce(operator.mul, [x for x in
                                     np.amax(pts, 0) - np.amin(pts, 0)
                                     if x > 0.01])


#######
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def draw(pset):
    # TODO
    return [4, 5]


def amatrix(p, n):
    return np.array(list(chunks(p, n))[:-1])


def bmatrix(p, n):
    return np.array(p[-n:]).T


def affine_eval(p, x, n):
    A = amatrix(p, n)
    b = bmatrix(p, n)
    return A.dot(x) + b


def contains(m, p):
    if m is not None:
        return all(eq[0] + eq[1:].dot(p) >= 0 for eq in m)
    else:
        return False


class Equation(object):

    def __init__(self, dom, pset, domnp):
        self._dom = dom
        self._pset = pset
        self._domnp = domnp

    @property
    def dom(self):
        return self._dom

    @dom.setter
    def dom(self, value):
        self._dom = value

    @property
    def pset(self):
        return self._pset

    @pset.setter
    def pset(self, value):
        self._pset = value

    @property
    def domnp(self):
        return self._domnp

    @domnp.setter
    def domnp(self, value):
        self._domnp = value


class PWASystem(object):
    OUT = 'out'

    def __init__(self, poly, sconstrs, psets):
        self.eqs = {index: Equation(dom, pset, np.array(dom))
                    for pset, index, dom
                    in zip(psets,
                           *zip(*sorted(partition(poly, sconstrs).items())))}
        self.eqs[PWASystem.OUT] = Equation(
            CDDMatrix([smul(-1, h) for h in poly]), None, None)
        self.n = len(sconstrs[0]) - 1

    @property
    def states(self):
        return sorted(self.eqs.keys())

    def evalf(self, x, t=0):
        eq = next(e for e in self.eqs.values() if contains(e.domnp, x))
        p = draw(eq.pset)
        return affine_eval(p, x, self.n)

    # There's actually no need for this
    def integrate(self, t, x0):
        return odeint(self.evalf, x0, t)

    def connected_poly(self, l1, l2):
        if l1 == PWASystem.OUT:
            return False

        xl2 = self.eqs[l2].dom
        eq1 = self.eqs[l1]
        for pset in eq1.pset.components():
            hull = CDDMatrix([[1] + list(affine_eval(p[1:], v[1:], self.n))
                              for p in vrep(pset)
                              for v in vrep(eq1.dom)], False)
            if l2 == PWASystem.OUT:
                if any(not pempty(pinters(hull, CDDMatrix([hs])))
                       for hs in xl2):
                    return True
            else:
                if not pempty(pinters(hull, xl2)):
                    return True

        return False

    def disconnect(self, l1, l2):
        Xl1 = self.eqs[l1].dom
        Xl2 = self.eqs[l2].dom
        if l1 == l2:
            trans = [CDDMatrix([[np.array(h).dot(v)] +
                                list(- np.array(h).dot(reshape_x(v)))
                                for v in vrep_pts(Xl1)])
                     for h in sample_h(Xl1.col_size - 1)]
            ps = [pinters(self.eqs[l1].pset, t) for t in trans]
            self.eqs[l1].pset = max([(p, volume(p)) for p in ps],
                                    key=lambda pair: pair[1])[0]
        else:
            self.eqs[l1].pset = pinters(self.eqs[l1].pset,
                                        punderset(Xl1, Xl2))

    def connected_dreal(self, l1, l2):
        if l1 == PWASystem.OUT:
            return False

        Xl1 = self.eqs[l1].dom
        Xl2 = self.eqs[l2].dom
        for Pl1 in self.eqs[l1].pset.components():
            if l2 == PWASystem.OUT:
                if any(dreal_check_sat(dreal_connect_smt(Xl1, Pl1, CDDMatrix([hs]))) for hs in Xl2):
                    return True
            else:
                smt = dreal_connect_smt(Xl1, Pl1, Xl2)
                #if l1 == "10" and l2 == "10":
                #    print smt
                if dreal_check_sat(smt):
                    return True

        return False

    connected = connected_dreal

FP_REGEXP = "[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

def dreal_find_p(smt):
    check, out = _dreal_check_sat(smt, verbose=True)
    if check:
        r = re.compile("p([0-9]+) : \[ ENTIRE \] = \[(%s), (%s)\]" % (FP_REGEXP, FP_REGEXP))
        p_tuples = sorted([(int(i), (float(a) + float(b)) / 2)
                           for i, a, b in r.findall(out)])
        return zip(*p_tuples)[1]
    else:
        return None


def dreal_check_sat(smt):
    return _dreal_check_sat(smt)[0]


def _dreal_check_sat(smt, verbose=False):
    t = tempfile.NamedTemporaryFile(suffix=".smt2", delete=False)
    t.write(smt)
    t.close()
    process = ["lib/dReal/bin/dReal"]
    if verbose:
        process.append("--model")
    process.append(t.name)
    ps = Popen(process, stdout=PIPE, stderr=PIPE,
               env=dict(os.environ, LD_LIBRARY_PATH="lib/dReal/lib"))
    out, err = ps.communicate()
    outlines = out.splitlines()
    if outlines[0].startswith("delta-sat") or outlines[-1].startswith("delta-sat"):
        return True, out
    elif outlines[0].startswith("unsat") or outlines[-1].startswith("unsat"):
        return False, out
    else:
        print smt
        print out
        print err
        raise Exception()


def dreal_linear(eq, prefix):
    return " ".join(["(* %s%d %f)" %
                     (prefix, i - 1, eq[i]) for i in range(1, len(eq))])


def dreal_poly(poly, prefix, slack=0):
    return "(and %s)" % " ".join(["(%s %f (+ %f %s))" %
                                  ("=" if i in poly.lin_set else "<=",
                                   slack, eq[0], dreal_linear(eq, prefix))
                                  for i, eq in enumerate(poly)])


def dreal_connect_smt(Xl1, Pl1, Xl2, PExcl=None, excl_slack=0):
    n = len(Xl1[0]) - 1
    if PExcl is None:
        PExcl = []
    out = StringIO()
    print >>out, "(set-logic QF_NRA)"
    for i in range(n):
        print >>out, "(declare-fun x%d () Real)" % i
        print >>out, "(declare-fun xn%d () Real)" % i

    for i in range(n * n + n):
        print >>out, "(declare-fun p%d () Real)" % i

    print >>out, "(assert %s)" % dreal_poly(Xl1, "x")
    print >>out, "(assert %s)" % dreal_poly(Xl2, "xn", 0.01)

    if len(Pl1) > 0:
        print >>out, "(assert %s)" % dreal_poly(Pl1, "p")

    if len(PExcl) > 0:
        print >>out, "(assert (not %s))" % dreal_poly(PExcl, "p", -excl_slack)

    for i in range(n):
        print >>out, "(assert (= xn%d (+ %s)))" % \
            (i,
             " ".join(["(* p%d x%d)" % (n * i + j, j) for j in range(n)] +
                      ["p%d" % (n * n + i)]))

    print >>out, "(check-sat)"
    print >>out, "(exit)"

    s = out.getvalue()
    out.close()
    return s


def punderset(Xl1, Xl2):
    return CDDMatrixUnion(
        [CDDMatrix([[-h[0]] + list(- np.array(h[1:]).dot(reshape_x(x[1:])))
                    for x in vrep(Xl1)])
         for h in Xl2])


def reshape_x(x):
    return np.hstack((scipy.linalg.block_diag(*[x for i in range(len(x))]),
                      np.identity(len(x))))


def sample_h(n):
    return np.vstack([np.identity(n), -np.identity(n)])


class PWATS(object):

    def __init__(self, pwa, init=[], ltl=None):
        self._pwa = pwa
        self._ltl = ltl
        self.ts = self.build_ts(pwa)
        self._init = [self.states.index(i) for i in init]

    # TODO take a look at multiprocessing. Should be easy to parallelize
    def build_ts(self, pwa):
        m = np.array([[1 if pwa.connected(l1, l2) else 0
                       for l2 in pwa.states] for l1 in pwa.states])
        return lil_matrix(m)

    @property
    def states(self):
        return self._pwa.states

    @property
    def init(self):
        return self._init

    def remove_link(self, i, j):
        if self.ts[i, j] == 0:
            return []

        self._pwa.disconnect(self.states[i], self.states[j])
        rem = self.update_connected(i)
        self.ts[i, j] = 0
        if (i, j) not in rem:
            rem.append((i, j))
        remall = rem + self.remove_blocking()
        return remall

    def update_connected(self, i):
        remove = [j for j in self.ts[i, :].nonzero()[1]
                  if not self._pwa.connected(self.states[i], self.states[j])]
        self.ts[i, remove] = 0
        return zip([i for x in remove], remove)

    def isblocking(self, i):
        return self.issink(i) and (self.ts[:, i].getnnz() > 0 or
                                   i in self.init)

    def issink(self, i):
        return self.ts[i, ].getnnz() == 0

    def remove_blocking(self):
        try:
            r = next(i for i in range(len(self.states))
                     if self.isblocking(i) and i not in self.init)
            removed = sum([self.remove_link(i, r)
                           for i in self.ts[:, r].nonzero()[0]], [])
            return removed + self.remove_blocking()
        except StopIteration:
            return []

    def modelcheck(self):
        ps = Popen('lib/nusmv/NuSMV', stdin=PIPE, stdout=PIPE)
        out = ps.communicate(self.toNUSMV())[0]
        try:
            return parse_nusmv(out)
        except:
            print self.toNUSMV()
            raise Exception()

    def isfeasible(self):
        ltl = self._ltl
        self._ltl = '! (' + ltl + ')'
        check, trace = self.modelcheck()
        self._ltl = ltl
        return not check

    def __eq__(self, other):
        # TODO
        # Should I compare the pwa? The paper implies it's not necessary, but
        # I'm not sure
        return (self.ts - other.ts).nnz == 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(str(self.ts))

    def copy(self):
        return copy.deepcopy(self)

    def toNUSMV(self):
        out = StringIO()
        print >>out, "MODULE main"
        print >>out, "VAR"
        print >>out, 'state : %s;' % nusmv_statelist(self.states)
        print >>out, 'ASSIGN'
        print >>out, 'init(state) := %s;' % \
            nusmv_statelist([self.states[i] for i in self._init])
        print >>out, 'next(state) := '
        print >>out, 'case'

        for i in range(len(self.states)):
            if not self.issink(i):
                print >>out, 'state = %s : %s;' % \
                    ("s" + self.states[i],
                     nusmv_statelist([self.states[j]
                                      for j in self.ts[i, ].nonzero()[1]]))

        print >>out, 'TRUE : state;'
        print >>out, 'esac;'
        if self._ltl is not None:
            print >>out, 'LTLSPEC %s' % self._ltl

        s = out.getvalue()
        out.close()
        return s


def nusmv_statelist(l):
    return '{%s}' % ', '.join(map(lambda x: "s" + x, l))


def parse_nusmv(out):
    if out.find('true') != -1:
        return True, []
    elif out.find('Parser error') != -1:
        print out
        raise Exception()
    else:
        lines = out.splitlines()
        start = next(i for i in range(len(lines))
                     if lines[i].startswith('Trace Type: Counterexample'))
        loop = next(i for i in range(len(lines))
                    if lines[i].startswith('-- Loop starts here'))

        p = re.compile('state = s([0,1]+)')
        matches = (p.search(line) for line in lines[start:])
        chain = [m.group(1) for m in matches if m is not None]
        if loop == len(lines) - 4:
            chain.append(chain[-1])

        trace = [(x, y) for x, y in zip(chain, chain[1:])]

        return False, trace


class Tree(object):

    def __init__(self, node, children=None):
        self._node = node
        if children is None:
            self._children = []
        else:
            self._children = children

    def __getitem__(self, i):
        return self._children[i]

    def isleave(self):
        return len(self._children) == 0

    @property
    def node(self):
        return self._node

    def add_child(self, tree):
        self._children.append(tree)

    def add_children(self, children):
        self._children.extend(children)

    def contains(self, item, f):
        if f(self._node, item):
            return True
        else:
            return any((child.contains(item, f) for child in self._children))

    def __str__(self):
        return self.pprint(0)

    def pprint(self, indent):
        return ''.join([(' ' * indent + '{}').format(l) for l in
                        ['{%s,\n' % self.node.path.__str__(),
                         '%s,\n' % self.node.path.__str__(),
                         '[\n%s\n' % ',\n'.join([x.pprint(indent + 2)
                                                 for x in self._children]),
                         ']}']
                        ])


class Node(object):

    def __init__(self, path, ts, feas, trace):
        self._path = path
        self._ts = ts
        self._feas = feas
        self._trace = trace

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        self._path = value

    @property
    def ts(self):
        return self._ts

    @ts.setter
    def ts(self, value):
        self._ts = value

    @property
    def feas(self):
        return self._feas

    @feas.setter
    def feas(self, value):
        self._feas = value

    @property
    def trace(self):
        return self._trace

    @trace.setter
    def trace(self, value):
        self._trace = value


def compare_nodes(a, b):
    return a.ts == b


def leaves(tree):
    if tree.isleave():
        yield tree
    else:
        for c in tree._children:
            for l in leaves(c):
                yield l


PATH = '\033[91m'
ENDC = '\033[0m'


def synthesize(ts, depth=-1):
    root = Tree(Node([], ts, None, None))
    children, feas = _synthesize(ts, [], depth, set([ts]))
    root.node.feas = feas
    root.add_children(children)
    return root


def _synthesize(t, path, depth, memo):
    #print len(memo)
    if any((t.isblocking(q) for q in t.init)) or not t.isfeasible():
        return [], False
    elif depth != 0:
        check, trace = t.modelcheck()
        itrace = [(t.states.index(i), t.states.index(j)) for i, j in trace]
        if check:
            print PATH + "FOOOOOOOOOO" + ENDC
            return [], True
        else:
            tnexts = [t.copy() for l in trace]
            removed = [[(tn.states[a], tn.states[b]) for a, b in tn.remove_link(i, j)]
                       for tn, (i, j) in zip(tnexts, itrace)]
            children = [Tree(Node(path + rem, tn, None, l))
                        for tn, rem, l in zip(tnexts, removed, itrace)
                        if tn not in memo]

            memo.update([c.node.ts for c in children])

            for child in children:
                cs, feas = _synthesize(child.node.ts, child.node.path,
                                       depth - 1, memo)
                child.node.feas = feas
                child.add_children(cs)

            return children, True
    else:
        return [], True

