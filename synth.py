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
from scipy.spatial import Delaunay
from math import factorial
from subprocess import Popen, PIPE
import copy
import re


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
    def obj_type(self):
        return self._m.obj_type

    @obj_type.setter
    def obj_type(self, value):
        self._m.obj_type = value

    @property
    def rep_type(self):
        return self._m.rep_type

    @rep_type.setter
    def rep_type(self, value):
        self._m.rep_type = value

    def extend(self, rows, linear=False):
        self._m.extend(rows, linear)

    def copy(self):
        return _CDDMatrix(self._m.copy())

    def __deepcopy__(self, memo):
        return self.copy()

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


def extend(a, b, linear=False):
    x = a.copy()
    x.extend(b, linear)
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

    def extend(self, x, linear=False):
        if isinstance(x, CDDMatrixUnion):
            for y in x.components():
                self._extend_single(y, linear)

        else:
            self._extend_single(x, linear)

        self.prune()

    def prune(self):
        self._components = [m for m in self.components() if not pempty(m)]

    def _extend_single(self, x, linear):
        for y in self.components():
            y.extend(x, linear)


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


def pinters(a, b):
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
    pts = vrep_pts(m)
    dt = Delaunay(pts)
    return sum(simp_vol(s) for s in dt.points[dt.simplices])


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
    return all(eq[0] + eq[1:].dot(p) >= 0 for eq in m)


class PWASystem(object):

    def __init__(self, poly, sconstrs, psets):
        self.eqs = {index: {'dom': dom,
                            'pset': pset,
                            'domnp': np.array(dom)}
                    for pset, index, dom
                    in zip(psets,
                           *zip(*sorted(partition(poly, sconstrs).items())))}
        self.n = len(sconstrs[0]) - 1

    @property
    def states(self):
        return sorted(self.eqs.keys())

    def evalf(self, x, t=0):
        eq = next(e for e in self.eqs.values() if contains(e['domnp'], x))
        p = draw(eq['pset'])
        return affine_eval(p, x, self.n)

    # There's actually no need for this
    def integrate(self, t, x0):
        return odeint(self.evalf, x0, t)

    def connected(self, l1, l2):
        xl2 = self.eqs[l2]['dom']
        eq1 = self.eqs[l1]
        for pset in eq1['pset'].components():
            hull = CDDMatrix([[1] + list(affine_eval(p[1:], v[1:], self.n))
                              for p in vrep(pset)
                              for v in vrep(eq1['dom'])], False)
            if not pempty(pinters(hull, xl2)):
                return True

        return False

    def disconnect(self, l1, l2):
        Xl1 = self.eqs[l1]['dom']
        Xl2 = self.eqs[l2]['dom']
        if l1 == l2:
            trans = [CDDMatrix([[np.array(h).dot(v)] +
                                (- np.array(h).dot(reshape_x(v)))
                                for v in vrep_pts(Xl1)])
                     for h in sample_h()]
            ps = [pinters(self.eqs[l1]['pset'], t) for t in trans]
            self.eqs[l1]['pset'] = max([(p, volume(p)) for p in ps],
                                       lambda pair: pair[1])[0]
        else:
            punderset = CDDMatrixUnion(
                [CDDMatrix([[-h[0]] + list(- np.array(h[1:]).dot(reshape_x(x[1:])))
                            for x in vrep(Xl1)])
                for h in Xl2])
            self.eqs[l1]['pset'] = pinters(self.eqs[l1]['pset'], punderset)


def reshape_x(x):
    return np.hstack((scipy.linalg.block_diag(*[x for i in range(len(x))]),
                      np.identity(len(x))))

def sample_h(n):
    return np.identity(n)



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
        return rem + self.remove_blocking()

    def update_connected(self, i):
        remove = [j for j in self.ts[i,:].nonzero()[1]
                  if not self._pwa.connected(self.states[i], self.states[j])]
        self.ts[i, remove] = 0
        return zip([i for x in remove], remove)

    def isblocking(self, i):
        return self.ts[i, ].getnnz() == 0 and (self.ts[:,i].getnnz() > 0 or
                                               i in self.init)

    def remove_blocking(self):
        try:
            r = next(i for i in range(len(self.states))
                     if self.isblocking(i) and i not in self.init)
            removed = sum([self.remove_link(i, r)
                           for i in self.ts[:,r].nonzero()[0]], [])
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

    def __eq__(self, other):
        # TODO
        # Should I compare the pwa? The paper implies it's not necessary, but
        # I'm not sure
        return (self.ts - other.ts).nnz == 0

    def __ne__(self, other):
        return not self.__eq__(other)

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
            if not self.isblocking(i):
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
        loopstate = p.search(lines[loop + 2]).group(1)

        trace = [(x, y) for x, y in zip(chain, chain[1:])] + \
            [(chain[-1], loopstate)]

        return False, trace


class Tree(object):

    def __init__(self, node, children=None):
        self._node = node
        if children is None:
            self._children = []
        else:
            self._children = children

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
        return '{%s, [%s]}' % (self._node.__str__(),
                               ', '.join([child.__str__()
                                          for child in self._children]))


def compare_nodes(a, b):
    return a[1] == b


def synthesize(ts):
    root = Tree(([], ts))
    stack = [root]

    while len(stack) > 0:
        cur = stack.pop()
        t = cur._node[1]
        if not any((t.isblocking(q) for q in t.init)):
            check, trace = t.modelcheck()
            itrace = [(t.states.index(i), t.states.index(j)) for i, j in trace]
            if not check:
                tnexts = [t.copy() for l in trace]
                removed = [tn.remove_link(i, j)
                           for tn, (i, j) in zip(tnexts, itrace)]
                children = [Tree((rem, tn)) for tn, rem in zip(tnexts, removed)
                            if not root.contains(tn, compare_nodes)]
                cur.add_children(children)
                stack.extend(children)

    return root
