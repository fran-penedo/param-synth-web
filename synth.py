import cdd
import numpy as np
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO
from numpy import mat
from scipy.integrate import odeint
from scipy.sparse import lil_matrix
import re


def CDDMatrix(rows, ineq=True):
    m = cdd.Matrix(rows)
    if not ineq:
        m.rep_type = cdd.RepType.GENERATOR
        m = cdd.Polyhedron(m).get_inequalities()
    m.rep_type = cdd.RepType.INEQUALITY
    return m


def vrep(m):
    return cdd.Polyhedron(m).get_generators()


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


def smul(s, v):
    return [s * i for i in v]


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)))


def pequal(m1, m2):
    p1 = cdd.Polyhedron(m1)
    p2 = cdd.Polyhedron(m2)

    return set(p1.get_generators()) == set(p2.get_generators())


def pempty(m):
    if len(m.copy().canonicalize()[0]) > 0:
        return True

    m.obj_type = cdd.LPObjType.MAX
    m.obj_func = [0 for i in range(m.col_size)]
    lp = cdd.LinProg(m)
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


#######
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def draw(pset):
    return [4, 5]


def amatrix(p, n):
    return mat(list(chunks(p, n))[:-1])


def bmatrix(p, n):
    return mat(p[-n:]).T


def affine_eval(p, x, n):
    A = amatrix(p, n)
    b = bmatrix(p, n)
    return (A.dot(x) + b).getA()[0]


def contains(m, p):
    return all(eq[0] + eq[1:].dot(p) >= 0 for eq in m.getA())


class PWASystem(object):

    def __init__(self, poly, sconstrs, psets):
        self.eqs = {index: {'dom': dom,
                            'pset': pset,
                            'domnp': mat(dom)}
                    for pset, index, dom
                    in zip(psets,
                           *zip(*sorted(partition(poly, sconstrs).items())))}
        self.n = len(sconstrs[0]) - 1

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
        pass


class PWATS(object):

    def __init__(self, pwa, init=None, ltl=None):
        self._pwa = pwa
        self._init = init
        self._ltl = ltl
        self.ts = self.build_ts(pwa)

    # TODO take a look at multiprocessing. Should be easy to parallelize
    def build_ts(self, pwa):
        m = np.array([[1 if pwa.connected(l1, l2) else 0
                       for l2 in pwa.states()] for l1 in pwa.states()])
        return lil_matrix(m)

    def states(self):
        return self._pwa.states()

    def remove_link(self, i, j):
        if self.ts[i, j] == 0:
            return

        self._pwa.disconnect(i, j)
        self.ts[i, j] = 0

    def isblocking(self, i):
        return self.ts[i, ].getnnz() == 0

    # TODO change parameter sets
    def remove_blocking(self):
        try:
            r = next(i for i in range(len(self.states()))
                     if self.isblocking(i))
            self.ts[:, r] = 0
            self.remove_blocking()
        except StopIteration:
            return

    def modelcheck(self):
        ps = Popen('lib/nusmv/NuSMV', stdin=PIPE, stdout=PIPE)
        out = ps.comunicate(self.toNUSMV())
        return parse_nusmv(out)

    def toNUSMV(self):
        out = StringIO()
        print >>out, "MODULE main"
        print >>out, "VAR"
        print >>out, 'state : %s;' % nusmv_statelist(self.states())
        print >>out, 'ASSIGN'
        print >>out, 'init(state) := %s;' % nusmv_statelist(self._init)
        print >>out, 'next(state) := '
        print >>out, 'case'

        for i in range(len(self.states())):
            if not self.isblocking(i):
                print >>out, 'state = %s : %s;' % \
                    ("s" + self.states()[i],
                     nusmv_statelist([self.states()[j]
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

