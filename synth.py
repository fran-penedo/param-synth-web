import cdd
import numpy as np
from numpy import mat
from scipy.integrate import odeint


class CDDMatrix(cdd.Matrix):


    def __init__(self, rows):
        cdd.Matrix.__init__(self, rows)
        self.rep_type = cdd.RepType.INEQUALITY


    def extend(self, rows):
        cdd.Matrix.extend(self, rows)
        return self


    def copy(self):
        return CDDMatrix(self)


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


def bin(s):
    return [s] if s<=1 else bin(s>>1) + [s&1]


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


def partition(constrs):
    part = {}
    for b in constr_it(len(constrs)):
        mat = CDDMatrix([smul(m, row) for row, m in zip(constrs, b)])
        if not pempty(mat):
            part[contobin(b)] = mat

    return part


def refine(part, constrs):
    newpart = partition(constrs)

    prod = ((a.copy().extend(b), a, i1, i2)
            for i1, a in part.items() for i2, b in newpart.items())

    return {i1 if pequal(c, a) else i1 + i2: c
            for c, a, i1, i2 in prod if not pempty(c)}


#######
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def draw(pset):
    return [4, 5]


def amatrix(p, n):
    return mat(list(chunks(p, n))[:-1])


def bmatrix(p, n):
    return mat(p[-n:]).T


def contains(m, p):
    return all(eq[0] + eq[1:].dot(p) >= 0 for eq in m.getA())


class PWASystem(object):


    def __init__(self, sconstrs, psets):
        self.eqs = {index: {'dom': dom, 'pset': pset, 'domnp': mat(dom)}
                    for pset, index, dom
                    in zip(psets, *zip(*partition(sconstrs).items()))}


    def evalf(self, x, t):
        eq = next(e for e in self.eqs.values() if contains(e['domnp'], x))
        p = draw(eq['pset'])
        try:
            n = len(x)
        except TypeError:
            n = 1
        A = amatrix(p, n)
        b = bmatrix(p, n)
        return (A.dot(x) + b).getA()[0]


    def integrate(self, t, x0):
        return odeint(self.evalf, x0, t)


