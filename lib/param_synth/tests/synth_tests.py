from nose.tools import *
from synth import *
import numpy as np
from numpy.testing import *


def bin_test():
    assert_equal(binp(5, 8), [0, 0, 0, 0, 0, 1, 0, 1])


def contains_test():
    m = np.array([[-1, 1, 0], [2, -1, 0],
              [-1, 0, 1], [2, 0, -1]])
    assert_true(contains(m, [1.5, 1.5]))
    assert_false(contains(m, [0, 0]))


def partition_test():
    constr = [[1, -1, 0], [2, -1, 0],
              [1, 0, -1], [2, 0, -1]]
    poly = CDDMatrix([[0, 1, 0], [0, 0, 1], [3, -1, 0], [3, 0, -1]])
    x = set(partition(poly, constr).keys())
    y = set(['0000', '0001', '0011', '0100', '0101', '0111',
             '1100', '1101', '1111'])
    assert_set_equal(x, y)


def pequal_test():
    a = CDDMatrix([[1, -1], [2, -2], [2, -1], [4, -2]])
    b = CDDMatrix([[2, -1]])

    assert_false(pequal(a, b))
    assert_false(pequal(b, a))

    c = CDDMatrix([[3, -3]])
    assert_true(pequal(a, c))
    assert_true(pequal(c, a))


def refine_test():
    constr = [[1, -1, 0], [2, -1, 0],
              [1, 0, -1], [2, 0, -1]]
    add = [[0, -1, 1]]
    poly = CDDMatrix([[0, 1, 0], [0, 0, 1], [3, -1, 0], [3, 0, -1]])
    x = set(refine(poly, partition(poly, constr), add).keys())
    y = set(['00000', '00001', '0001', '0011', '0100',
             '01010', '01011', '0111',
             '1100', '1101', '11111', '11110'])
    assert_set_equal(x, y)


def volume_test():
    pts = [[1,a,b,c,d]
           for a in [0,1] for b in [0,1] for c in [0,1] for d in [0,1]]
    assert_almost_equal(volume(CDDMatrix(pts, False)), 1)

def pwa_test():
    constr = [[1, -1, 0], [2, -1, 0],
              [1, 0, -1], [2, 0, -1]]
    poly = CDDMatrix([[0, 1, 0], [0, 0, 1], [3, -1, 0], [3, 0, -1]])
    psets = [None for i in range(9)]
    psets[4] = 'P'
    pwa = PWASystem(poly, constr, psets)
    eq = pwa.eqs['0101']
    assert_equal(eq.pset, 'P')


def integrate_test():
    constr = [[5, 1]]
    psets = [None, None]
    poly = CDDMatrix([[1, -7], [1, 200000]], False)
    pwa = PWASystem(poly, constr, psets)

    assert_equal(pwa.evalf(3, 0), 17)

    assert_allclose(list(pwa.integrate(np.array([0, 2]), np.array([3]))[:,0]),
                      [3, (17 * np.exp(8) - 5)/4], rtol=1e-2)


def pwats_test():
    constr = [[0, 1], [-1, 1]]
    psets = [CDDMatrix([[1, 1, 0.5]], False),
             CDDMatrix([[1, 1, 1]], False),
             CDDMatrix([[1, -2, 2.5]], False)]
    poly = CDDMatrix([[1, -1], [1, 2]], False)
    psets = map(CDDMatrixUnion, psets)
    pwa = PWASystem(poly, constr, psets)
    ts = PWATS(pwa)

    #print pwa.states
    desired = np.array([[1, 1, 0, 0], [0, 0, 1, 0], [1, 1, 0, 1], [0, 0, 0, 0]])
    assert_array_equal(ts.ts.toarray(), desired)


def pwats_toNUSMV_test():
    constr = [[0, 1], [-1, 1]]
    psets = [CDDMatrix([[1, 1, 0.5]], False),
             CDDMatrix([[1, 1, 1]], False),
             CDDMatrix([[1, -1, 1.5]], False)]
    poly = CDDMatrix([[1, -1], [1, 2]], False)
    psets = map(CDDMatrixUnion, psets)
    pwa = PWASystem(poly, constr, psets)
    ts = PWATS(pwa, init=["00"], ltl='F state = s10')

    #print ts.toNUSMV()


def nusmv_statelist_test():
    assert_equal(nusmv_statelist(["01", "10"]), '{s01, s10}')
    assert_equal(nusmv_statelist(["01"]), '{s01}')

def parse_nusmv_test():
    with open('tests/parse_nusmv_test.txt', 'r') as f:
        check, trace = parse_nusmv(f.read())

    assert_false(check)
    assert_list_equal(trace,
                      [('00', '10'), ('10', '11'), ('11', '00')])

def pwa_disconnect_test():
    constr = [[0, 1], [-1, 1]]
    psets = [CDDMatrix([[1, 1, 0.5]], False),
             CDDMatrix([[1, 1, 1]], False),
             CDDMatrix([[1, -1, 1.5]], False)]
    poly = CDDMatrix([[1, -1], [1, 2]], False)
    psets = map(CDDMatrixUnion, psets)
    pwa = PWASystem(poly, constr, psets)

    assert_true(pwa.connected('00', '10'))
    pwa.disconnect('00', '10')
    assert_false(pwa.connected('00', '10'))


def synthesize_test():
    constr = [[0, 1], [-1, 1]]
    psets = [CDDMatrix([[1, 1, 0.5]], False),
             CDDMatrix([[1, 1, 1]], False),
             CDDMatrix([[1, -1, 1.5]], False)]
    poly = CDDMatrix([[1, -1], [1, 2]], False)
    psets = map(CDDMatrixUnion, psets)
    pwa = PWASystem(poly, constr, psets)
    ts = PWATS(pwa, init=["00"], ltl='G state = s00')

    tree = synthesize(ts)
    #print tree


def dreal_connect_smt_test():
    Xl1 = CDDMatrix([[1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]], False)
    Pl1 = CDDMatrix([[1, 1, 0, 0, 1, 0, 0], [1, 1, 0, 0, 1, 1, 0],
                    [1, 1, 0, 0, 1, 0, 1], [1, 1, 0, 0, 1, 1, 1]], False)
    Xl2 = CDDMatrix([[1, 1, 0], [1, 1, 1], [1, 2, 0], [1, 2, 1]], False)

    print dreal_connect_smt(Xl1, Pl1, Xl2)


def dreal_find_p_test():
    Xl1 = CDDMatrix([[1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]], False)
    Pl1 = CDDMatrix([[1, 1, 0, 0, 1, 0, 0], [1, 1, 0, 0, 1, 1, 0],
                    [1, 1, 0, 0, 1, 0, 1], [1, 1, 0, 0, 1, 1, 1]], False)
    Xl2 = CDDMatrix([[1, 1, 0], [1, 1, 1], [1, 2, 0], [1, 2, 1]], False)

    smt = dreal_connect_smt(Xl1, Pl1, Xl2)
    p = dreal_find_p(smt)
    print p

def pwa_connected_dreal_test():
    constr = [[0, 1], [-1, 1]]
    psets = [CDDMatrix([[1, 1, 0.5]], False),
             CDDMatrix([[1, 1, 1]], False),
             CDDMatrix([[1, -1, 1.5]], False)]
    poly = CDDMatrix([[1, -1], [1, 2]], False)
    psets = map(CDDMatrixUnion, psets)
    pwa = PWASystem(poly, constr, psets)

    assert_true(pwa.connected_dreal('00', '10'))
    pwa.disconnect('00', '10')
    assert_false(pwa.connected_dreal('00', '10'))

