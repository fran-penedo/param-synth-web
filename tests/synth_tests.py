from nose.tools import *
from synth import *
import numpy as np
from numpy.testing import *


def bin_test():
    assert_equal(binp(5, 8), [0, 0, 0, 0, 0, 1, 0, 1])


def contains_test():
    m = np.mat([[-1, 1, 0], [2, -1, 0],
              [-1, 0, 1], [2, 0, -1]])
    assert_true(contains(m, [1.5, 1.5]))
    assert_false(contains(m, [0, 0]))


def partition_test():
    constr = [[1, -1, 0], [2, -1, 0],
              [1, 0, -1], [2, 0, -1]]
    x = set(partition(constr).keys())
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
    x = set(refine(partition(constr), add).keys())
    y = set(['00000', '00001', '0001', '0011', '0100',
             '01010', '01011', '0111',
             '1100', '1101', '11111', '11110'])
    assert_set_equal(x, y)


def pwa_test():
    constr = [[1, -1, 0], [2, -1, 0],
              [1, 0, -1], [2, 0, -1]]
    psets = [None for i in range(9)]
    psets[4] = 'P'
    pwa = PWASystem(constr, psets)
    eq = pwa.eqs['0101']
    assert_equal(eq['pset'], 'P')


def integrate_test():
    constr = [[5, 1]]
    psets = [None]
    pwa = PWASystem(constr, psets)

    assert_equal(pwa.evalf(3, 0), 17)

    print pwa.integrate(np.array([0, 2]), np.array([3]))[:,0]
    print [3, (17 * np.exp(8) - 5)/4]

    assert_allclose(list(pwa.integrate(np.array([0, 2]), np.array([3]))[:,0]),
                      [3, (17 * np.exp(8) - 5)/4], rtol=1e-2)
