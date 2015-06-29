from nose.tools import *
import synth as s
from synth import Matrix
import cdd


def bin_test():
    assert_equal(s.binp(5, 8), [0, 0, 0, 0, 0, 1, 0, 1])


def partition_test():
    constr = [[1, -1, 0], [2, -1, 0],
              [1, 0, -1], [2, 0, -1]]
    x = set(s.partition(constr).keys())
    y = set(['0000', '0001', '0011', '0100', '0101', '0111',
             '1100', '1101', '1111'])
    assert_set_equal(x, y)


def pequal_test():
    a = Matrix([[1, -1], [2, -2], [2, -1], [4, -2]])
    b = Matrix([[2, -1]])

    assert_false(s.pequal(a, b))
    assert_false(s.pequal(b, a))

    c = Matrix([[3, -3]])
    assert_true(s.pequal(a, c))
    assert_true(s.pequal(c, a))


def refine_test():
    constr = [[1, -1, 0], [2, -1, 0],
              [1, 0, -1], [2, 0, -1]]
    add = [[0, -1, 1]]
    x = set(s.refine(s.partition(constr), add).keys())
    y = set(['00000', '00001', '0001', '0011', '0100',
             '01010', '01011', '0111',
             '1100', '1101', '11111', '11110'])
    assert_set_equal(x, y)
