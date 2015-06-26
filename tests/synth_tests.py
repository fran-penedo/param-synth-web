from nose.tools import *
import synth as s


def bin_test():
    assert_equal(s.binp(5, 8), [0, 0, 0, 0, 0, 1, 0, 1])


def partition_test():
    constr = [[1, -1, 0], [2, -1, 0],
              [1, 0, -1], [2, 0, -1]]
    p = zip(*s.partition(constr))[0]
    x = set([''.join(['1' if i == 1 else '0' for i in b]) for b in p])
    y = set(['0000', '0001', '0011', '0100', '0101', '0111',
             '1100', '1101', '1111'])
    assert_set_equal(x, y)

