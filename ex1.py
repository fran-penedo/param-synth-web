
from synth import *

def build_pwa():
    A = [1, 0.9, -0.05, 0.05, 0.9]
    bvs = [[-1.5, -1.5], [-1.5, 1.5], [1.5, -1.5], [1.5, 1.5]]
    P = CDDMatrixUnion(CDDMatrix([A + v for v in bvs], False))
    poly = CDDMatrix([[0, 1, 0], [0, 0, 1], [3, -1, 0], [3, 0, -1]])
    constr = [[1, -1, 0], [2, -1, 0],
              [1, 0, -1], [2, 0, -1]]
    pwa = PWASystem(poly, constr, [P for i in range(9)])

    return pwa

def test():
    pwa = build_pwa()
    ts = PWATS(pwa, init=["1111"],
               ltl='G ((F state = s1111) & (F state = s0001) & \
                    (F state = s1100) & (F state = s0100) & \
                    (F state = s0000) & (F state = s1101) & \
                    (F state = s0011) & (F state = s0111) & \
                    ! (state = s0101))')
    tree = synthesize(ts)
    return tree


if __name__ == '__main__':
    print test()



