import cdd

def foo():
    mat = cdd.Matrix([[1.0, 0, -1], [0, 0, 1], [2, -1, 0], [0, 1, 0]])
    mat.rep_type = cdd.RepType.INEQUALITY
    poly = cdd.Polyhedron(mat)
    print poly

    ext = poly.get_generators()
    print ext

    mat = cdd.Matrix([[1, -1], [-2, 1]])
    mat.obj_type = cdd.LPObjType.MAX
    mat.obj_func = [0, 0]
    lp = cdd.LinProg(mat)
    lp.solve()
    print lp.status == cdd.LPStatusType.INCONSISTENT
    print lp.primal_solution


    mat = cdd.Matrix([[1, -1], [2, -1]])
    mat.obj_type = cdd.LPObjType.MAX
    mat.obj_func = [0, 0]
    lp = cdd.LinProg(mat)
    lp.solve()
    print lp.status == cdd.LPStatusType.INCONSISTENT
    print lp.primal_solution


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


def partition(constr):
    part = []
    for b in constr_it(len(constr)):
        mat = cdd.Matrix([[a * m for a in row] for row, m in zip(constr, b)])
        mat.obj_type = cdd.LPObjType.MAX
        mat.obj_func = [0 for i in constr[0]]
        lp = cdd.LinProg(mat)
        lp.solve()
        if lp.status != cdd.LPStatusType.INCONSISTENT:
            part.append((b, mat))

    return part

