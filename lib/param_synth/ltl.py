
EXPR = 0
NOT = 1
AND = 2
OR = 3
NEXT = 4
ALWAYS = 5
EVENTUALLY = 6


class Signal(object):

    def __init__(self, labels, f, bounds=None):
        self._labels = labels
        self._f = f
        self.bounds = bounds

    def signal(self, model, t):
        vs = map(lambda l: model.getVarByName(l(t)), self._labels)
        if any(var is None for var in vs):
            return None
        else:
            return self._f(vs)


class Formula(object):

    """Docstring for Formula. """

    def __init__(self, operator, args, bounds=[0, 0]):
        self.op = operator
        self.args = args
        self.bnd = bounds

    def _hexpr(self):
        return 0

    def _hnot(self):
        return 0

    def _hand(self):
        return max(map(lambda f: f.horizon(), self.args))

    def _hor(self):
        return self._hand()

    def _halways(self):
        return self.bnd[1] + self.args[0].horizon()

    def _hnext(self):
        return 1 + self.args[0].horizon()

    def _heventually(self):
        return self._halways()

    def horizon(self):
        return {
            EXPR: self._hexpr,
            NOT: self._hnot,
            AND: self._hand,
            OR: self._hor,
            NEXT: self._hnext,
            ALWAYS: self._halways,
            EVENTUALLY: self._heventually
        }[self.op]()


def robustness(formula, model, t=0):
    return {
        EXPR: lambda: formula.args[0].signal(model, t),
        NOT: lambda: -robustness(formula.args[0], model, t),
        AND: lambda: min(map(lambda f: robustness(f, model, t), formula.args)),
        OR: lambda: max(map(lambda f: robustness(f, model, t), formula.args)),
        NEXT: lambda: robustness(formula.args[0], model, t + 1),
        ALWAYS: lambda: min(map(
            lambda j: robustness(formula.args[0], model, t + j),
            range(formula.bnd[0], formula.bnd[1] + 1))),
        EVENTUALLY: lambda: max(map(
            lambda j: robustness(formula.args[0], model, t + j),
            range(formula.bnd[0], formula.bnd[1] + 1)))
    }[formula.op]()
