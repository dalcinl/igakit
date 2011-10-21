import numpy as np
from igakit.transform import transform
from igakit import igalib

_api = {
    0 : igalib.bsp,
    1 : igalib.crv,
    2 : igalib.srf,
    3 : igalib.vol,
    }

__all__ = ['NURBS', 'transform']

class NURBS(object):

    def __init__(self, control, knots, weights=None):
        #
        if (isinstance(control, (list, tuple))
            and len(control) == 2):
            assert weights is None
            control, weights = control
        #
        control = np.asarray(control, dtype='d')
        assert control.ndim >= 2
        nsd = control.shape[-1]
        assert nsd >= 2 and nsd <= 4
        if nsd < 4:
            shape = control.shape[:-1]
            Cw = np.zeros(shape+(4,), dtype='d')
            Cw[...,0:nsd] = control
            w = np.ones(shape, dtype='d')
            Cw[...,-1] = w
            if weights is not None:
                w = np.asarray(weights, dtype='d')
                assert w.shape == shape
                Cw *= w[...,np.newaxis]
            control = Cw
        else:
            assert weights is None
            control = np.ascontiguousarray(control)
        #
        knots = [np.asarray(k, dtype='d') for k in knots]
        assert len(knots) >= 1
        assert len(knots) <= 3
        for k in knots:
            assert k.ndim == 1
            assert k.size >= 4
        knots = tuple(knots)
        #
        assert control.ndim - 1 == len(knots)
        N = [n-1 for n in control.shape[:-1]]
        M = [len(k)-1 for k in knots]
        for n, m in zip(N, M):
            degree = m - n - 1
            assert degree >= 1
        #
        self._cntrl = control
        self._knots = knots

    @property
    def control(self):
        return self._cntrl

    @property
    def knots(self):
        return self._knots

    @property
    def dim(self):
        return self.control.ndim-1

    @property
    def shape(self):
        return self.control.shape[:-1]

    @property
    def degree(self):
        N = [n-1 for n in self.shape]
        M = [len(k)-1 for k in self.knots]
        return tuple([m-n-1 for (m, n) in zip(M, N)])

    @property
    def points(self):
        Cw = self.control
        w = self.weights
        return Cw[...,:-1] / w[...,np.newaxis]

    @property
    def weights(self):
        return self.control[...,-1]

    #

    def copy(self):
        nrb = NURBS.__new__(type(self))
        nrb._cntrl = nrb.control.copy()
        nrb._knots = tuple([k.copy() for k in self.knots])

    def clone(self):
        nrb = NURBS.__new__(type(self))
        nrb._cntrl = self.control
        nrb._knots = self.knots
        return nrb

    def transform(self, trans):
        if not isinstance(trans, transform):
            trans = transform(trans)
        Cw = trans(self.control)
        self._cntrl = np.ascontiguousarray(Cw)
        return self

    def transpose(self, *axes):
        if not axes:
            axes = range(self.dim)[::-1]
        else:
            assert len(axes) == self.dim
        caxes = list(axes)+[self.dim]
        control = self.control.transpose(caxes).copy()
        knots = [self.knots[i] for i in axes]
        #
        self._cntrl = control
        self._knots = tuple(knots)
        return self

    def reverse(self, *axes):
        def CntRev(C, axis):
            dim = C.ndim - 1
            index = list(np.index_exp[:,:,:][:dim])
            index[axis] = slice(None, None, -1)
            return C[index]
        def KntRev(p, U):
            U0, U1 = U[p], U[-p-1]
            return (U1+U0)-U
        allaxes = range(self.dim)
        if not axes: axes = allaxes
        control = self.control
        knots = list(self.knots)
        degree = self.degree
        for axis in axes:
            axis = allaxes[axis]
            control = CntRev(control, axis)
            knots[axis] = KntRev(degree[axis], knots[axis])
        control = control.copy()
        #
        self._cntrl = control
        self._knots = tuple(knots)
        return self

    def refine(self, *uvw):
        assert len(uvw) == self.dim
        def arg(U):
            if U is None: U = []
            return np.asarray(U, dtype='d').ravel()
        uvw = [arg(U) for U in uvw]
        #
        arglist = []
        for p, U in zip(self.degree, self.knots):
            arglist.extend([p,U])
        arglist.append(self.control)
        arglist.extend(uvw)
        #
        RefineKnotVector = _api[self.dim].RefineKnotVector
        result = RefineKnotVector(*arglist)
        control = result[-1]
        knots = result[:-1]
        #
        self._cntrl = control
        self._knots = knots
        return self

    def elevate(self, *rst):
        assert len(rst) == self.dim
        def arg(t):
            if t is None: t = 0
            return int(t)
        rst = [arg(t) for t in rst]
        for t in rst: assert t >= 0
        #
        arglist = []
        for p, U in zip(self.degree, self.knots):
            arglist.extend([p, U])
        arglist.append(self.control)
        arglist.extend(rst)
        #
        DegreeElevate = _api[self.dim].DegreeElevate
        result = DegreeElevate(*arglist)
        control = result[-1]
        knots = result[:-1]
        #
        self._cntrl = control
        self._knots = knots
        return self

    def extract(self, axis, value):
        assert self.dim > 1
        bsp = _api[0]
        axis = range(self.dim)[axis]

        p = self.degree[axis]
        U = self.knots[axis]
        u = float(value)
        assert u >= U[p]
        assert u <= U[-1-p]

        span = bsp.FindKnotSpan(p, U, u)
        mult = bsp.Multiplicity(p, U, u, span)
        times = p - mult
        if mult < p:
            times = p - mult
            uvw = [None] * self.dim
            uvw[axis] = [value] * times
            nrb = self.clone().refine(*uvw)
            control = nrb.control
            knots = list(nrb.knots)
            U = nrb.knots[axis]
            span = bsp.FindKnotSpan(p, U, u)
        else:
            control = self.control
            knots = list(self.knots)

        offset = span - p
        if u < U[span]: offset =  0
        if u > U[span]: offset = -1
        index = list(np.index_exp[:,:,:][:self.dim])
        index[axis] = offset
        control = control[index].copy()
        del knots[axis]

        nrb = NURBS.__new__(type(self))
        nrb._cntrl = control
        nrb._knots = tuple(knots)
        return nrb

    def boundary(self, axis, side):
        assert self.dim > 1
        assert side in (0, 1)
        axis = range(self.dim)[axis]
        p = self.degree[axis]
        U = self.knots[axis]
        if side == 0:
            value = U[p]
        else:
            value = U[-1-p]
        return self.extract(axis, value)

    def evaluate(self, *uvw):
        assert len(uvw) == self.dim
        uvw = [np.asarray(U, dtype='d') for U in uvw]
        #
        arglist = []
        for p, U in zip(self.degree, self.knots):
            arglist.extend([p, U])
        arglist.append(self.control)
        arglist.extend(uvw)
        #
        Evaluate = _api[self.dim].Evaluate
        Cw = Evaluate(*arglist)
        C = Cw[...,:-1] / Cw[...,-1,np.newaxis]
        return C

    #

    def plot(self):
        pass
