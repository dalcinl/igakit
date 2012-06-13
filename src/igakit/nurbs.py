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
    """
    An abstract NURBS object class.

    This NURBS class allows for the definition of B-spline or NURBS
    curves/surfaces/volumes by specifying a control point array, a
    sequence of knot vectors and optional rational weights.

    Parameters
    ----------
    control : array_like or 2-tuple of array_like
        Control points and optional rational weights.
    knots : sequence of array_like
        Knot vectors. The number of knot vectors will define what kind
        of NURBS object is created (1=curve, 2=surface, 3=volume).
    weights : array_like, optional
       Rational weights. If weights are omitted, the object will be
       non-rational (B-spline).

    Attributes
    ----------
    dim : int
        Parametric dimension of the NURBS object {1,2,3}
    shape : tuple of ints
        Number of control points in each parametric dimension.
    degree : tuple of ints
        Polynomial degrees in each parametric dimension.
    knots : tuple of numpy.ndarray
        Knot vectors in each parametric dimension.
    control : numpy.ndarray
        Control points in homogeneous 4D space (includes rational weights).
    points : numpy.ndarray
        Control points projected into 3D space.
    weigths : numpy.ndarray
        Rational weigths.

    Examples
    --------

    Create a quarter circle NURBS curve with 2D control points and
    rational weigths and check error:

    >>> C = [[0, 1], [1, 1], [1, 0]] # 3x2 grid of 2D control points
    >>> w = [1, np.sqrt(2)/2, 1]     # rational weigths
    >>> U = [0,0,0, 1,1,1]           # knot vector
    >>> crv = NURBS(C, [U], w)
    >>> u = np.linspace(0,1,1000)
    >>> xyz = crv.evaluate(u)
    >>> x, y, z = xyz.T
    >>> r = np.sqrt(x**2+y**2)
    >>> np.allclose(r, 1, rtol=0, atol=1e-15)
    True
    >>> np.allclose(z, 0, rtol=0, atol=1e-15)
    True

    Create a quarter circle NURBS curve with homogeneous 4D control
    points and check error:

    >>> wgt = np.sqrt(2)/2
    >>> Cw = np.zeros((3,4))
    >>> Cw[0,:] = [0.0, 1.0, 0.0, 1.0]
    >>> Cw[1,:] = [wgt, wgt, 0.0, wgt]
    >>> Cw[2,:] = [1.0, 0.0, 0.0, 1.0]
    >>> crv = NURBS(Cw, [U])
    >>> u = np.linspace(0,1,1000)
    >>> xyz = crv.evaluate(u)
    >>> x, y, z = xyz.T
    >>> r = np.sqrt(x**2+y**2)
    >>> np.allclose(r, 1, rtol=0, atol=1e-15)
    True
    >>> np.allclose(z, 0, rtol=0, atol=1e-15)
    True

    Create a random B-spline curve:

    >>> C = np.random.rand(3,3) # 3D control points
    >>> U = [0,0,0, 1,1,1]      # knot vector
    >>> crv = NURBS(C, [U])
    >>> crv.dim
    1
    >>> crv.shape
    (3,)
    >>> crv.degree
    (2,)
    >>> np.allclose(crv.knots[0], U, rtol=0, atol=1e-15)
    True
    >>> np.allclose(crv.points,   C, rtol=0, atol=1e-15)
    True
    >>> np.allclose(crv.weights,  1, rtol=0, atol=1e-15)
    True

    Create a random B-spline surface:

    >>> C = np.random.rand(3,2,3) # 3x2 grid of 3D control points
    >>> U = [0,0,0, 1,1,1]        # knot vector
    >>> V = [0,0, 1,1]            # knot vector
    >>> srf = NURBS(C, [U,V])
    >>> srf.dim
    2
    >>> srf.shape
    (3, 2)
    >>> srf.degree
    (2, 1)
    >>> np.allclose(srf.knots[0], U, rtol=0, atol=1e-15)
    True
    >>> np.allclose(srf.knots[1], V, rtol=0, atol=1e-15)
    True
    >>> np.allclose(srf.points,   C, rtol=0, atol=1e-15)
    True
    >>> np.allclose(srf.weights,  1, rtol=0, atol=1e-15)
    True

    Create a random B-spline volume:

    >>> C = np.random.rand(3,2,7,3) # 3x2x7 grid of 3D control points
    >>> U = [0,0,0, 1,1,1]          # knot vector
    >>> V = [0,0, 1,1]              # knot vector
    >>> W = [0]*4+[0.25, 0.5, 0.5]+[1]*4
    >>> vol = NURBS(C, [U,V,W])
    >>> vol.dim
    3
    >>> vol.shape
    (3, 2, 7)
    >>> vol.degree
    (2, 1, 3)
    >>> np.allclose(vol.knots[0], U, rtol=0, atol=1e-15)
    True
    >>> np.allclose(vol.knots[1], V, rtol=0, atol=1e-15)
    True
    >>> np.allclose(vol.knots[2], W, rtol=0, atol=1e-15)
    True
    >>> np.allclose(vol.points,   C, rtol=0, atol=1e-15)
    True
    >>> np.allclose(vol.weights,  1, rtol=0, atol=1e-15)
    True

    """

    def __init__(self, control, knots, weights=None):
        """
        Create a NURBS object.
        """
        #
        if (isinstance(control, tuple) and len(control) == 2):
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
        knots = tuple(np.asarray(k, dtype='d') for k in knots)
        assert len(knots) >= 1
        assert len(knots) <= 3
        for k in knots:
            assert k.ndim == 1
            assert k.size >= 4
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
    def dim(self):
        """
        Parametric dimension of the NURBS object {1,2,3}.
        """
        return self.control.ndim-1

    @property
    def shape(self):
        """
        Shape of the control point net. Also the number of
        basis functions in each parametric direction.
        """
        return self.control.shape[:-1]

    @property
    def degree(self):
        """
        Polynomial degrees for each parametric direction.
        """
        N = [n-1 for n in self.shape]
        M = [len(k)-1 for k in self.knots]
        return tuple(m-n-1 for (m, n) in zip(M, N))

    @property
    def knots(self):
        """
        Knot vectors for each parametric dimension.
        """
        return self._knots

    @property
    def control(self):
        """
        Control point grid in 4D space.
        """
        return self._cntrl

    @property
    def points(self):
        """
        Control point grid projected into 3D space.
        """
        Cw = self.control
        w = self.weights
        return Cw[...,:-1] / w[...,np.newaxis]

    @property
    def weights(self):
        """
        Rational weight grid.
        """
        return self.control[...,-1]

    #

    def copy(self):
        """
        Copy a NURBS object.

        Returns a new instace of the NURBS objects with copies of the
        control points and knot vectors. Modifying the knot vector or
        control points of the returned object WILL NOT affect this
        object.

        Examples
        --------

        Create a random curve, copy the curve, change the control points,
        demonstrate that now c1 and c2 are different

        >>> c1 = NURBS(np.random.rand(5,2),[[0,0,1,2,3,4,4]])
        >>> c2 = c1.copy()
        >>> c2.control[2,:] = [1.0,1.0,0.0,1.0]
        >>> (abs(c2.control-c1.control)).max() < 1.0e-15
        False

        """
        nrb = NURBS.__new__(type(self))
        nrb._cntrl = self.control.copy()
        nrb._knots = tuple(k.copy() for k in self.knots)
        return nrb

    def clone(self):
        """
        Clone a NURBS object.

        Returns a new instace of the NURBS objects with references to
        the control points and knot vectors of this NURBS
        object. Modifying the knot vector or control points of the
        returned object WILL affect this object.

        Examples
        --------

        Create a random curve, copy the curve, change the control points,
        demonstrate that changing c2 changes c1

        >>> c1 = NURBS(np.random.rand(5,2),[[0,0,1,2,3,4,4]])
        >>> c2 = c1.clone()
        >>> c2.control[2,:] = [1.0,1.0,0.0,1.0]
        >>> (abs(c2.control-c1.control)).max() < 1.0e-15
        True

        """
        nrb = NURBS.__new__(type(self))
        nrb._cntrl = self.control
        nrb._knots = self.knots
        return nrb

    #

    def transform(self, trans):
        """
        Apply a scaling, rotation, or a translation to a NURBS object.

        A NURBS object can be scaled, rotated, or translated by
        applying the tranformation to the control points. To contruct
        composite transformations, consult the docstrings in
        transform.py.

        Parameters
        ----------
        trans : array_like
              a matrix or transformation which scales, rotates, and/or
              translates a NURBS object.

        """
        if not isinstance(trans, transform):
            trans = transform(trans)
        Cw = trans(self.control)
        self._cntrl = np.ascontiguousarray(Cw)
        return self

    def translate(self, displ, axis=None):
        t = transform().translate(displ, axis)
        return self.transform(t)

    def move(self, displ, axis=None):
        t = transform().move(displ, axis)
        return self.transform(t)

    def scale(self, scale, axis=None):
        t = transform().scale(scale, axis)
        return self.transform(t)

    def rotate(self, angle, axis=2):
        t = transform().rotate(angle, axis)
        return self.transform(t)

    #

    def transpose(self, axes=None):
        """
        Permute the axes of a NURBS object.

        Permute parametric axes with the given ordering and adjust the
        control points accordingly.

        Parameters
        ----------
        axes : sequence of ints, optional
           By default, reverse order of the axes, otherwise permute
           the axes according to the values given.

        Examples
        --------

        Create a random B-spline volume, swap the 2nd and 3rd axes,
        evaluate and check.

        >>> C = np.random.rand(4,3,2,3)
        >>> U = [0,0,0,0,1,1,1,1]
        >>> V = [0,0,0,1,1,1]
        >>> W = [0,0,1,1]
        >>> vol1 = NURBS(C, [U,V,W])
        >>> vol2 = vol1.clone().transpose([0,2,1])
        >>> u = 0.25; v = 0.50; w = 0.75
        >>> xyz1 = vol1.evaluate(u,v,w)
        >>> xyz2 = vol2.evaluate(u,w,v)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=1e-15)
        True
        >>> vol3 = vol1.clone().transpose()
        >>> xyz3 = vol3.evaluate(w,v,u)
        >>> np.allclose(xyz1, xyz3, rtol=0, atol=1e-15)
        True

        """
        if axes is None:
            axes = range(self.dim)[::-1]
            axes = tuple(axes)
        else:
            axes = tuple(axes)
            assert len(axes) == self.dim
        kaxes = list(axes)
        caxes = kaxes+[self.dim]
        control = self.control.transpose(caxes).copy()
        knots = [self.knots[i] for i in kaxes]
        #
        self._cntrl = np.ascontiguousarray(control)
        self._knots = tuple(knots)
        return self

    def swap(self, axis1, axis2):
        """
        Interchange two parametric axes of NURBS object.

        Parameters
        ----------
        axis1 : int
            First axis.
        axis2 : int
            Second axis.

        Examples
        --------

        Create a random B-spline volume, swap the 2nd and 3rd axes,
        swap the first and last axes, evaluate and check.

        >>> C = np.random.rand(4,3,2,3)
        >>> U = [0,0,0,0,1,1,1,1]
        >>> V = [0,0,0,1,1,1]
        >>> W = [0,0,1,1]
        >>> vol1 = NURBS(C, [U,V,W])
        >>> vol2 = vol1.clone().swap(1,2)
        >>> vol3 = vol1.clone().swap(0,-1)
        >>> u = 0.25; v = 0.50; w = 0.75
        >>> xyz1 = vol1.evaluate(u,v,w)
        >>> xyz2 = vol2.evaluate(u,w,v)
        >>> xyz3 = vol3.evaluate(w,v,u)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=1e-15)
        True
        >>> np.allclose(xyz1, xyz3, rtol=0, atol=1e-15)
        True

        """
        allaxes = range(self.dim)
        ax1, ax2 = allaxes[axis1], allaxes[axis2]
        control = self.control.swapaxes(ax1, ax2).copy()
        knots = list(self.knots)
        knots[ax1], knots[ax2] = knots[ax2], knots[ax1]
        #
        self._cntrl = np.ascontiguousarray(control)
        self._knots = tuple(knots)
        return self

    def reverse(self, *axes):
        """
        Reverse the parametric orientation of a NURBS object along a
        specified axis.

        Given an axis or axes, reverse the parametric orientation of
        the NURBS object. If no axes are given, reverses the
        parametric direction of all axes.

        Parameters
        ----------
        axes : int
              axis indices to reverse separated by commas

        Examples
        --------

        Create a curve, copy it, reverse the copy, evaluate at
        equivalent parametric points, verify the point is the same.

        >>> c1 = NURBS(np.random.rand(6,2),[[0,0,0,0.25,0.75,0.75,1,1,1]])
        >>> c2 = c1.copy()
        >>> c2 = c2.reverse()
        >>> u = 0.3
        >>> (abs(c1.evaluate(u)-c2.evaluate(1.0-u))).max() < 1.0e-15
        True

        """
        def CntRev(C, axis):
            dim = C.ndim - 1
            index = [slice(None, None, None)] * dim
            index[axis] = slice(None, None, -1)
            return C[tuple(index)]
        def KntRev(p, U):
            U0, U1 = U[p], U[-p-1]
            return (U1+U0)-U[::-1]
        allaxes = range(self.dim)
        if not axes: axes = allaxes
        control = self.control
        knots = list(self.knots)
        degree = self.degree
        for axis in axes:
            axis = allaxes[axis]
            control = CntRev(control, axis)
            knots[axis] = KntRev(degree[axis], knots[axis])
        #
        self._cntrl = np.ascontiguousarray(control)
        self._knots = tuple(knots)
        return self

    #

    def remap(self, axis, start, end):
        """
        Linear reparametrization of knot vectors.

        Parameters
        ----------
        axis  : int
        start : float or None
        end   : float or None


        Examples
        --------
        
        >>> c0 = NURBS(np.random.rand(6,3),[[0,0,0,0.25,0.75,0.75,1,1,1]])
        >>> v0 = c0.evaluate([0,0.5,1])
        >>> c1 = c0.copy().remap(0, -2, 2)
        >>> c1.knots[0].tolist()
        [-2.0, -2.0, -2.0, -1.0, 1.0, 1.0, 2.0, 2.0, 2.0]
        >>> v1 = c1.evaluate([-2,0.0,2])
        >>> np.allclose(v0, v1, rtol=0, atol=1e-15)
        True
        >>> c2 = c0.copy().remap(0, None, 2)
        >>> c2.knots[0].tolist()
        [0.0, 0.0, 0.0, 0.5, 1.5, 1.5, 2.0, 2.0, 2.0]
        >>> v2 = c2.evaluate([0,1.0,2])
        >>> np.allclose(v0, v2, rtol=0, atol=1e-15)
        True
        >>> c3 = c0.copy().remap(0, -1, None)
        >>> c3.knots[0].tolist()
        [-1.0, -1.0, -1.0, -0.5, 0.5, 0.5, 1.0, 1.0, 1.0]
        >>> v3 = c3.evaluate([-1,0.0,1])
        >>> np.allclose(v0, v3, rtol=0, atol=1e-15)
        True

        """
        axis = range(self.dim)[axis]
        p = self.degree[axis]
        U = self.knots[axis]
        a, b = U[p], U[-p-1]
        c, d = start, end
        if c is None: c = a
        if d is None: d = b
        c, d = float(c), float(d)
        if (a == c) and (b == d): return self
        #
        assert c < d
        U = U.copy()
        U -= a
        U *= (d-c)/(b-a)
        U += c
        knots = list(self.knots)
        knots[axis] = U
        #
        self._knots = tuple(knots) 
        return self

    def insert(self, axis, value, times=1):
        """
        Insert a single knot value multiple times.

        Parameters
        ----------
        axis : int
        value: float
        times: int

        Examples
        --------

        Create a random curve, insert knots, check error:

        >>> C = np.random.rand(5,3)
        >>> U = [0,0,0,0,0.5,1,1,1,1]
        >>> c1 = NURBS(C, [U])
        >>> c2 = c1.clone().insert(0, 0.25)
        >>> c3 = c2.clone().insert(0, 0.50, 2)
        >>> c4 = c3.clone().insert(0, 0.75, 3)
        >>> u = np.linspace(0,1,100)
        >>> xyz1 = c1.evaluate(u)
        >>> xyz2 = c2.evaluate(u)
        >>> xyz3 = c3.evaluate(u)
        >>> xyz4 = c4.evaluate(u)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=1e-15)
        True
        >>> np.allclose(xyz1, xyz3, rtol=0, atol=1e-15)
        True
        >>> np.allclose(xyz1, xyz4, rtol=0, atol=1e-15)
        True

        Create a random surface, insert knots, check error:

        >>> C = np.random.rand(4,3,3)
        >>> U = [0,0,0,0,1,1,1,1]; V = [0,0,0,1,1,1]
        >>> s1 = NURBS(C, [U,V])
        >>> s1.shape
        (4, 3)
        >>> s2 = s1.clone().insert(0, 0.25).insert(1, 0.75, 2)
        >>> s2.shape
        (5, 5)
        >>> u = v = np.linspace(0,1,100)
        >>> xyz1 = s1.evaluate(u, v)
        >>> xyz2 = s2.evaluate(u, v)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=1e-15)
        True

        """
        axis = range(self.dim)[axis]
        control = self.control.view()
        knots = list(self.knots)
        p = self.degree[axis]
        U = knots[axis]
        assert U[p] <= value <= U[-p-1]
        #
        mult = _api[0].FindMult(p, U, value)
        mult = min(mult, p)
        if times is None: times = p-mult
        assert times + mult <= p
        if times == 0: return self
        #
        InsertKnot = _api[0].InsertKnot
        Pw = np.rollaxis(control, axis, 0)
        shape = Pw.shape
        Pw = Pw.reshape((shape[0], -1))
        V, Qw = InsertKnot(p, U, Pw, value, times)
        Qw.shape = (Qw.shape[0], ) + shape[1:]
        control = np.rollaxis(Qw, 0, axis+1)
        knots[axis] = V
        #
        self._cntrl = np.ascontiguousarray(control)
        self._knots = tuple(knots)
        return self

    def remove(self, axis, value, times=1, deviation=1.0e-9):
        r"""
        Remove a single knot value multiple times.

        Parameters
        ----------
        axis : int
        value: float
        times: int
        deviation: float

        Examples
        --------

        Create a random curve, insert knots,
        remove knots, check error:

        >>> C = np.random.rand(5,3)
        >>> U = [0,0,0,0,0.5,1,1,1,1]
        >>> c1 = NURBS(C, [U])
        >>> c1.shape
        (5,)
        >>> c2 = c1.clone().insert(0, 0.25) \
        ...                .remove(0, 0.25) \
        ...                .remove(0, 0.25)
        >>> c2.shape
        (5,)
        >>> c3 = c2.clone().insert(0, 0.50, 2) \
        ...                .remove(0, 0.50, 2)
        >>> c3.shape
        (5,)
        >>> c4 = c3.clone().insert(0, 0.75, 3) \
        ...                .remove(0, 0.75, 1) \
        ...                .remove(0, 0.75, 2) \
        >>> c3.shape
        (5,)
        >>> u = np.linspace(0,1,100)
        >>> xyz1 = c1.evaluate(u)
        >>> xyz2 = c2.evaluate(u)
        >>> xyz3 = c3.evaluate(u)
        >>> xyz4 = c4.evaluate(u)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=1e-15)
        True
        >>> np.allclose(xyz1, xyz3, rtol=0, atol=1e-15)
        True
        >>> np.allclose(xyz1, xyz4, rtol=0, atol=1e-15)
        True

        >>> c0 = c1.remove(0, 0).remove(0, 1)
        >>> np.allclose(c1.control, c1.control, rtol=0, atol=1e-15)
        True

        """
        axis = range(self.dim)[axis]
        control = self.control.view()
        knots = list(self.knots)
        p = self.degree[axis]
        U = knots[axis]
        assert U[p] <= value <= U[-p-1]
        #
        mult = _api[0].FindMult(p, U, value)
        if times is None: times = mult
        times = min(times, mult, p)
        assert times >= 0
        if times == 0: return self
        if value <= U[p]: return self
        if value >= U[-p-1]: return self
        #
        RemoveKnot = _api[0].RemoveKnot
        Pw = np.rollaxis(control, axis, 0)
        shape = Pw.shape
        Pw = Pw.reshape((shape[0], -1))
        wmin = Pw[:,-1].min()
        Pmax = []
        for i in range(Pw.shape[0]):
            Pmax.append(np.linalg.norm(Pw[i,:-1]))
        Pmax = max(Pmax)
        TOL = deviation*wmin/(1.+Pmax)
        t, V, Qw = RemoveKnot(p, U, Pw, value, TOL, times)
        if t > 0: V = V[:-t].copy()
        if t > 0: Qw = Qw[:-t,:].copy()
        Qw.shape = (Qw.shape[0], ) + shape[1:]
        control = np.rollaxis(Qw, 0, axis+1)
        knots[axis] = V
        #
        self._cntrl = np.ascontiguousarray(control)
        self._knots = tuple(knots)
        return self

    def clamp(self, *axes):
        """

        >>> C = np.random.rand(3,3)
        >>> U = [0,0,0,1,1,1]
        >>> c1 = NURBS(C, [U])
        >>> c1.knots[0].tolist()
        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        >>> c2 = c1.clone().unclamp()
        >>> c2.knots[0].tolist()
        [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
        >>> c3 = c2.clone().clamp()
        >>> c3.knots[0].tolist()
        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

        >>> C = np.random.rand(4,4)
        >>> U = [0,0,0,0.5,1,1,1]
        >>> c1 = NURBS(C, [U])
        >>> c1.knots[0].tolist()
        [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]
        >>> c2 = c1.clone().unclamp()
        >>> c2.knots[0].tolist()
        [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
        >>> c3 = c2.clone().clamp()
        >>> c3.knots[0].tolist()
        [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]

        """
        allaxes = range(self.dim)
        if not axes: axes = allaxes
        #
        Clamp = _api[0].Clamp
        control = self.control
        knots = list(self.knots)
        degree = self.degree
        for axis in axes:
            axis = allaxes[axis]
            p = degree[axis]
            U = knots[axis]
            Pw = np.rollaxis(control, axis, 0).copy()
            shape = Pw.shape
            Pw.shape = (shape[0], -1)
            V, Qw = Clamp(p, U, Pw)
            Qw.shape = shape
            control = np.rollaxis(Qw, 0, axis+1)
            knots[axis] = V
        #
        self._cntrl = np.ascontiguousarray(control)
        self._knots = tuple(knots)
        return self

    def unclamp(self, *axes):
        """

        Examples
        --------

        Create a random curve, unclamp, check error:

        >>> C = np.random.rand(3,3)
        >>> U = [0,0,0,1,1,1]
        >>> c1 = NURBS(C, [U])
        >>> c1.knots[0].tolist()
        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        >>> c2 = c1.clone().unclamp()
        >>> c2.knots[0].tolist()
        [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
        >>> u = np.linspace(0,1,100)
        >>> xyz1 = c1.evaluate(u)
        >>> xyz2 = c2.evaluate(u)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=1e-15)
        True

        """
        allaxes = range(self.dim)
        if not axes: axes = allaxes
        #
        Unclamp = _api[0].Unclamp
        control = self.control
        knots = list(self.knots)
        degree = self.degree
        #
        for axis in axes:
            axis = allaxes[axis]
            p = self.degree[axis]
            U = knots[axis]
            Pw = np.rollaxis(control, axis, 0).copy()
            shape = Pw.shape
            Pw.shape = (shape[0], -1)
            V, Qw = Unclamp(p, U, Pw)
            Qw.shape = shape
            control = np.rollaxis(Qw, 0, axis+1)
            knots[axis] = V
        #
        self._cntrl = np.ascontiguousarray(control)
        self._knots = tuple(knots)
        return self

    def refine(self, u, *vw, **kargs):
        """
        Knot refine a NURBS object.

        Given a list of knots to insert in each parameter direction,
        refine the curve by knot refinement. The routine both refines
        the NURBS object in place and returns the object.

        Parameters
        ----------
        u, v, w : float or array_like or None
            Knots to insert in each parameter direction

        Examples
        --------

        Create a random surface, knot refine, check error:

        >>> C = np.random.rand(4,3,3)
        >>> U = [0,0,0,0,1,1,1,1]; V = [0,0,0,1,1,1]
        >>> s1 = NURBS(C, [U,V])
        >>> s1.shape
        (4, 3)
        >>> u = [0.25, 0.50, 0.75, 0.75]
        >>> v = [0.33, 0.33, 0.67, 0.67]
        >>> s2 = s1.clone().refine(u, v)
        >>> s3 = s1.clone().refine(u, axis=0).refine(u, axis=1)
        >>> s2.shape
        (8, 7)
        >>> s3.shape
        (8, 7)
        >>> u = v = np.linspace(0,1,100)
        >>> xyz1 = s1.evaluate(u, v)
        >>> xyz2 = s2.evaluate(u, v)
        >>> xyz3 = s3.evaluate(u, v)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=1e-15)
        True
        >>> np.allclose(xyz1, xyz3, rtol=0, atol=1e-15)
        True

        """
        axis = kargs.pop('axis', None)
        if axis is None:
            uvw = (u,) + vw
            assert len(uvw) == self.dim
            for ax, u in enumerate(uvw):
                self.refine(u, axis=ax)
            return self
        #
        def Arg(p, U, u):
            u.sort(kind='heapsort')
            assert u[ 0] >= U[p]
            assert u[-1] <= U[-p-1]
            tmp = np.concatenate((u, U[1:-1]))
            try:  np_unique = np.unique1d
            except AttributeError: np_unique = np.unique
            uu, i = np_unique(tmp, return_inverse=True)
            s = np.bincount(i)
            assert s.max() <= p
            return u
        #
        assert len(vw) == 0
        axis = range(self.dim)[axis]
        if u is None: return self
        u = np.asarray(u, dtype='d').ravel()
        if u.size == 0: return self
        p = self.degree[axis]
        U = self.knots[axis]
        u = Arg(p, U, u)
        control = self.control.view()
        knots = list(self.knots)
        #
        RefineKnotVector = _api[0].RefineKnotVector
        Pw = np.rollaxis(control, axis, 0)
        shape = Pw.shape
        Pw = Pw.reshape((shape[0], -1))
        V, Qw = RefineKnotVector(p, U, Pw, u)
        Qw.shape = Qw.shape[:1] + shape[1:]
        control = np.rollaxis(Qw, 0, axis+1)
        knots[axis] = V
        #
        self._cntrl = np.ascontiguousarray(control)
        self._knots = tuple(knots)
        return self

    def elevate(self, t, *sr, **kargs):
        """
        Degree elevate a NURBS object.

        Given a list of polynomial degrees to elevate in each
        parameter direction, refine the curve. The routine both
        refines the NURBS object in place and returns the object.

        Parameters
        ----------
        t, s, r : int or None
            Polynomial orders to elevate by in each parametric direction.

        Examples
        --------

        Create a random curve, degree elevate, check error:

        >>> C = np.random.rand(3,3)
        >>> U = [0,0,0,1,1,1]
        >>> c1 = NURBS(C, [U])
        >>> c1.degree
        (2,)
        >>> c2 = c1.clone().elevate(2)
        >>> c2.degree
        (4,)
        >>> u = np.linspace(0,1,100)
        >>> xyz1 = c1.evaluate(u)
        >>> xyz2 = c2.evaluate(u)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=1e-15)
        True

        Create a random surface, degree elevate, check error:

        >>> C = np.random.rand(3,3,3)
        >>> U = [0,0,0,1,1,1]; V = [0,0,0.5,1,1]
        >>> s1 = NURBS(C, [U,V])
        >>> s1.degree
        (2, 1)
        >>> s2 = s1.clone().elevate(1, axis=0).elevate(1, axis=1)
        >>> s2.degree
        (3, 2)
        >>> u = v = np.linspace(0,1,100)
        >>> xyz1 = s1.evaluate(u, v)
        >>> xyz2 = s2.evaluate(u, v)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=1e-15)
        True

        """
        axis = kargs.pop('axis', None)
        if axis is None:
            tsr = (t,) + sr
            assert len(tsr) == self.dim
            for ax, t in enumerate(tsr):
                self.elevate(t, axis=ax)
            return self
        #
        assert len(sr) == 0
        axis = range(self.dim)[axis]
        if t is None: return self
        t = int(t)
        assert t >= 0
        if t == 0: return self
        p = self.degree[axis]
        U = self.knots[axis]
        control = self.control.view()
        knots = list(self.knots)
        #
        DegreeElevate = _api[0].DegreeElevate
        Pw = np.rollaxis(control, axis, 0)
        shape = Pw.shape
        Pw = Pw.reshape((shape[0], -1))
        V, Qw = DegreeElevate(p, U, Pw, t)
        Qw.shape = Qw.shape[:1] + shape[1:]
        control = np.rollaxis(Qw, 0, axis+1)
        knots[axis] = V
        #
        self._cntrl = np.ascontiguousarray(control)
        self._knots = tuple(knots)
        return self

    #

    def slice(self, axis, start, end):
        """
        Parameters
        ----------
        axis: int
        start, end: float

        Examples
        --------

        Create a random curve, slice, check error:

        >>> C = np.random.rand(5,3)
        >>> U = [0,0,0,0,0.5,1,1,1,1]
        >>> crv = NURBS(C,[U])
        >>> sub = crv.slice(0,0.5,0.75)
        >>> u = np.linspace(0.5,0.75,100)
        >>> xyz1 = crv.evaluate(u)
        >>> xyz2 = sub.evaluate(u)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=1e-15)
        True

        Create a random surface, slice along first axis,
        check error:

        >>> C = np.random.rand(5,4,3)
        >>> U = [0,0,0,0,0.5,1,1,1,1]; V = U[1:-1];
        >>> srf = NURBS(C,[U,V])
        >>> sub = srf.slice(0,1./3,2./3)
        >>> u = np.linspace(1./3,2./3,100)
        >>> v = np.linspace(0,1,100)
        >>> xyz1 = srf.evaluate(u,v)
        >>> xyz2 = sub.evaluate(u,v)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=2e-15)
        True

        Create a random volume, slice along first and second axes,
        check error:

        >>> C = np.random.rand(4,3,2,3)
        >>> U = [0,0,0,0,1,1,1,1]; V = U[1:-1]; W = V[1:-1];
        >>> vol = NURBS(C,[U,V,W])
        >>> sub = vol.slice(0,1./3,2./3).slice(1,0.25,0.75)
        >>> sub = sub.slice(0,1./3,2./3) # no-op
        >>> sub = sub.slice(1,0.25,0.75) # no-op
        >>> sub = sub.slice(2,0,1)       # no-op
        >>> u = np.linspace(1./3,2./3,100)
        >>> v = np.linspace(0.25,0.75,100)
        >>> w = np.linspace(0,1,100)
        >>> xyz1 = vol.evaluate(u,v,w)
        >>> xyz2 = sub.evaluate(u,v,w)
        >>> np.allclose(xyz1, xyz2, rtol=0, atol=2e-15)
        True

        """
        dim = self.dim
        axis = range(dim)[axis]
        p = self.degree[axis]
        U = self.knots[axis]
        if start is None:
            start = U[p]
        if end is None:
            end = U[-p-1]
        start = float(start)
        end = float(end)
        assert start >= U[p]
        assert end <= U[-p-1]
        assert start < end
        #
        FindSpanMult = _api[0].FindSpanMult
        u0 = start
        k0, s0 = FindSpanMult(p,U,u0)
        t0 = max(0, p-s0)
        u1 = end
        k1, s1 = FindSpanMult(p,U,u1)
        t1 = max(0, p-s1)
        u = np.repeat([u0,u1],[t0,t1]).astype('d')
        #
        nrb = self.clone()
        uvw = [None] * dim
        uvw[axis] = u
        nrb.refine(*uvw)
        control = nrb.control
        knots = nrb.knots
        #
        U = knots[axis]
        i0 = U.searchsorted(u0, 'r') - 1
        i1 = U.searchsorted(u1, 'l')
        index = [slice(None)] * dim
        index[axis] = slice(i0-p, i1)
        control = control[index].copy()
        Ul = U[i0].repeat(p)
        Uc = U[i0:i1+1]
        Ur = U[i1].repeat(p)
        knots = list(nrb.knots)
        knots[axis] = np.concatenate([Ul, Uc, Ur])
        #
        nrb = NURBS.__new__(type(self))
        nrb._cntrl = np.ascontiguousarray(control)
        nrb._knots = tuple(knots)
        return nrb

    def extract(self, axis, value):
        """
        Extract lower dimensional NURBS object.

        Examples
        --------

        Create a random volume, extract a surface along the 3rd
        parametric direction at w=0.3, further extract a curve from
        the surface along the 1st parametric direction at u=0.5,
        compare evaluations at equivalent points.

        >>> C = np.random.rand(4,3,2,3)
        >>> U = [0,0,0,0,1,1,1,1]; V = U[1:-1]; W = V[1:-1];
        >>> u = 0.5; v = 0.75; w = 0.3;
        >>> vol = NURBS(C,[U,V,W])
        >>> vol.shape
        (4, 3, 2)
        >>> vol.degree
        (3, 2, 1)
        >>> srf = vol.extract(2,w)
        >>> srf.shape
        (4, 3)
        >>> srf.degree
        (3, 2)
        >>> crv = srf.extract(0,u)
        >>> crv.shape
        (3,)
        >>> crv.degree
        (2,)
        >>> p1 = vol.evaluate(u,v,w)
        >>> p2 = srf.evaluate(u,v)
        >>> p3 = crv.evaluate(v)
        >>> np.allclose(p1, p2, rtol=0, atol=1e-15)
        True
        >>> np.allclose(p2, p3, rtol=0, atol=1e-15)
        True

        """
        assert self.dim > 1
        axis = range(self.dim)[axis]
        p = self.degree[axis]
        U = self.knots[axis]
        u = float(value)
        assert u >= U[p]
        assert u <= U[-p-1]
        control = self.control.view()
        knots = list(self.knots)
        #
        Extract = _api[0].Extract
        Pw = np.rollaxis(control, axis, 0)
        shape = Pw.shape
        Pw = Pw.reshape((shape[0], -1))
        Cw = Extract(p,U,Pw,u)
        control = Cw.reshape(shape[1:])
        del knots[axis]
        #
        nrb = NURBS.__new__(NURBS)
        nrb._cntrl = np.ascontiguousarray(control)
        nrb._knots = tuple(knots)
        return nrb

    def boundary(self, axis, side):
        """
        Extract the boundary of a NURBS object along the specified
        axis and side.

        Parameters
        ----------
        axis : int
              index of axis along which to extract boundary
        side : int
              side of axis from which to extract the boundary

        """
        assert self.dim > 1
        assert side in (0, 1)
        axis = range(self.dim)[axis]
        p = self.degree[axis]
        U = self.knots[axis]
        if side == 0:
            value = U[p]
        else:
            value = U[-p-1]
        return self.extract(axis, value)

    #

    def evaluate(self, u, *vw):
        """
        Evaluate the NURBS object at the given parametric values.

        Parameters
        ----------
        u, v, w : float or array_like

        Examples
        --------

        >>> C = [[-1,0],[0,1],[1,0]]
        >>> U = [0,0,0,1,1,1]
        >>> crv = NURBS(C, [U])
        >>> crv.evaluate(0.5).tolist()
        [0.0, 0.5, 0.0]
        >>> crv.evaluate([0.5]).tolist()
        [[0.0, 0.5, 0.0]]
        >>> crv.evaluate([0, 1]).tolist()
        [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]

        """
        uvw = (u,) + vw
        assert len(uvw) == self.dim
        uvw = [np.asarray(a, dtype='d') for a in uvw]
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
        #
        newshape = list(C.shape)
        remove = [i for (i, a) in enumerate(uvw) if not a.ndim]
        for i in reversed(remove): del newshape[i]
        C.shape = newshape
        return C

    #

    def plot(self):
        pass
