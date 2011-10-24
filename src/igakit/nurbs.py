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
    curves/surfaces/volumes by specifying a control point array and a
    list of knot vectors.

    Attributes
    ----------
    control : numpy.ndarray
          control point array in 4D space (includes rational weights)
    degree : array of int
          array of polynomial degrees
    dim : int
          parametric dimension of NURBS object {1,2,3}
    knots : array of numpy.ndarray
          array of knot vectors, one for each parametric dim
    points : numpy.ndarray
         control point array projected into 3D space
    shape : array of int
          array of number of control points in each parametric dimension

    """
    def __init__(self, control, knots, weights=None):
        """
        Creates a NURBS object
        
        Given a control point array and a list of knot vectors,
        creates a NURBS curve/surface/volume. The parametric dimension
        will be inferred from the number of knot vectors specified. If
        weights are omitted, the object will be polynomial (B-spline).

        Parameters
        ----------
        control : array_like
              two dimensional array of control points. The first 
              dimension is the number of control points. The second 
              dimension is the spatial dimension of the object (could be 
              4D where the last dimension is the weight)
        knots : list of array_like
              list of knot vectors. Number of knot vectors will define 
              what kind of NURBS object is created (1=curve,2=surface,
              3=volume)
        weights :
              optional weights


        Examples
        --------

        Initialize a random B-spline curve

        >>> U = [0,0,0, 1,1,1]
        >>> XYZ = np.random.rand(3,3)
        >>> crv = NURBS(XYZ,[U])

        Initialize a quarter circle NURBS curve and check maximum error at
        1000 points

        >>> XYZ[0,:] = [0.0, 1.0, 0.0]
        >>> XYZ[1,:] = [1.0, 1.0, 0.0]
        >>> XYZ[2,:] = [1.0, 0.0, 0.0]
        >>> W = np.asarray([1.0,0.5*np.sqrt(2.0),1.0])
        >>> crv = NURBS(XYZ,[U],W)
        >>> u = np.linspace(0,1.0,1000,endpoint=True)
        >>> xy = crv.evaluate(u)
        >>> xy[:,2] = abs(1.0-np.sqrt(xy[:,0]**2+xy[:,1]**2))
        >>> xy[:,2].max() < 1.0e-15
        True

        Initialize a quarter circle NURBS curve with projected control points
        and check maximum error of 1000 points

        >>> wgt = 0.5*np.sqrt(2.0)
        >>> XYZW = np.zeros((3,4))
        >>> XYZW[0,:] = [0.0, 1.0, 0.0, 1.0]
        >>> XYZW[1,:] = [wgt, wgt, 0.0, wgt]
        >>> XYZW[2,:] = [1.0, 0.0, 0.0, 1.0]
        >>> crv = NURBS(XYZW,[U])
        >>> u = np.linspace(0,1.0,1000,endpoint=True)
        >>> xy = crv.evaluate(u)
        >>> xy[:,2] = abs(1.0-np.sqrt(xy[:,0]**2+xy[:,1]**2))
        >>> xy[:,2].max() < 1.0e-15
        True

        Initialize a random B-spline surface 

        >>> V = [0,0, 1,1]
        >>> XYZ = np.random.rand(3,2,3)
        >>> srf = NURBS(XYZ,[U,V])

        Initialize a random B-spline volume

        >>> W = [0,0,0,0, 0.25, 0.5,0.5, 1,1,1,1]
        >>> XYZ = np.random.rand(3,2,7,3)
        >>> vol = NURBS(XYZ,[U,V,W])

        """
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
        """
        Copies a NURBS object.

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
        nrb._knots = tuple([k.copy() for k in self.knots])
        return nrb

    def clone(self):
        """
        Clones a NURBS object.

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

    def transpose(self, *axes):
        """
        Swap axes of a NURBS object

        Given the axes to swap (or none if a NURBS surface), this
        swaps the two given parametric axes and adjusts the control
        points accordingly. Modified the object in place and returns a
        reference.

        Parameters
        ----------
        axes : int
             new axes order separated by commas

        Examples
        --------
        
        Create a volume, copy the volume, swap the 2nd and 3rd axes,
        evalute each at symmetric locations and verify the point is
        the same.

        >>> v1 = NURBS(np.random.rand(4,3,2,3),[ [0,0,0,0,1,1,1,1], [0,0,0,1,1,1], [0,0,1,1] ])
        >>> v2 = v1.copy()
        >>> v2 = v2.transpose(0,2,1)
        >>> u = 0.25; v = 0.5; w = 0.75
        >>> (abs(v1.evaluate(u,v,w)-v2.evaluate(u,w,v))).max() < 1.0e-15
        True

        """
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
        """
        Knot refine a NURBS object.

        Given a list of knots to insert in each parameter direction,
        refine the curve by knot refinement. The routine both refines
        the NURBS object in place and returns the object.

        Parameters
        ----------
        uvw : array_like
              a list of knots to insert in each parameter direction

        Examples
        --------

        Create a random surface, copy the surface, knot refine the
        copy, check maximum error at 100 points

        >>> s1 = NURBS(np.random.rand(4,3,3),[ [0,0,0,0,1,1,1,1], [0,0,0,1,1,1] ])
        >>> s2 = s1.copy()
        >>> s2 = s2.refine([0.25, 0.5, 0.75, 0.75], [0.33, 0.33, 0.67, 0.67])
        >>> u = np.linspace(0.0,1.0,10,endpoint=True)
        >>> (abs(s1.evaluate(u,u)-s2.evaluate(u,u))).max() < 1.0e-15
        True

        """
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
        """
        Degree elevate a NURBS object.

        Given a list of polynomial degrees to elevate in each
        parameter direction, refine the curve. The routine both
        refines the NURBS object in place and returns the object.

        Parameters
        ----------
        rst : int
              polynomial orders to elevate by in each parametric direction

        Examples
        --------

        Create a random curve, copy the curve, degree elevate the
        copy, check maximum error at 100 points

        >>> c1 = NURBS(np.random.rand(2,3),[ [0,0,1,1] ])
        >>> c2 = c1.copy()
        >>> c2 = c2.elevate(1)
        >>> u = np.linspace(0.0,1.0,100,endpoint=True)
        >>> (abs(c1.evaluate(u)-c2.evaluate(u))).max() < 1.0e-15
        True

        Create a random surface, copy the surface, degree elevate the
        copy, check maximum error at 10000 points

        >>> s1 = NURBS(np.random.rand(3,3,3),[ [0,0,0,1,1,1], [0,0,0.5,1,1] ])
        >>> s2 = s1.copy()
        >>> s2 = s2.elevate(1,1)
        >>> u = np.linspace(0.0,1.0,100,endpoint=True)
        >>> (abs(s1.evaluate(u,u)-s2.evaluate(u,u))).max() < 1.0e-8
        True

        .. note: Lisandro, this is strange to me that the curves and
        surfaces have pointwise error on the order of 10^-9 after
        degree elevation.

        """
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
