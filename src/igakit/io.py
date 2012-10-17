from igakit.nurbs import NURBS
import numpy as np

__all__ = ['PetIGA']

class PetIGA(object):

    """
    PetIGA_ Reader/Writer

    .. _PetIGA: https://bitbucket.org/dalcinl/petiga

    """

    IGA_ID = 1211299
    VEC_ID = 1211214

    def __init__(self, integer='i4', real='f8', scalar='f8'):
        I, R, S = integer, real, scalar
        I = np.dtype(I).newbyteorder('>')
        R = np.dtype(R).newbyteorder('>')
        S = np.dtype(S).newbyteorder('>')
        self.dtypes = (I, R, S)

    def write(self, filename, nurbs,
              geometry=True, nsd=None):
        """
        Parameters
        ----------
        filename : string
        nurbs : NURBS
        geometry : bool, optional
        nsd : int, optional

        """
        IGA_ID = self.IGA_ID
        VEC_ID = self.VEC_ID
        I, R, S = self.dtypes
        def _write(fh, T, a):
            return np.asarray(a, T).ravel('f').tofile(fh)
        if geometry:
            dim = nurbs.dim
            if nsd is None: nsd = dim
            assert dim <= nsd <= 3
            descr = 1
        else:
            assert nsd is None
            descr = 0
        fh = open(filename, 'wb')
        _write(fh, I, IGA_ID)
        _write(fh, I, descr)
        _write(fh, I, nurbs.dim)
        for p, U in zip(nurbs.degree, nurbs.knots):
            _write(fh, I, p)
            _write(fh, I, len(U))
            _write(fh, R, U)
        if abs(descr) >= 1:
            idx = list(range(nsd))+[3]
            Cw = nurbs.control[..., idx]
            Cw = np.rollaxis(Cw, -1)
            _write(fh, I, nsd)
            _write(fh, I, VEC_ID)
            _write(fh, I, Cw.size)
            _write(fh, S, Cw)
        fh.flush()
        fh.close()

    def read(self, filename):
        """
        Parameters
        ----------
        filename : string

        Returns
        -------
        nurbs : NURBS

        """
        IGA_ID = self.IGA_ID
        VEC_ID = self.VEC_ID
        I, R, S = self.dtypes
        def _read(fh, T, n=None):
            a = np.fromfile(fh, T, n or 1)
            if n is None: a.shape = ()
            return a
        fh = open(filename, 'rb')
        iga_id = _read(fh, I)
        assert iga_id == IGA_ID
        descr = _read(fh, I)
        dim = _read(fh, I)
        knots, sizes = [], []
        for i in range(dim):
            p = _read(fh, I)
            m = _read(fh, I)
            U = _read(fh, R, m)
            knots.append(U)
            sizes.append(m-p-1)
        if abs(descr) >= 1:
            nsd = _read(fh, I)
            assert dim <= nsd <= 3
            vec_id = _read(fh, I)
            assert vec_id == VEC_ID
            n  = _read(fh, I)
            A = _read(fh, S, n)
            shape = [nsd+1] + sizes
            A = A.reshape(shape, order='f')
            A = np.rollaxis(A, 0, A.ndim)
            shape = sizes + [4]
            Cw = np.zeros(shape, dtype=A.dtype)
            Cw[..., :nsd] = A[..., :-1]
            Cw[...,   -1] = A[...,  -1]
        else:
            Cw = None
        nurbs = NURBS(knots, Cw)
        fh.close()
        return nurbs
