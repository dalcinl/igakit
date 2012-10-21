import numpy as np
from igakit.nurbs import NURBS

__all__ = ['PetIGA', 'VTK']

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
            return np.asarray(a, T).tofile(fh)
        #
        dim = nurbs.dim
        knots = nurbs.knots
        degree = nurbs.degree
        if geometry:
            if nsd is None: nsd = dim
            assert dim <= nsd <= 3
            idx = list(range(nsd))+[3]
            Cw = nurbs.control[...,idx]
            Cw = np.rollaxis(Cw, -1).ravel('f')
            descr = 1
        else:
            assert nsd is None
            Cw = None
            descr = 0
        #
        fh = open(filename, 'wb')
        _write(fh, I, IGA_ID)
        _write(fh, I, descr)
        _write(fh, I, nurbs.dim)
        for p, U in zip(degree, knots):
            _write(fh, I, p)
            _write(fh, I, len(U))
            _write(fh, R, U)
        if geometry:
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

        """
        IGA_ID = self.IGA_ID
        VEC_ID = self.VEC_ID
        I, R, S = self.dtypes
        def _read(fh, T, n=None):
            if n is None:
                a = np.fromfile(fh, T, 1)[0]
            else:
                a = np.fromfile(fh, T, n)
            return a.astype(T.newbyteorder('='))
        #
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
            n = _read(fh, I)
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
        fh.close()
        return NURBS(knots, Cw)

    def write_vec(self, filename, array, nurbs=None):
        VEC_ID = self.VEC_ID
        I, R, S = self.dtypes
        def _write(fh, T, a):
            return np.asarray(a, T).tofile(fh)
        #
        A = np.asarray(array)
        if nurbs is not None:
            shape = nurbs.shape + (-1,)
            A = A.reshape(shape)
            A = np.rollaxis(A, -1)
            A = A.ravel('f')
        #
        fh = open(filename, 'wb')
        _write(fh, I, VEC_ID)
        _write(fh, I, A.size)
        _write(fh, S, A)
        fh.flush()
        fh.close()

    def read_vec(self, filename, nurbs=None):
        VEC_ID = self.VEC_ID
        I, R, S = self.dtypes
        def _read(fh, T, n=None):
            if n is None:
                a = np.fromfile(fh, T, 1)[0]
            else:
                a = np.fromfile(fh, T, n)
            return a.astype(T.newbyteorder('='))
        #
        fh = open(filename, 'rb')
        vec_id = _read(fh, I)
        assert vec_id == VEC_ID
        n = _read(fh, I)
        A = _read(fh, S, n)
        fh.close()
        #
        if nurbs is not None:
            shape = (-1,) + nurbs.shape
            A = A.reshape(shape, order='f')
            A = np.rollaxis(A, 0, A.ndim)
            A = A.squeeze()
            A = np.ascontiguousarray(A)
        return A


class VTK(object):

    """
    VTK Writer
    """

    title = 'VTK Data'

    def __init__(self):
        pass

    def write(self, filename, nurbs,
              geometry=True, sampler=None,
              scalars=(), vectors=()):
        """
        Parameters
        ----------
        filename : string
        nurbs : NURBS
        geometry : bool, optional
        sampler : callable, optional
        scalars : dict or sequence of 2-tuple, optional
        vectors : dict or sequence or 2-tuple, optional

        """
        try:  unique = np.unique1d
        except AttributeError: unique = np.unique
        if sampler is None:
            sampler = lambda U: U
        dim  = nurbs.dim
        degs = nurbs.degree
        knts = nurbs.knots
        uvw = [sampler(unique(U[p:-p]))
               for p, U in zip(degs, knts)]
        flag = bool(scalars or vectors)
        out = nurbs.evaluate(*uvw, fields=flag)
        if flag: C, F = out
        else:    C, F = out, out[..., 0:0]

        dimensions = C.shape[:-1] + (1,)*(3-dim)
        coordinates = uvw + [np.zeros(1)]*(3-dim)
        points = np.rollaxis(C, -1).ravel('f')
        points.shape = (-1, 3)
        fields = np.rollaxis(F, -1).ravel('f')
        fields.shape = (len(points), -1)

        if isinstance(scalars, dict):
            keys = sorted(scalars.keys())
            scalars = [(k, scalars[k]) for k in keys]
        else:
            scalars = list(scalars)
        for i, (name, index) in enumerate(scalars):
            array = np.zeros((len(points), 1), dtype='d')
            array[:,0] = fields[:,index]
            scalars[i] = (name, array)

        if isinstance(vectors, dict):
            keys = sorted(vectors.keys())
            vectors = [(k, vectors[k]) for k in keys]
        else:
            vectors = list(vectors)
        for i, (name, index) in enumerate(vectors):
            array = np.zeros((len(points), 3), dtype='d')
            array[:,:len(index)] = fields[:,index]
            vectors[i] = (name, array)

        fh = open(filename, 'wb')

        header = '# vtk DataFile Version %d.%d'
        version = (2, 0)
        fh.write(header % version)
        fh.write('\n')
        title = self.title
        fh.write(title[:255])
        fh.write('\n')

        format = 'BINARY'
        fh.write(format)
        fh.write('\n')

        if geometry:
            dataset_type = 'STRUCTURED_GRID'
            fh.write('DATASET %s' % dataset_type);
            fh.write('\n')
            fh.write('DIMENSIONS %d %d %d' % dimensions)
            fh.write('\n')
            fh.write('POINTS %d %s' % (len(points), 'double'))
            fh.write('\n')
            points.astype('>d').tofile(fh)
            fh.write('\n')
        else:
            dataset_type = 'RECTILINEAR_GRID'
            fh.write('DATASET %s' % dataset_type);
            fh.write('\n')
            fh.write('DIMENSIONS %d %d %d' % dimensions)
            fh.write('\n')
            for X, array in zip("XYZ", coordinates):
                label = X+'_COORDINATES'
                fh.write('%s %s %s' % (label, len(array), 'double'))
                fh.write('\n')
                array.astype('>d').tofile(fh)
                fh.write('\n')

        if (not scalars and
            not vectors):
            fh.flush()
            fh.close()
            return

        data_type = 'POINT_DATA'
        fh.write('%s %d' % (data_type, len(points)))
        fh.write('\n')

        for i, (name, array) in enumerate(scalars):
            attr_type = 'SCALARS'
            attr_name = name or (attr_type.lower() + str(i))
            attr_name = attr_name.replace(' ', '_')
            fh.write('%s %s %s' %(attr_type, attr_name, 'double'))
            fh.write('\n')
            lookup_table = 'default'
            lookup_table = lookup_table.replace(' ', '_')
            fh.write('LOOKUP_TABLE %s' % lookup_table)
            fh.write('\n')
            array.astype('>d').tofile(fh)
            fh.write('\n')

        for i, (name, array) in enumerate(vectors):
            attr_type = 'VECTORS'
            attr_name = name or (attr_type.lower() + str(i))
            attr_name = attr_name.replace(' ', '_')
            fh.write('%s %s %s' %(attr_type, attr_name, 'double'))
            fh.write('\n')
            array.astype('>d').tofile(fh)
            fh.write('\n')

        fh.flush()
        fh.close()
