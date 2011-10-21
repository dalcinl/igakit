import numpy as np

def controlpoints(nurbs, **kwargs):
    plt = get_backend()
    #
    C = nurbs.points
    x, y, z = C.T
    #
    options = dict(kwargs)
    if 'color' not in options:
        options['color'] = (1,0,0)
    if 'mode' not in options:
        options['mode'] = 'sphere'
    #
    pts = plt.points3d(x, y, z, **options)
    if 'scale_factor' not in options:
        try:
            pts.glyph.glyph.scale_factor *= 0.125
        except AttributeError:
            pass
    return pts

def controlgrid(nurbs, **kwargs):
    plt = get_backend()
    #
    C = nurbs.points
    x, y, z = C.T
    #
    options = dict(kwargs)
    if 'color' not in options:
        options['color'] = (0,0,1)
    if 'representation' not in options:
        options['representation'] = 'wireframe'
        options['tube_radius'] = None
    grd = plt.grid3d(x, y, z, **options)
    return grd


def knotpoints(nurbs, **kwargs):
    plt = get_backend()
    #
    uvw = [np.unique1d(U[p:-p]) for (p, U)
           in zip(nurbs.degree, nurbs.knots)]
    C = nurbs.evaluate(*uvw)
    x, y, z = C.T
    #
    options = dict(kwargs)
    if 'color' not in options:
        options['color'] = (0,1,0)
    if 'mode' not in options:
        options['mode'] = 'cube'
    #
    pts = plt.points3d(x, y, z, **options)
    if 'scale_factor' not in options:
        try:
            pts.glyph.glyph.scale_factor *= 0.1
        except AttributeError:
            pass
    return pts

def knotgrid(nurbs, **kwargs):
    plt = get_backend()
    #
    resol = plt._resolution[1]
    uvw = [np.unique1d(U[p:-p]) for (p, U)
           in zip(nurbs.degree, nurbs.knots)]
    lines = []
    for i in range(nurbs.dim):
        u = uvw[i];
        a = np.linspace(u[0], u[-1], resol)
        abc = list(uvw)
        abc[i] = a
        C = nurbs.evaluate(*abc)
        C = np.rollaxis(C, i, -1)
        C = C.reshape((-1, a.size, 3))
        lines.extend(C)
    #
    options = dict(kwargs)
    if 'color' not in options:
        options['color'] = (0,1,0)
    if 'representation' not in options:
        options['representation'] = 'wireframe'
        options['tube_radius'] = None
    #
    lins = []
    for C in lines:
        x, y, z = C.T
        lin = plt.plot3d(x, y, z, **options)
        lins.append(lin)
    return lins

def curve(nurbs, **kwargs):
    plt = get_backend()
    if nurbs.dim == 1:
        resol =plt._resolution[1]
        p, U = nurbs.degree[0], nurbs.knots[0]
        u = np.linspace(U[p], U[-p-1], resol)
        C = nurbs.evaluate(u)
        x, y, z = C.T
        options = dict(kwargs)
        if 'color' not in options:
            options['color'] = (1,1,0)
        if 'representation' not in options:
            options['representation'] = 'wireframe'
        crv = plt.plot3d(x, y, z, **options)
        return crv
    else:
        for axis in range(nurbs.dim):
            for side in range(2):
                bnd = nurbs.boundary(axis, side)
                crv = curve(bnd, **kwargs)
        return None # XXX

def surface(nurbs, **kwargs):
    plt = get_backend()
    if nurbs.dim == 2:
        resol = plt._resolution[2]
        uvw = [np.linspace(U[p], U[-p-1], resol)
               for (p, U) in zip(nurbs.degree, nurbs.knots)]
        C = nurbs.evaluate(*uvw)
        x, y, z = C.T
        options = dict(kwargs)
        if 'color' not in options:
            options['color'] = (1,1,0)
        srf = plt.mesh(x, y, z, **options)
        return srf
    elif nurbs.dim == 3:
        for axis in range(nurbs.dim):
            for side in range(2):
                bnd = nurbs.boundary(axis, side)
                srf = surface(bnd, **kwargs)
        return None # XXX

def volume(nurbs, **kwargs):
    plt = get_backend()
    if nurbs.dim == 3:
        surface(nurbs, **kwargs)
        return None # XXX

# -----

_backend_current = None
_backend_modules = {'mpl': None, 'myv': None}

def use_backend(backend):
    global _backend_current
    try:
        if backend == 'matplotlib': backend = 'mpl'
        if backend == 'mayavi':     backend = 'myv'
        module = _backend_modules[backend]
    except KeyError:
        raise ValueError("unknown backend '%s'" % backend)
    if module is None:
        modulename = 'igakit.plot_' + backend
        module = __import__(modulename, fromlist=[None])
        _backend_modules[backend] = module
    _backend_current = module

def get_backend():
    if _backend_current is None:
        use_backend('mpl')
    return _backend_current

# -----

import sys

class Module(type(sys)):
  def __init__(self, module):
      super(type(self), self).__init__(module.__name__)
      self.__dict__.update(module.__dict__)
      self.module = module
  def __getattr__(self, name):
      module = self.module
      try:
          return getattr(module, name)
      except AttributeError:
          backend = module.get_backend()
          return getattr(backend, name)
sys.modules[__name__] = Module(sys.modules[__name__])

del sys, Module

# -----
