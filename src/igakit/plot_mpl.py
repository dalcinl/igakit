import numpy as np
import matplotlib.pyplot as mpl
import mpl_toolkits.mplot3d as _mplot3d
import matplotlib.colors as _colors

class colors(object): pass
colors = colors()
colors.__dict__.update(_colors.cnames)

show = mpl.show

def figure(*args, **kwargs):
    fig = mpl.figure(*args)
    ax = fig.add_subplot(111, projection='3d', **kwargs)
    return ax

close = mpl.close
save = mpl.savefig

def title(*args, **kwargs):
    mpl.gca().set_title(*args, **kwargs)

def xlabel(*args, **kwargs):
    mpl.gca().set_xlabel(*args, **kwargs)

def ylabel(*args, **kwargs):
    mpl.gca().set_ylabel(*args, **kwargs)

def zlabel(*args, **kwargs):
    mpl.gca().set_zlabel(*args, **kwargs)

colorbar = mpl.colorbar

def points3d(x, y ,z, **kwargs):
    x, y, z = x.T, y.T, z.T
    m = kwargs.pop('mode', 'sphere')
    if m == 'sphere':
        kwargs['marker'] = 'o'
    if m == 'cube':
        kwargs['marker'] = 's'
    ax = mpl.gca()
    x, y, z = [c.ravel() for c in (x, y, z)]
    pts = ax.scatter(x, y, z, s=30, **kwargs)

def plot3d(x, y ,z, **kwargs):
    x, y, z = x.T, y.T, z.T
    r = kwargs.pop('representation', None)
    t = kwargs.pop('tube_radius', None)
    ax = mpl.gca()
    lns = ax.plot(x, y, z, **kwargs)
    return lns

def mesh(x, y, z, **kwargs):
    x, y, z = x.T, y.T, z.T
    r = kwargs.pop('representation', 'surface')
    kwargs['rstride'] = 1
    kwargs['cstride'] = 1
    ax = mpl.gca()
    if r == 'surface':
        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = 0
        kwargs['antialiased'] = True
        srf = ax.plot_surface(x, y, z, **kwargs)
    elif r == 'wireframe':
        srf = ax.plot_wireframe(x, y, z, **kwargs)
    return srf

def grid3d(x, y ,z, **kwargs):
    x, y, z = x.T, y.T, z.T
    r = kwargs.pop('representation', None)
    t = kwargs.pop('tube_radius', None)
    ax = mpl.gca()
    if x.ndim == 1:
        x, y, z = [c.ravel() for c in (x, y, z)]
        grd = ax.plot(x, y, z, **kwargs)
        return grd
    if x.ndim == 2:
        kwargs['rstride'] = 1
        kwargs['cstride'] = 1
        grd = ax.plot_wireframe(x, y, z, **kwargs)
        return grd
    if x.ndim == 3:
        I, J, K = [np.arange(s) for s in x.shape]
        for j in J:
            for k in K:
                X = x[:,j,k].ravel()
                Y = y[:,j,k].ravel()
                Z = z[:,j,k].ravel()
                ax.plot(X, Y, Z, **kwargs)
        kwargs['rstride'] = 1
        kwargs['cstride'] = 1
        for i in I:
            X = x[i,...]
            Y = y[i,...]
            Z = z[i,...]
            ax.plot_wireframe(X, Y, Z, **kwargs)

_resolution = { 1:128, 2:16 }
