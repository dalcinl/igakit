import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.cm as _cm
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
    if m == 'cone':
        kwargs['marker'] = '^'
    ax = mpl.gca()
    x, y, z = [c.ravel() for c in (x, y, z)]
    pts = ax.scatter(x, y, z, s=30, **kwargs)

def line3d(lines, **kwargs):
    r = kwargs.pop('representation', None)
    t = kwargs.pop('tube_radius', None)
    ax = mpl.gca()
    for (x, y ,z) in lines:
        x, y, z = x.T, y.T, z.T
        if x.ndim == 1:
            options = dict(kwargs)
            x, y, z = (c.ravel() for c in (x, y, z))
            lns = ax.plot(x, y, z, **options)
        if x.ndim == 2:
            options = dict(kwargs)
            options['rstride'] = 1
            options['cstride'] = 1
            lns = ax.plot_wireframe(x, y, z, **options)
        if x.ndim == 3:
            options = dict(kwargs)
            options['rstride'] = 1
            options['cstride'] = 1
            I, J, K = [np.arange(s) for s in x.shape]
            for i in I:
                X, Y, Z = (c[i,...] for c in (x, y, z))
                lns = ax.plot_wireframe(X, Y, Z, **options)
            options = dict(kwargs)
            for j in J:
                for k in K:
                    X, Y, Z = (c[:,j,k].ravel() for c in (x, y, z))
                    lns = ax.plot(X, Y, Z, **options)

def surf3d(surfs=(), **kwargs):
    r = kwargs.pop('representation', None)
    t = kwargs.pop('tube_radius', None)
    ax = mpl.gca()
    for (x, y ,z) in surfs:
        x, y, z = x.T, y.T, z.T
        if x.ndim == 2:
            options = dict(kwargs)
            options['rstride'] = 1
            options['cstride'] = 1
            if 'linewidth' not in options:
                options['linewidth'] = 0
            if ('color' not in options and
                'cmap'  not in options):
                options['cmap'] = _cm.jet
            options['shade'] = True
            options['antialiased'] = True
            srf = ax.plot_surface(x, y, z, **options)
        if x.ndim == 3:
            continue # XXX
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
