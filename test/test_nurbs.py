import bootstrap
import numpy as np
from igakit.nurbs import NURBS, transform

def make_crv():
    C = [[ 6.0, 0.0, 6.0],
         [-5.5, 0.5, 5.5],
         [-5.0, 1.0,-5.0],
         [ 4.5, 1.5,-4.5],
         [ 4.0, 2.0, 4.0],
         [-3.5, 2.5, 3.5],
         [-3.0, 3.0,-3.0],
         [ 2.5, 3.5,-2.5],
         [ 2.0, 4.0, 2.0],
         [-1.5, 4.5, 1.5],
         [-1.0, 5.0,-1.0],
         [ 0.5, 5.5,-0.5],
         [ 0.0, 6.0, 0.0],]
    U = [0, 0, 0, 0,
         .1, .2, .3, .4, .5, .6, .7, .8, .9,
         1, 1, 1, 1,]
    nrb = NURBS([U], C)
    return nrb

def make_srf():
    C = np.zeros((3,5,5))
    C[:,:,0] = [[ 0.0,  3.0,  5.0,  8.0, 10.0],
                [ 0.0,  0.0,  0.0,  0.0,  0.0],
                [ 2.0,  2.0,  7.0,  7.0,  8.0],]
    C[:,:,1] = [[ 0.0,  3.0,  5.0,  8.0, 10.0],
                [ 3.0,  3.0,  3.0,  3.0,  3.0],
                [ 0.0,  0.0,  5.0,  5.0,  7.0],]
    C[:,:,2] = [[ 0.0,  3.0,  5.0,  8.0, 10.0],
                [ 5.0,  5.0,  5.0,  5.0,  5.0],
                [ 0.0,  0.0,  5.0,  5.0,  7.0],]
    C[:,:,3] = [[ 0.0,  3.0,  5.0,  8.0, 10.0],
                [ 8.0,  8.0,  8.0,  8.0,  8.0],
                [ 5.0,  5.0,  8.0,  8.0, 10.0],]
    C[:,:,4] = [[ 0.0,  3.0,  5.0,  8.0, 10.0],
                [10.0, 10.0, 10.0, 10.0, 10.0],
                [ 5.0,  5.0,  8.0,  8.0, 10.0],]
    C = C.transpose()
    U = [0, 0, 0, 1/3., 2/3., 1, 1, 1,]
    V = [0, 0, 0, 1/3., 2/3., 1, 1, 1,]
    nrb = NURBS([U,V], C)
    return nrb

def make_vol():
    C = np.zeros((3,2,3,4), dtype='d')
    C[0,0,:,0:2] = [0.0, 0.5]
    C[1,0,:,0:2] = [0.5, 0.5]
    C[2,0,:,0:2] = [0.5, 0.0]
    C[0,1,:,0:2] = [0.0, 1.0]
    C[1,1,:,0:2] = [1.0, 1.0]
    C[2,1,:,0:2] = [1.0, 0.0]
    C[:,:,0,2] = 0.0
    C[:,:,1,2] = 0.5
    C[:,:,2,2] = 1.0
    C[:,:,:,3] = 1.0
    C[1,:,:,:] *= np.sqrt(2)/2
    U = [0,0,0,     1,1,1]
    V = [0,0,         1,1]
    W = [0,0,   0.5,  1,1]
    nrb = NURBS([U,V,W],C)
    return nrb


def test_vol_bnd(plot=False):
    nrb = make_vol()
    result =[[None]*2 for i in range(3)]
    for side in (0, 1):
        for axis in range(3):
            bnd = nrb.boundary(axis, side)
            result[axis][side] = bnd
            if plot: bnd.plot()
    if plot:
        plt.figure()
        for side in (0, 1):
            for axis in range(3):
                result[axis][side].plot()
        plt.title('Boundary')
    return result

def test_vol_ext(plot=False):
    nrb = make_vol()
    if plot: plt.figure()
    for axis in range(3):
        for value in [0, 0.25, 0.5, 0.75, 1]:
            ext = nrb.extract(axis, value)
            if plot: ext.plot()
    if plot: plt.title('Extract')

def test_vol_degelv(plot=False):
    nrb = make_vol()
    assert nrb.degree == (2, 1, 1)
    for i, t in enumerate((0, 1, 1)):
        nrb.elevate(i, t)
    assert nrb.degree == (2, 2, 2)
    if plot:
        plt.figure()
        nrb.plot()
        plt.title('Degree Elevation')
    return nrb

def test_vol_kntins(plot=False):
    nrb = make_vol()
    assert nrb.shape == (3, 2, 3)
    nrb.refine(0, [0.25, 0.5, 0.75])
    nrb.refine(1, [0.25, 0.5, 0.75])
    nrb.refine(2, [0.25,      0.75])
    assert nrb.shape == (3+3, 2+3, 3+2)
    if plot:
        plt.figure()
        nrb.plot()
        plt.title('Knot Refinement')
    return nrb

def test_vol_degelv_kntins(plot=False):
    nrb = make_vol()
    for i, t in enumerate((0, 1, 1)):
        nrb.elevate(i, t)
    for i, u in enumerate(([0.25, 0.5, 0.75],
                           [0.25, 0.5, 0.75],
                           [0.25,      0.75])):
        nrb.refine(i, u)
    if plot:
        plt.figure()
        nrb.plot()
        plt.title('Deg Elev + Knt Ref')
    return nrb

def test_evaluate():
    crv = make_crv()
    srf = make_srf()
    vol = make_vol()
    for nrb in (crv, srf, vol):
        C = nrb.evaluate(*nrb.knots)
        C2, D2 = nrb.evaluate(*nrb.knots, **dict(fields=nrb.points))
        C3, D3 = nrb.evaluate(*nrb.knots, **dict(fields=nrb.weights))
        assert np.allclose(C, C2)
        assert np.allclose(C, C3)
        assert np.allclose(C, D2)
        assert D3.shape[:-1] == C3.shape[:-1]

# ---

if __name__ == '__main__':
    try:
        from matplotlib import pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        PLOT=1
    except ImportError:
        PLOT=0
    crv = make_crv()
    srf = make_srf()
    vol = make_vol()
    test_vol_bnd(PLOT)
    test_vol_ext(PLOT)
    test_vol_degelv(PLOT)
    test_vol_kntins(PLOT)
    test_vol_degelv_kntins(PLOT)
    test_evaluate()
    plt.show()
