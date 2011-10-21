import bootstrap
import numpy as np
from igakit.igalib import vol

def test_vol_ki(VERB=0, PLOT=0):
    if VERB: print(vol.RefineKnotVector.__doc__)

    px = 2
    py = 1
    pz = 1

    Ux = np.asarray([0,0,0, 1,1,1], dtype=float)
    Uy = np.asarray([0,0,     1,1], dtype=float)
    Uz = np.asarray([0,0,     1,1], dtype=float)

    nx = len(Ux)-1-(px+1)
    ny = len(Uy)-1-(py+1)
    nz = len(Uz)-1-(pz+1)

    Pw = np.zeros((nx+1,ny+1,nz+1,4), dtype=float)

    Pw[0,0,0,:] = [0.0, 0.5, 0.0, 1.0]
    Pw[1,0,0,:] = [0.5, 0.5, 0.0, 1.0]
    Pw[2,0,0,:] = [0.5, 0.0, 0.0, 1.0]
    Pw[0,1,0,:] = [0.0, 1.0, 0.0, 1.0]
    Pw[1,1,0,:] = [1.0, 1.0, 0.0, 1.0]
    Pw[2,1,0,:] = [1.0, 0.0, 0.0, 1.0]

    Pw[0,0,1,:] = [0.0, 0.5, 1.0, 1.0]
    Pw[1,0,1,:] = [0.5, 0.5, 1.0, 1.0]
    Pw[2,0,1,:] = [0.5, 0.0, 1.0, 1.0]
    Pw[0,1,1,:] = [0.0, 1.0, 1.0, 1.0]
    Pw[1,1,1,:] = [1.0, 1.0, 1.0, 1.0]
    Pw[2,1,1,:] = [1.0, 0.0, 1.0, 1.0]

    Pw[1,:,:,:] *= np.sqrt(2)/2

    X = np.asarray([], dtype=float)
    Y = np.asarray([], dtype=float)
    Z = np.asarray([], dtype=float)
    Vx, Vy, Vz, Qw = vol.RefineKnotVector(px,Ux,
                                          py,Uy,
                                          pz,Uz,
                                          Pw,
                                          X,Y,Z)
    assert np.allclose(Ux, Vx)
    assert np.allclose(Uy, Vy)
    assert np.allclose(Uz, Vz)
    assert np.allclose(Pw, Qw)

    X = np.asarray([0.25, 0.50])
    Y = np.asarray([0.50, 0.75])
    Z = np.asarray([0.25, 0.75])
    Vx, Vy, Vz, Qw = vol.RefineKnotVector(px,Ux,
                                          py,Uy,
                                          pz,Uz,
                                          Pw,
                                          X,Y,Z)
    x = np.linspace(Vx[0], Vx[-1], 15)
    y = np.linspace(Vy[0], Vy[-1], 5)
    z = np.linspace(Vy[0], Vy[-1], 4)
    Cw = vol.Evaluate(px,Vx,py,Vy,pz,Vz,Qw,x,y,z)

    Q = Qw[:,:,:,:3] / Qw[:,:,:,3,None]
    C = Cw[:,:,:,:3] / Cw[:,:,:,3,None]

    if not PLOT: return
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.title("Volume - Knot Insertion")

    x = Q[:,:,:,0].flatten()
    y = Q[:,:,:,1].flatten()
    z = Q[:,:,:,2].flatten()
    ax.scatter(x,y,z, c='r')

    x = C[:,:,:,0].flatten()
    y = C[:,:,:,1].flatten()
    z = C[:,:,:,2].flatten()
    ax.scatter(x,y,z, c='b')

    t = np.linspace(0,np.pi/2,100)
    a = np.cos(t)
    b = np.sin(t)
    ax.plot(1.0*a,1.0*b,0,'-k')
    ax.plot(1.0*a,1.0*b,1,'-k')
    ax.plot(0.5*a,0.5*b,0,'-k')
    ax.plot(0.5*a,0.5*b,1,'-k')

    plt.axis("equal")

def test_vol_de(VERB=0, PLOT=0):
    if VERB: print(vol.DegreeElevate.__doc__)

    px = 2
    py = 1
    pz = 1

    Ux = np.asarray([0,0,0, 1,1,1], dtype=float)
    Uy = np.asarray([0,0,     1,1], dtype=float)
    Uz = np.asarray([0,0,     1,1], dtype=float)

    nx = len(Ux)-1-(px+1)
    ny = len(Uy)-1-(py+1)
    nz = len(Uz)-1-(pz+1)

    Pw = np.zeros((nx+1,ny+1,nz+1,4), dtype=float)

    Pw[0,0,0,:] = [0.0, 0.5, 0.0, 1.0]
    Pw[1,0,0,:] = [0.5, 0.5, 0.0, 1.0]
    Pw[2,0,0,:] = [0.5, 0.0, 0.0, 1.0]
    Pw[0,1,0,:] = [0.0, 1.0, 0.0, 1.0]
    Pw[1,1,0,:] = [1.0, 1.0, 0.0, 1.0]
    Pw[2,1,0,:] = [1.0, 0.0, 0.0, 1.0]

    Pw[0,0,1,:] = [0.0, 0.5, 1.0, 1.0]
    Pw[1,0,1,:] = [0.5, 0.5, 1.0, 1.0]
    Pw[2,0,1,:] = [0.5, 0.0, 1.0, 1.0]
    Pw[0,1,1,:] = [0.0, 1.0, 1.0, 1.0]
    Pw[1,1,1,:] = [1.0, 1.0, 1.0, 1.0]
    Pw[2,1,1,:] = [1.0, 0.0, 1.0, 1.0]

    Pw[1,:,:,:] *= np.sqrt(2)/2

    Vx, Vy, Vz, Qw = vol.DegreeElevate(px,Ux,
                                       py,Uy,
                                       pz,Uz,
                                       Pw,
                                       0,0,0)
    assert np.allclose(Ux, Vx)
    assert np.allclose(Uy, Vy)
    assert np.allclose(Uz, Vz)
    assert np.allclose(Pw, Qw)

    tx, ty, tz = 1, 2, 2
    Vx, Vy, Vz, Qw = vol.DegreeElevate(px,Ux,
                                       py,Uy,
                                       pz,Uz,
                                       Pw,
                                       tx,ty,tz)
    qx = px + tx
    qy = py + ty
    qz = pz + tz

    x = np.linspace(Vx[0], Vx[-1], 15)
    y = np.linspace(Vy[0], Vy[-1], 5)
    z = np.linspace(Vz[0], Vz[-1], 4)
    Cw = vol.Evaluate(qx,Vx,qy,Vy,qz,Vz,Qw,x,y,z)

    Q = Qw[:,:,:,:3] / Qw[:,:,:,3,None]
    C = Cw[:,:,:,:3] / Cw[:,:,:,3,None]

    if not PLOT: return
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.title("Volume - Degree Elevation")

    x = Q[:,:,:,0].flatten()
    y = Q[:,:,:,1].flatten()
    z = Q[:,:,:,2].flatten()
    ax.scatter(x,y,z, c='r')

    x = C[:,:,:,0].flatten()
    y = C[:,:,:,1].flatten()
    z = C[:,:,:,2].flatten()
    ax.scatter(x,y,z, c='b')

    t = np.linspace(0,np.pi/2,100)
    a = np.cos(t)
    b = np.sin(t)
    ax.plot(1.0*a,1.0*b,0,'-k')
    ax.plot(1.0*a,1.0*b,1,'-k')
    ax.plot(0.5*a,0.5*b,0,'-k')
    ax.plot(0.5*a,0.5*b,1,'-k')

    plt.axis("equal")


if __name__ == '__main__':
    VERB=1
    try:
        from matplotlib import pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        PLOT=1
    except ImportError:
        PLOT=0
    if 1: test_vol_ki(VERB=VERB, PLOT=PLOT)
    if 1: test_vol_de(VERB=VERB, PLOT=PLOT)
    if PLOT: plt.show()
