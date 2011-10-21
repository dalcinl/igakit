import bootstrap
import numpy as np
from igakit.igalib import srf

def test_srf_ev(VERB=0, PLOT=0):
    if VERB: print(srf.Evaluate.__doc__)

    p = 2
    q = 1
    U = np.asarray([0,0,0, 1,1,1], dtype=float)
    V = np.asarray([0,0,     1,1], dtype=float)
    n = len(U)-1-(p+1)
    m = len(V)-1-(q+1)
    Pw = np.zeros((n+1,m+1,3))
    Pw[0,0,:] = [0.0, 0.5, 1.0]
    Pw[1,0,:] = [0.5, 0.5, 1.0]
    Pw[2,0,:] = [0.5, 0.0, 1.0]
    Pw[0,1,:] = [0.0, 1.0, 1.0]
    Pw[1,1,:] = [1.0, 1.0, 1.0]
    Pw[2,1,:] = [1.0, 0.0, 1.0]
    Pw[1,:,:] *= np.sqrt(2)/2

    u = np.linspace(U[0], U[-1], 31)
    v = np.linspace(V[0], V[-1], 10)

    Cw = srf.Evaluate(p,U,q,V,Pw,u,v)
    Dw = srf.Evaluate(p,U,q,V,Pw,U,V)

    P = Pw[:,:,:2] / Pw[:,:,2, None]
    C = Cw[:,:,:2] / Cw[:,:,2, None]
    D = Dw[:,:,:2] / Dw[:,:,2, None]

    if not PLOT: return
    plt.figure()
    plt.title("Surface - Evaluation")

    x1 = C[:,:,0]
    y1 = C[:,:,1]
    plt.plot(x1,y1,'.b')

    x2 = D[:,:,0]
    y2 = D[:,:,1]
    plt.plot(x2,y2,'og')

    x0 = P[:,:,0]
    y0 = P[:,:,1]
    plt.plot(x0,y0,'sr')

    t = np.linspace(0,np.pi/2,100)
    a = np.cos(t)
    b = np.sin(t)
    plt.plot(1.0*a,1.0*b,'-k')
    plt.plot(0.5*a,0.5*b,'-k')

    plt.axis("equal")

def test_srf_ki(VERB=0, PLOT=0):
    if VERB: print(srf.RefineKnotVector.__doc__)

    p = 2
    q = 1

    U = np.asarray([0,0,0, 1,1,1], dtype=float)
    V = np.asarray([0,0,     1,1], dtype=float)

    n = len(U)-1-(p+1)
    m = len(V)-1-(q+1)
    Pw = np.zeros((n+1,m+1,3))
    Pw[0,0,:] = [0.0, 0.5, 1.0]
    Pw[1,0,:] = [0.5, 0.5, 1.0]
    Pw[2,0,:] = [0.5, 0.0, 1.0]
    Pw[0,1,:] = [0.0, 1.0, 1.0]
    Pw[1,1,:] = [1.0, 1.0, 1.0]
    Pw[2,1,:] = [1.0, 0.0, 1.0]
    Pw[1,:,:] *= np.sqrt(2)/2

    X = np.asarray([])
    Y = np.asarray([])
    Ubar, Vbar, Qw = srf.RefineKnotVector(p,U,q,V,Pw,X,Y)
    assert np.allclose(U,  Ubar)
    assert np.allclose(V,  Vbar)
    assert np.allclose(Pw, Qw)

    X = np.asarray([.25, 0.5])#;X = np.asarray([])
    Y = np.asarray([0.5, .75])#;Y = np.asarray([])
    Ubar, Vbar, Qw = srf.RefineKnotVector(p,U,q,V,Pw,X,Y)

    u = np.linspace(Ubar[0], Ubar[-1], 31)
    v = np.linspace(Vbar[0], Vbar[-1], 10)
    Cw = srf.Evaluate(p,Ubar,q,Vbar,Qw,u,v)

    Q = Qw[:,:,:2] / Qw[:,:,2, None]
    C = Cw[:,:,:2] / Cw[:,:,2, None]

    if not PLOT: return
    plt.figure()
    plt.title("Surface - Knot Insertion")

    x = Q[:,:,0]
    y = Q[:,:,1]
    plt.plot(x,y,'sr')

    x = C[:,:,0]
    y = C[:,:,1]
    plt.plot(x,y,'.b')

    t = np.linspace(0,np.pi/2,100)
    a = np.cos(t)
    b = np.sin(t)
    plt.plot(1.0*a,1.0*b,'-k')
    plt.plot(0.5*a,0.5*b,'-k')

    plt.axis("equal")

def test_srf_de(VERB=0, PLOT=0):
    if VERB: print(srf.DegreeElevate.__doc__)

    p = 2
    q = 1
    U = np.asarray([0,0,0, 1,1,1], dtype=float)
    V = np.asarray([0,0,     1,1], dtype=float)
    n = len(U)-1-(p+1)
    m = len(V)-1-(q+1)
    Pw = np.zeros((n+1,m+1,3))
    Pw[0,0,:] = [0.0, 0.5, 1.0]
    Pw[1,0,:] = [0.5, 0.5, 1.0]
    Pw[2,0,:] = [0.5, 0.0, 1.0]
    Pw[0,1,:] = [0.0, 1.0, 1.0]
    Pw[1,1,:] = [1.0, 1.0, 1.0]
    Pw[2,1,:] = [1.0, 0.0, 1.0]
    Pw[1,:,:] *= np.sqrt(2)/2

    #X = np.asarray([.25])
    #Y = np.asarray([.75])
    #U, V, Pw = srf.RefineKnotVector(p,U,q,V,Pw,X,Y)

    Uh, Vh, Qw = srf.DegreeElevate(p,U,q,V,Pw,0,0)
    Uh, Vh, Qw = srf.DegreeElevate(p,U,q,V,Pw,1,0)
    Uh, Vh, Qw = srf.DegreeElevate(p,U,q,V,Pw,0,1)
    r, s = 1, 2
    Uh, Vh, Qw = srf.DegreeElevate(p,U,q,V,Pw,r,s)
    ph = p+r; qh = q+s;

    u = np.linspace(Uh[0], Uh[-1], 31)
    v = np.linspace(Vh[0], Vh[-1], 10)
    Cw = srf.Evaluate(ph,Uh,qh,Vh,Qw,u,v)

    Q = Qw[:,:,:2] / Qw[:,:,2, None]
    C = Cw[:,:,:2] / Cw[:,:,2, None]

    if not PLOT: return
    plt.figure()
    plt.title("Surface - Degree Elevation")

    x = Q[:,:,0]
    y = Q[:,:,1]
    plt.plot(x,y,'sr')

    x = C[:,:,0]
    y = C[:,:,1]
    plt.plot(x,y,'.b')

    t = np.linspace(0,np.pi/2,100)
    a = np.cos(t)
    b = np.sin(t)
    plt.plot(1.0*a,1.0*b,'-k')
    plt.plot(0.5*a,0.5*b,'-k')

    plt.axis("equal")

if __name__ == '__main__':
    VERB=1
    try:
        from matplotlib import pylab as plt
        PLOT=1
    except ImportError:
        PLOT=0
    if 1: test_srf_ev(VERB=VERB, PLOT=PLOT)
    if 1: test_srf_ki(VERB=VERB, PLOT=PLOT)
    if 1: test_srf_de(VERB=VERB, PLOT=PLOT)
    if PLOT: plt.show()
