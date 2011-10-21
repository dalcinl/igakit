import bootstrap
import numpy as np
from igakit.igalib import crv

def test_crv_ki(VERB=0, PLOT=0):
    if VERB: print(crv.RefineKnotVector.__doc__)

    p = 2
    U = np.asarray([0,0,0, 1,1,1], dtype=float)
    n = len(U)-1-(p+1)
    Pw = np.zeros((n+1,3))
    Pw[0,:] = [0.0, 1.0, 1.0]
    Pw[1,:] = [1.0, 1.0, 1.0]
    Pw[2,:] = [1.0, 0.0, 1.0]
    Pw[1,:] *= np.sqrt(2)/2

    X = np.asarray([0.25,0.5,0.5,.8,.9])
    Ubar, Qw = crv.RefineKnotVector(p,U,Pw,X)

    u = np.linspace(U[0], U[-1], 31)
    Cw = crv.Evaluate(p,Ubar,Qw,u)
    Dw = crv.Evaluate(p,Ubar,Qw,Ubar)

    P = Qw[:,:2] / Qw[:,2, None]
    C = Cw[:,:2] / Cw[:,2, None]
    D = Dw[:,:2] / Dw[:,2, None]

    if not PLOT: return
    plt.figure()
    plt.title("Curve - Knot Insertion")

    x1 = C[:,0]
    y1 = C[:,1]
    plt.plot(x1,y1,'.b')

    x2 = D[:,0]
    y2 = D[:,1]
    plt.plot(x2,y2,'og')

    x0 = P[:,0]
    y0 = P[:,1]
    plt.plot(x0,y0,'sr')

    t = np.linspace(0,np.pi/2,100)
    a = np.cos(t)
    b = np.sin(t)
    plt.plot(a,b,'-k')

    plt.axis("equal")

def test_crv_kr(VERB=0, PLOT=0):
    if VERB: print(crv.RemoveKnot.__doc__)

    p = 2
    U = np.asarray([0,0,0, 1,1,1], dtype=float)
    n = len(U)-1-(p+1)
    Pw = np.zeros((n+1,3))
    Pw[0,:] = [0.0, 1.0, 1.0]
    Pw[1,:] = [1.0, 1.0, 1.0]
    Pw[2,:] = [1.0, 0.0, 1.0]
    Pw[1,:] *= np.sqrt(2)/2
    U0, P0 = U, Pw
    u = np.asarray([0.25, 0.5,0.5, 0.75])
    U, Pw = crv.RefineKnotVector(p,U,Pw,u)

    Ub, Qw = U,Pw
    for uu in u:
        t, Ub, Qw = crv.RemoveKnot(p,Ub,Qw,uu)
        if t > 0: Ub = Ub[:-t].copy()
        if t > 0: Qw = Qw[:-t,:].copy()
    assert np.allclose(Ub,U0)
    assert np.allclose(Qw,P0)

    u = np.linspace(U[0], U[-1], 31)
    Cw = crv.Evaluate(p,U,Pw,u)
    Dw = crv.Evaluate(p,Ub,Qw,u)

    P = Pw[:,:2] / Pw[:,2, None]
    Q = Qw[:,:2] / Qw[:,2, None]

    C = Cw[:,:2] / Cw[:,2, None]
    D = Dw[:,:2] / Dw[:,2, None]

    if not PLOT: return
    t = np.linspace(0,np.pi/2,100)
    a = np.cos(t)
    b = np.sin(t)

    plt.figure()

    plt.subplot(121)
    plt.title("Curve - Knot Removal")
    plt.plot(a,b,'-k')
    x = C[:,0]
    y = C[:,1]
    plt.plot(x,y,'ob')
    x = P[:,0]
    y = P[:,1]
    p = plt.plot(x,y,'sr')
    plt.axis("equal")
    plt.axis(xmax=1.1,ymax=1.1)

    plt.subplot(122)
    plt.title("Curve - Knot Removal")
    plt.plot(a,b,'-k')
    x = D[:,0]
    y = D[:,1]
    plt.plot(x,y,'ob')
    x = Q[:,0]
    y = Q[:,1]
    plt.plot(x,y,'sr')
    plt.axis("equal")
    plt.axis(xmax=1.1,ymax=1.1)

def test_crv_de(VERB=0, PLOT=0):
    if VERB: print(crv.DegreeElevate.__doc__)

    p = 2
    U = np.asarray([0,0,0, 1,1,1], dtype=float)
    n = len(U)-1-(p+1)
    Pw = np.zeros((n+1,3))
    Pw[0,:] = [0.0, 1.0, 1.0]
    Pw[1,:] = [1.0, 1.0, 1.0]
    Pw[2,:] = [1.0, 0.0, 1.0]
    Pw[1,:] *= np.sqrt(2)/2

    X = np.asarray([0.5])
    U, Pw = crv.RefineKnotVector(p,U,Pw,X)
    #t = 1
    #U, Pw = crv.DegreeElevate(p,U,Pw,t)
    #p = p + t

    t = 2
    Uh, Qw = crv.DegreeElevate(p,U,Pw,t)
    ph = p + t

    u = np.linspace(U[0], U[-1], 31)
    Cw = crv.Evaluate(ph,Uh,Qw,u)
    Dw = crv.Evaluate(ph,Uh,Qw,Uh)

    P = Qw[:,:2] / Qw[:,2, None]
    C = Cw[:,:2] / Cw[:,2, None]
    D = Dw[:,:2] / Dw[:,2, None]

    if not PLOT: return
    plt.figure()
    plt.title("Curve - Degree Elevation")

    x1 = C[:,0]
    y1 = C[:,1]
    plt.plot(x1,y1,'.b')

    x2 = D[:,0]
    y2 = D[:,1]
    plt.plot(x2,y2,'og')

    x0 = P[:,0]
    y0 = P[:,1]
    plt.plot(x0,y0,'sr')

    t = np.linspace(0,np.pi/2,100)
    a = np.cos(t)
    b = np.sin(t)
    plt.plot(a,b,'-k')

    plt.axis("equal")

if __name__ == '__main__':
    VERB=1
    try:
        from matplotlib import pylab as plt
        PLOT=1
    except ImportError:
        PLOT=0
    if 1: test_crv_ki(VERB=VERB,PLOT=PLOT)
    if 1: test_crv_kr(VERB=VERB,PLOT=PLOT)
    if 1: test_crv_de(VERB=VERB,PLOT=PLOT)
    if PLOT: plt.show()
