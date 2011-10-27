import numpy as np
from igakit.nurbs import NURBS
from igakit.transform import transform

# -----

def line(p0=(0,0), p1=(1,0)):
    """
    p0         p1
    o-----------o
        +--> u
    """
    p0 = np.asarray(p0, dtype='d')
    p1 = np.asarray(p1, dtype='d')
    points = np.zeros((2,3), dtype='d')
    points[0,:p0.size] = p0
    points[1,:p1.size] = p1
    knots = [0,0,1,1]
    return NURBS(points, [knots])

def circle(radius=1, center=None, angle=None):
    Pi = np.pi
    #
    if angle is None:
        start, end = 0, 2*Pi
    elif isinstance(angle, (tuple, list)):
        start, end = angle
        if start is None: start = 0
        if end is None: end = 2*Pi
    else:
        start, end = 0, angle
    quadrants = (0.0, Pi/2, Pi, 3*Pi/2)
    sweep = end - start
    if sweep < 0: sweep = 2*Pi + sweep
    spans = np.searchsorted(quadrants, abs(sweep))
    #
    alpha = sweep/(2*spans)
    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)
    tan_a = np.tan(alpha)
    x = radius*cos_a
    y = radius*sin_a
    w = cos_a
    m = x + y*tan_a
    Ca = np.array([[  x, -y, 0, 1],
                   [w*m,  0, 0, w],
                   [  x,  y, 0, 1]],
                  dtype='d')
    #
    Cw = np.empty((2*spans+1,4), dtype='d')
    R = transform().rotate(alpha+start, 2)
    Cw[0:3,:] = R(Ca)
    if spans > 1:
        R = transform().rotate(2*alpha, 2)
        for i in range(1, spans):
            n = 2*i+1
            Cw[n:n+2,:] = R(Cw[n-2:n,:])
    if center is not None:
        T = transform().translate(center)
        Cw = T(Cw)
    #
    a, b = 0, 1
    U = np.empty(2*(spans+1)+2, dtype='d')
    U[0], U[-1] = a, b
    U[1:-1] = np.linspace(a,b,spans+1).repeat(2)
    #
    return NURBS(Cw, [U])

def bilinear(p00=(-0.5,-0.5),
             p01=(+0.5,-0.5),
             p10=(-0.5,+0.5),
             p11=(+0.5,+0.5)):
    """
    o------------------o
    |p10  v         p11|
    |     ^            |
    |     |            |
    |p00  +----> u  p01|
    o------------------o
    """
    p00 = np.asarray(p00, dtype='d')
    p01 = np.asarray(p01, dtype='d')
    p10 = np.asarray(p10, dtype='d')
    p11 = np.asarray(p11, dtype='d')
    points = np.zeros((2,2,3), dtype='d')
    points[0,0,:p00.size] = p00
    points[0,1,:p01.size] = p01
    points[1,0,:p10.size] = p10
    points[1,1,:p11.size] = p11
    knots = [0,0,1,1]
    return NURBS(points, [knots]*2)

def trilinear(points=None):
    """
       +---------+
      /         /|
     /         / |
    +---------+  |
    |         |  |
    |  +      |  +
    |         | /
    |         |/
    +---------+
    """
    if points is None:
        points = np.mgrid[-.5:+.5:2j,-.5:+.5:2j,-.5:+.5:2j]
        points = np.rollaxis(points, 0, 4)
    else:
        points = np.asarray(points, dtype='d')
        if points.shape == (3,2,2,2):
            points = np.rollaxis(points, 0, 4)
    knots = [0,0,1,1]
    return NURBS(grid, [knots]*3)

# -----

def extrude(nrb, displ, axis=None):
    T = transform().translate(displ, axis)
    Cw = np.empty(nrb.shape+(2,4))
    Cw[...,0,:] = nrb.control
    Cw[...,1,:] = T(nrb.control)
    UVW = nrb.knots + ([0,0,1,1],)
    return NURBS(Cw, UVW)

def revolve(nrb, point, axis, angle=None):
    point = np.asarray(point, dtype='d')
    assert point.ndim in (0, 1)
    assert point.size <= 3
    axis = np.asarray(axis, dtype='d')
    assert axis.ndim in (0, 1)
    assert 1 <= axis.size <= 3
    if axis.ndim == 0:
        v = np.zeros(3, dtype='d')
        axis = (0,1,2)[axis]
        v[axis] = 1
    else:
        v = np.zeros(3, dtype='d')
        v[:axis.size] = axis
        norm_axis = np.linalg.norm(v)
        assert norm_axis > 0
        v /= norm_axis
    #
    n = [v[1], -v[0], 0]    # n = cross(v, z)
    gamma = np.arccos(v[2]) # cos_gamma = dot(v, z)
    T = transform().translate(-point).rotate(gamma, n)
    nrb = nrb.clone().transform(T)

    Cw = nrb.control
    X, Y, Z, W = (Cw[...,i] for i in range(4))
    rho = np.hypot(X, Y)
    alpha = np.arctan2(Y, X)
    alpha[alpha<0] += 2*np.pi
    sines = np.sin(alpha)
    cosines = np.cos(alpha)

    arc = circle(angle=angle)
    Aw = arc.control

    Qw = np.empty(nrb.shape + arc.shape + (4,))
    UVW = nrb.knots + arc.knots

    dot = np.dot
    for idx in np.ndindex(nrb.shape):
        z = Z[idx]
        w = W[idx]
        r = rho[idx]
        # for the sake of speed, inline
        # the transformation matrix
        # M = Rz(alpha)*Tz(z)*Sxy(rho)
        r_sin_a = r*sines[idx]
        r_cos_a = r*cosines[idx]
        M = np.identity(4, dtype='d')
        M[0,0] = r_cos_a; M[0,1] = -r_sin_a
        M[1,0] = r_sin_a; M[1,1] =  r_cos_a
        M[2,3] = z

        Qi = Qw[idx]
        Qi[...] = dot(Aw, M.T)
        Qi[...,3] *= w

    return NURBS(Qw, UVW).transform(T.invert())

# -----
