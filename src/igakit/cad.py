import numpy as np
from igakit.nurbs import NURBS
from igakit.transform import transform

# -----

Pi = np.pi
deg2rad = np.deg2rad
rad2deg = np.rad2deg
radians = np.radians
degrees = np.degrees

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
    """
    Construct a NURBS circular arc or full circle

    Parameters
    ----------
    radius : float, optional
    center : array_like, optional
    angle : float or 2-tuple of floats, optional

    Examples
    --------

    >>> crv = circle()
    >>> crv.shape
    (9,)
    >>> P = crv.evaluate([0, 0.25, 0.5, 0.75, 1])
    >>> assert np.allclose(P[0], ( 1,  0, 0))
    >>> assert np.allclose(P[1], ( 0,  1, 0))
    >>> assert np.allclose(P[2], (-1,  0, 0))
    >>> assert np.allclose(P[3], ( 0, -1, 0))
    >>> assert np.allclose(P[4], ( 1,  0, 0))

    >>> crv = circle(angle=3*Pi/2)
    >>> crv.shape
    (7,)
    >>> P = crv.evaluate([0, 1/3., 2/3., 1])
    >>> assert np.allclose(P[0], ( 1,  0, 0))
    >>> assert np.allclose(P[1], ( 0,  1, 0))
    >>> assert np.allclose(P[2], (-1,  0, 0))
    >>> assert np.allclose(P[3], ( 0, -1, 0))

    >>> crv = circle(radius=2, center=(1,1), angle=(Pi/2,-Pi/2))
    >>> crv.shape
    (5,)
    >>> P = crv.evaluate([0, 0.5, 1])
    >>> assert np.allclose(P[0], (1,  3, 0))
    >>> assert np.allclose(P[1], (3,  1, 0))
    >>> assert np.allclose(P[2], (1, -1, 0))

    >>> crv = circle(radius=3, center=2, angle=Pi/2)
    >>> crv.shape
    (3,)
    >>> P = crv.evaluate([0, 1])
    >>> assert np.allclose(P[0], ( 5, 0, 0))
    >>> assert np.allclose(P[1], ( 2, 3, 0))

    """
    if angle is None:
        # Full circle, 4 knot spans, 9 control points
        spans = 4
        Cw = np.zeros((9,4), dtype='d')
        Cw[:,:2] = [[ 1, 0], [ 1, 1], [ 0, 1],
                    [-1, 1], [-1, 0], [-1,-1],
                    [ 0,-1], [ 1,-1], [ 1, 0]]
        wm = np.sqrt(2)/2
        Cw[:,3] = 1; Cw[1::2,:] *= wm
    else:
        Pi = np.pi # inline numpy.pi
        # Determine start and end angles
        if isinstance(angle, (tuple, list)):
            start, end = angle
            if start is None: start = 0
            if end is None: end = 2*Pi
        else:
            start, end = 0, angle
        # Compute sweep and number knot spans
        sweep = end - start
        quadrants = (0.0, Pi/2, Pi, 3*Pi/2)
        spans = np.searchsorted(quadrants, abs(sweep))
        # Construct a single-segment NURBS circular arc
        # centered at the origin and bisected by +X axis
        alpha = sweep/(2*spans)
        sin_a = np.sin(alpha)
        cos_a = np.cos(alpha)
        tan_a = np.tan(alpha)
        x = radius*cos_a
        y = radius*sin_a
        wm = cos_a
        xm = x + y*tan_a
        Ca = [[    x, -y, 0,  1],
              [wm*xm,  0, 0, wm],
              [    x,  y, 0,  1]]
        # Compute control points by successive rotation
        # of the controls points in the first segment
        Cw = np.empty((2*spans+1,4), dtype='d')
        R = transform().rotate(alpha+start, 2)
        Cw[0:3,:] = R(Ca)
        if spans > 1:
            R = transform().rotate(2*alpha, 2)
            for i in range(1, spans):
                n = 2*i+1
                Cw[n:n+2,:] = R(Cw[n-2:n,:])
    # Translate control points to center
    if center is not None:
        T = transform().translate(center)
        Cw = T(Cw)
    # Compute knot vector in the range [0,1]
    a, b = 0, 1
    U = np.empty(2*(spans+1)+2, dtype='d')
    U[0], U[-1] = a, b
    U[1:-1] = np.linspace(a,b,spans+1).repeat(2)
    # Return the new NURBS object
    return NURBS(Cw, [U])

def linear(points=None):
    """
    p[0]         p[1]
    o------------o
       +----> u
    """
    if points is None:
        points = np.zeros((2,2), dtype='d')
        points[0,0] = -0.5
        points[1,0] = +0.5
    else:
        points = np.asarray(points, dtype='d')
        assert points.shape[:-1] == (2,)
    knots = [0,0,1,1]
    return NURBS(points, [knots])

def bilinear(points=None):
    """
    p[0,1]       p[1,1]
    o------------o
    |  v         |
    |  ^         |
    |  |         |
    |  +----> u  |
    o------------o
    p[0,0]       p[1,0]
    """
    if points is None:
        s = slice(-0.5, +0.5, 2j)
        x, y = np.ogrid[s, s]
        points = np.zeros((2,2,2), dtype='d')
        points[...,0] = x
        points[...,1] = y
    else:
        points = np.array(points, dtype='d')
        assert points.shape[:-1] == (2,2)
    knots = [0,0,1,1]
    return NURBS(points, [knots]*2)

def trilinear(points=None):
    """
       p[0,1,1]     p[1,1,1]
       o------------o
      /|           /|
     / |          / |          w
    o------------o  |          ^  v
    | p[0,0,1]   | p[1,0,1]    | /
    |  |         |  |          |/
    |  o-------- | -o          +----> u
    | / p[0,1,0] | / p[1,1,0]
    |/           |/
    o------------o
    p[0,0,0]     p[1,0,0]
    """
    if points is None:
        s = slice(-0.5, +0.5, 2j)
        x, y, z = np.ogrid[s, s, s]
        points = np.zeros((2,2,2,3), dtype='d')
        points[...,0] = x
        points[...,1] = y
        points[...,2] = z
    else:
        points = np.array(points, dtype='d')
        assert points.shape[:-1] == (2,2,2)
    knots = [0,0,1,1]
    return NURBS(points, [knots]*3)

# -----

def compat(*nurbs):
    """

    Parameters
    ----------
    nurbs: sequence of NURBS


    Returns
    -------
    nurbs: list of NURBS

    """
    #
    def SameBounds(nurbs):
        degree = [nrb.degree for nrb in nurbs]
        knots  = [nrb.knots  for nrb in nurbs]
        shape  = (len(knots), len(knots[0]), 2)
        bounds = np.zeros(shape, dtype='d')
        for i, (degs, knts) in enumerate(zip(degree, knots)):
            for j, (p, U) in enumerate(zip(degs, knts)):
                bounds[i,j,:] = U[p], U[-p-1]
        umin = bounds[...,0].min(axis=0)
        umax = bounds[...,1].max(axis=0)
        limits = zip(umin, umax)
        for i, bnds in enumerate(bounds):
            knts = knots[i] = list(knots[i])
            for j, (ab, cd) in enumerate(zip(bnds, limits)):
                # u in [a, b] -> s in [d, d]
                a, b = ab; c, d = cd;
                if (a == c) and (b == d): continue
                # U <- (d-c)/(b-a) * (U-a) + c
                U = knts[j] = knts[j].copy()
                U -= a; U *= (d-c)/(b-a); U += c
        for i, knts in enumerate(knots):
            nurbs[i]._knots = tuple(knts)
    #
    def SameDegree(nurbs):
        # Ensure same degree by degree elevation
        degree = np.row_stack([nrb.degree for nrb in nurbs])
        degmax = degree.max(axis=0)
        times = degmax - degree
        for i, rst in enumerate(times):
            if np.any(rst):
                nurbs[i].elevate(*rst)
    #
    def MergeKnots(nurbs):
        degree = [nrb.degree for nrb in nurbs]
        knots  = [nrb.knots  for nrb in nurbs]
        kvalues = []
        for degs, knts in zip(zip(*degree), zip(*knots)):
            breaks = []; mults = [];
            for (p, U) in zip(degs, knts):
                # knot vector -> breaks & multiplicities
                s = np.empty(0, dtype='i')
                u, i = np.unique(U[p+1:-p-1], return_inverse=True)
                if i.size: s = np.bincount(i)
                breaks.append(u); mults.append(s)
            # Merge breaks and multiplicities
            u = np.unique(np.concatenate(breaks))
            s = np.zeros(u.size, dtype='i')
            for (ui, si) in zip(breaks, mults):
                mask = np.in1d(u, ui)
                s[mask] = np.maximum(s[mask], si)
            # Detemine knots to insert
            kvals = []
            for ui in breaks:
                vi = np.setdiff1d(u, ui) # breaks to insert
                ti = s[np.in1d(u, vi)]   # multiplicities
                kv = np.repeat(vi, ti)
                kvals.append(kv)
            kvalues.append(kvals)
        # Apply knot refinement
        for i, uvw in enumerate(zip(*kvalues)):
            if np.any(u.size for u in uvw):
                nurbs[i].refine(*uvw)
    #
    if len(nurbs) == 1:
        if not isinstance(nurbs[0], NURBS):
            nurbs = nurbs[0]
    nurbs = [nrb.clone() for nrb in nurbs]
    if len(nurbs) < 2: return nurbs
    assert (min(nrb.dim for nrb in nurbs) ==
            max(nrb.dim for nrb in nurbs))
    SameBounds(nurbs)
    SameDegree(nurbs)
    MergeKnots(nurbs)
    #
    return nurbs

# -----

def extrude(nrb, displ, axis=None):
    """
    Construct a NURBS surface/volume by
    extruding a NURBS curve/surface.

    Parameters
    ----------
    nrb : NURBS
    displ : array_like or float
    axis : array_like or int, optional

    Example
    -------

    >>> crv = circle()
    >>> srf = extrude(crv, displ=1, axis=2)

    >>> srf = bilinear()
    >>> vol = extrude(srf, displ=1, axis=2)

    """
    assert nrb.dim <= 2
    T = transform().translate(displ, axis)
    Cw = np.empty(nrb.shape+(2,4))
    Cw[...,0,:] = nrb.control
    Cw[...,1,:] = T(nrb.control)
    UVW = nrb.knots + ([0,0,1,1],)
    return NURBS(Cw, UVW)

def revolve(nrb, point, axis, angle=None):
    """
    Construct a NURBS surface/volume by
    revolving a NURBS curve/surface.

    Parameters
    ----------
    nrb : NURBS
    point : array_like
    axis : array_like or int
    angle : float, optional

    Example
    -------

    >>> crv = line(1,2)
    >>> srf = revolve(crv,  point=0, axis=2, angle=[Pi/2,2*Pi])
    >>> vol = revolve(srf, point=3, axis=1, angle=-Pi/2)

    """
    assert nrb.dim <= 2
    point = np.asarray(point, dtype='d')
    assert point.ndim in (0, 1)
    assert point.size <= 3
    axis = np.asarray(axis)
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
    # Transform the NURBS object to a new reference frame
    # (O,X,Y,Z) centered at point and z-oriented with axis
    n = [v[1], -v[0], 0]    # n = cross(v, z)
    gamma = np.arccos(v[2]) # cos_gamma = dot(v, z)
    T = transform().translate(-point).rotate(gamma, n)
    nrb = nrb.clone().transform(T)
    # Map cartesian coordinates (x,y,z) to cylindrical coordinates
    # (rho,theta,z) with theta in [0,2*pi] and precompute sines and
    # cosines of theta angles.
    Cw = nrb.control
    X, Y, Z, W = (Cw[...,i] for i in range(4))
    rho = np.hypot(X, Y)
    theta = np.arctan2(Y, X); theta[theta<0] += 2*np.pi
    sines, cosines = np.sin(theta), np.cos(theta)
    # Create a circular arc in the XY plane
    arc = circle(angle=angle)
    Aw = arc.control
    # Allocate control points and knots of the result
    Qw = np.empty(nrb.shape + arc.shape + (4,))
    UVW = nrb.knots + arc.knots
    # Loop over all control points of the NURBS object
    # to revolve taking advantage of NumPy nd-indexing
    dot = np.dot # inline numpy.dot
    zeros = np.zeros # inline numpy.zeros
    for idx in np.ndindex(nrb.shape):
        z = Z[idx]
        w = W[idx]
        r = rho[idx]
        r_sin_a = r*sines[idx]
        r_cos_a = r*cosines[idx]
        # for the sake of speed, inline
        # the transformation matrix
        # M = Rz(theta)*Tz(z)*Sxy(rho)
        M = zeros((4,4))
        M[0,0] = r_cos_a; M[0,1] = -r_sin_a
        M[1,0] = r_sin_a; M[1,1] =  r_cos_a
        M[2,3] = z
        M[3,3] = 1
        # Compute new 4D control points by transforming the
        # arc control point and tensor-product the weights
        Qi = Qw[idx]
        Qi[...] = dot(Aw, M.T)
        Qi[...,3] *= w
    # Create the new NURBS object and map
    # back to the original reference frame
    return NURBS(Qw, UVW).transform(T.invert())

def ruled(nrb1, nrb2):
    """
    Construct a ruled surface/volume
    between two NURBS curves/surfaces.

    Parameters
    ----------
    nrb1, nrb2 : NURBS

    """
    assert nrb1.dim == nrb2.dim
    assert nrb1.dim <= 2
    assert nrb2.dim <= 2
    nrb1, nrb2 = compat(nrb1, nrb2)
    Cw = np.zeros(nrb1.shape+(2,4),dtype='d')
    Cw[...,0,:] = nrb1.control
    Cw[...,1,:] = nrb2.control
    UVW = nrb1.knots + ([0,0,1,1],)
    return NURBS(Cw, UVW)

def sweep(section, trajectory):
    """
    Construct the translational sweep of a section
    curve/surface along a trajectory curve.

    S(u,v) = C(u) + T(v)

    V(u,v,w) = S(u,v) + T(w)

    Parameters
    ----------

    section : NURBS
        Section curve/surface
    trajectory : NURBS
        Trajectory curve

    """
    assert 1 <= section.dim <= 2
    assert trajectory.dim == 1
    Cs, ws = section.points, section.weights
    Ct, wt = trajectory.points, trajectory.weights
    C = Cs[...,np.newaxis,:] + Ct
    w = ws[...,np.newaxis] * wt
    UVW = section.knots + trajectory.knots
    return NURBS((C, w), UVW)

# -----
