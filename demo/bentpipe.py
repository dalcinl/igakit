import numpy as np
from igakit.nurbs import NURBS, transform
import igakit.plotting as plt

L  = 3.0
R0 = 1.0
R1 = 2.0
Rb = R0+R1

circle = np.zeros((9,4), dtype='d')
circle[:,:2] = [[ 1, 0],
                [ 1, 1],
                [ 0, 1],
                [-1, 1],
                [-1, 0],
                [-1,-1],
                [ 0,-1],
                [ 1,-1],
                [ 1, 0]]
circle[:,-1] = 1
circle[1::2,:] *= np.sqrt(2)/2

annulus = np.zeros((9,2,4), dtype='d')
t0 = transform().scale(R0)
t1 = transform().scale(R1)
annulus[:,0,:] = t0(circle)
annulus[:,1,:] = t1(circle)

bentpipe = np.zeros((9,2,5,4), dtype='d')
t0 = transform().move(L, 2)
t1 = transform().move(L/2, 2)
t2 = transform().move(0, 2)
t3 = transform()
t3.scale([np.sqrt(2), 1, np.sqrt(2)])
t3.rotate(-np.pi/4, 1).move(-Rb, 2)
t4 = transform()
t4.rotate(-np.pi/2, 1).move([Rb, 0, -Rb])
bentpipe[:,:,0,:] = t0(annulus)
bentpipe[:,:,1,:] = t1(annulus)
bentpipe[:,:,2,:] = t2(annulus)
bentpipe[:,:,3,:] = t3(annulus)
bentpipe[:,:,4,:] = t4(annulus)

Pw = bentpipe
U = [0,0,0, 1,1, 2,2, 3,3, 4,4,4]
V = [0,0, 1,1]
W = [0,0,0, 1,1, 2,2,2]
nrb = NURBS(Pw, [U,V,W])
assert nrb.degree == (2,1,2)

plt.use_backend('mayavi')

plt.figure()
plt.controlpoints(nrb)
plt.controlgrid(nrb)
plt.knotpoints(nrb)
plt.knotgrid(nrb)

plt.figure()
plt.knotpoints(nrb, color=(0,0,1))
plt.knotgrid(nrb,   color=(0,0,1))
plt.volume(nrb)

plt.show()
