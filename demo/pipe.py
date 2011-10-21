import numpy as np
from igakit.nurbs import NURBS
import igakit.plotting as plt

circle = np.array(
    [[ 1, 0],
     [ 1, 1],
     [ 0, 1],
     [-1, 1],
     [-1, 0],
     [-1,-1],
     [ 0,-1],
     [ 1,-1],
     [ 1, 0]],
    dtype='d')

R0 = 0.5
R1 = 1.0
L  = 4.0

C = np.zeros((9,2,2,4), dtype='d')
C[:,0,0,:2] = circle*R0
C[:,1,0,:2] = circle*R1
C[:,:,0, 2] = 0
C[:,:,1,:2] = C[:,:,0,:2]
C[:,:,1, 2] = L

C[:,:,:,-1] = 1
C[1::2,:,:,:] *= np.sqrt(2)/2


U = [0,0,0, 1,1, 2,2, 3,3, 4,4,4]
V = [0,0, 1,1]
W = [0,0, 1,1]

pipe = NURBS(C, [U,V,W])

#pipe1 = pipe.clone().elevate(0,1,1)
#pipe2 = pipe.clone().refine([],[0.5],[0.25, 0.5, 0.75])

plt.use_backend('mayavi')

plt.figure()
plt.controlpoints(pipe)
plt.controlgrid(pipe)

plt.figure()
plt.knotpoints(pipe)
plt.knotgrid(pipe)

plt.figure()
plt.curve(pipe)
plt.figure()
plt.surface(pipe)

plt.figure()
plt.controlpoints(pipe)
plt.controlgrid(pipe)
plt.knotpoints(pipe)
plt.knotgrid(pipe)
plt.surface(pipe)

plt.show()
