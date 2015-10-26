__author__ = 'kamal'


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri

X = np.load('./mydatars.npy')
# My data points are strictly positive. This doesn't work if I don't center about the origin.
X -= X.mean(axis=0)

rad = np.linalg.norm(X, axis=1)
zen = np.arccos(X[:,-1] / rad)
azi = np.arctan2(X[:,1], X[:,0])

tris = mtri.Triangulation(zen, azi)

fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X[:,0], X[:,1], X[:,2], triangles=tris.triangles, cmap=plt.cm.bone)
plt.show()


