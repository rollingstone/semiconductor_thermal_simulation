__author__ = 'kamal'

import numpy as np

from matplotlib import pyplot as plt


nx, ny = (5, 5)

x = np.linspace(0, 1,  nx)
y = np.linspace(0, 1,  ny)


xv, yv = np.meshgrid(x,y)

xv
yv


vxy = x,y
xv, yv = np.meshgrid(vxy)

xv
yv


zz = xv * xv + yv

plt.contourf(zz)


a3 = np.linspace(0 ,8, 8)
a3


a3 = np.reshape(a3, (2,2,2))
a3


def ff(a,b):
    a = a*b


v1 = [1,2,3]
v2 = 2


ff(v1, v2)

v1
v2

##

v = np.zeros([6,6])
v

v11 = np.ones([3,3])
v11

xn = 3
yn = 3


v[3:,3:] = v11
v


a = [1,2,3]

v = np.ones([4,4,4])

v0 = np.zeros([2,2,2])

v
v[0:2,0:2,0:2] = v0

v





from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

n_angles = 36
n_radii = 8

# An array of radii
# Does not include radius r=0, this is to eliminate duplicate points
radii = np.linspace(0.125, 1.0, n_radii)

# An array of angles
angles = np.linspace(0, 2*np.pi, n_angles, endpoint=False)

# Repeat all angles for each radius
angles = np.repeat(angles[...,np.newaxis], n_radii, axis=1)

# Convert polar (radii, angles) coords to cartesian (x, y) coords
# (0, 0) is added here. There are no duplicate points in the (x, y) plane
x = np.append(0, (radii*np.cos(angles)).flatten())
y = np.append(0, (radii*np.sin(angles)).flatten())

# Pringle surface
z = np.sin(-x*y)

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)

plt.show()