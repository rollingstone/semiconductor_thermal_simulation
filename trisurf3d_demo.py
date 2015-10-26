from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np


xx  = np.linspace(0, 9,10)
yy = np.linspace(0, 9,10)
zz  = np.linspace(0, 9,10)

xx[np.all(xx > 0 and xx < 9)] = 0
yy[np.all(yy > 0 and yy < 9)] = 0
zz[np.all(zz > 0 and zz < 9)] = 0


yy[0] = 1
yy[9] = 1

zz[0] = 1
zz[9] = 1

x,y,z = np.meshgrid(xx,yy,zz)



x = x.reshape((1,-1))
y = y.reshape((1,-1))
z = z.reshape((1,-1))

ff = plt.figure(990099)
plt.clf()
plt.cla()
ax = plt.gca(projection='3d')
ax.plot_trisurf(x,y,z)
plt.show()



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
