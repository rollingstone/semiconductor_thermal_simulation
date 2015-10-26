from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

n_angles = 36
n_radii = 8

plt.ioff()

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
plt.ion()
ax = fig.gca(projection='3d')

ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)

plt.show()


v = np.linspace(0,9,10)
xx, yy = np.meshgrid(v,v)

xx = xx.flatten()
yy = yy.flatten()
zz = np.ones(yy.shape) * 10

plt.ion()
plt.figure(99001)
plt.clf()
ax = plt.gca(projection='3d')
ax.plot_trisurf(xx,yy,zz, cmap=cm.coolwarm, alpha = 0.5, linewidth=0.1)
plt.show()




