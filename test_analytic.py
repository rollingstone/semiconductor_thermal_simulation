__author__ = 'kamal'



import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import numpy as np
import time as tm


plt.ioff()


def plot2DSurfaceAndContour(X, Y, Z, figure_id = 1001, xlim = None, ylim = None, zlim = None, xlabel = 'X', ylabel = 'Y', zlabel = 'Z', cmap = cm.coolwarm , cticks = None, cbarLabel = ''):
    fig = plt.figure(figure_id)
    plt.ioff()
    plt.clf()
    plt.cla()
    ax = fig.gca(projection='3d')


    if xlim is None:
        xlim = (np.min(X), np.max(X))
    if ylim is None:
        ylim = (np.min(Y), np.max(Y))

    zm = np.max(Z)
    if zlim is None:
        zlim = (-zm, zm)

    zmin = np.min(Z)

    if cticks is None:
        ct = np.linspace(zmin, zm , 30)
    else:
        ct = cticks


    ph =ax.plot_surface(X, Y, Z, rstride=4, cstride=4, alpha=0.5, cmap= cmap)
    cset = ax.contourf(X, Y, Z, zdir='z', offset=zlim[0], alpha=0.7, cmap=cmap)
    cset = ax.contourf(X, Y, Z, zdir='x', offset=xlim[0], alpha=0.7, cmap=cmap)
    cset = ax.contourf(X, Y, Z, zdir='y', offset=ylim[1], alpha=0.7, cmap=cmap)

    ax.set_xlabel(xlabel)
    ax.set_xlim(*xlim)
    ax.set_ylabel(ylabel)
    ax.set_ylim(*xlim)
    ax.set_zlabel(zlabel)
    ax.set_zlim(*zlim)
    # plt.colorbar(ph, cmap = cmap, ticks = ct, bounds = ct)
    cbar = plt.colorbar(ph, cmap = cmap)
    cbar.set_label(cbarLabel, rotation = 270)

    # plt.colorbar(ph, ticks = ct)

    h = plt.draw()
    tm.sleep((0.05))

    return h



X, Y, Z = axes3d.get_test_data(0.03)

plot2DSurfaceAndContour(X,Y,Z, figure_id= 1001, xlim=(-40,40), ylim = (-40,40), zlim=(-100,100))

plt.figure(100)

plt.show()

def Thermal_Diffusion(X,Y, A = 1, pos0 = (0,0), t = 0):

    if t == 0:
        zz = X.copy()
        zz[:] = 0
        return A

    return A/t * np.exp(-0.001*((X - pos0[1])**2 + (Y - pos0[0])**2)/t)


# L = 100
#
# x2,y2 = np.meshgrid(np.linspace(0, L-1,L), np.linspace(0, L-1,L))
#
# z2 = Thermal_Diffusion(x2,y2, A= 100,  pos0=(L/2,L/2),t = 0.1)
#
# zmin  = np.min(z2)
# zmax = np.max(z2)
#
# ct = np.linspace(zmin, zmax, 30)
#
#
# for i in xrange(100):
#     z2 = Thermal_Diffusion(x2,y2, A= 100,  pos0=(L/2,L/2),t = 0.01+ float(i) * 1)
#     plot2DSurfaceAndContour(x2,y2,z2, zlim = (zmin, zmax), figure_id=20001, cticks= ct)



def B(n, f0, a, b):
    val  = f0 * (1 - (-1)**n) /(n * np.pi * np.sinh(n*np.pi * b/a))
    return val

def meu_val(n, a, x):
    return n * np.pi/a *x


def rectAnalyticSolution(u,  diriCond  = 0):
    ys, xs = u.shape

    yy, xx = np.meshgrid(np.linspace(0, ys-1, ys), np.linspace(0, xs-1, xs))

    # xx = np.reshape(xx,(1,-1))
    # yy = np.reshape(yy,(1,-1))

    xx = np.linspace(0, xs-1, xs)
    yy = np.linspace(0, ys-1, ys)

    # yy = np.reshape(yy,(1,-1))

    yBCPos = ys - 1
    yBCValue = diriCond

    n = np.linspace(1, 201, 200)

    a = xs - 1
    b = ys - 1

    for y in yy:
        for x in xx:
            val = B(n, diriCond, a, b) * np.sin(meu_val(n,a,x)) * np.sinh(meu_val(n, b, y))
            u[x ,y] = np.nansum(val)

    return np.meshgrid(yy, xx)

L = 100
u = np.zeros((L,L), dtype= np.float) + 10.0
yy, xx = rectAnalyticSolution(u, 100)

plot2DSurfaceAndContour(xx,yy, u, figure_id= 90001, cbarLabel='Temperature (C)')





xi = np.array([0., 0.5, 1.0])/2
yi = np.array([0., 0.5, 1.0])/2
zi = np.array([[0., 1.0, 2.0],
               [0., 1.0, 2.0],
               [-0.1, 1.0, 2.0]])/2

v = np.linspace(-.1, 2.0, 15, endpoint=True)
plt.contour(xi, yi, zi, v, linewidths=0.5, colors='k')
plt.contourf(xi, yi, zi, v, cmap=plt.cm.jet)
x = plt.colorbar(ticks=v)
print x
plt.show()