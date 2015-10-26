__author__ = 'Marzuk Kamal'


import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import itertools
# import Generator as G

from matplotlib.tri import triangulation
from scipy.spatial import ConvexHull

from tempfile import NamedTemporaryFile, TemporaryFile

import MKMultiMaterialSystem as mm
reload(mm)


def u_total_v3d(x, y, z, m, n, p, A0 = 10,  a = 1, b = 1, c = 1,  D = 1 , t_star = 1, t =0):

    # sval = (1 - (1 - (np.array(m, dtype = np.int) & 1))) * (1-(-1)**p)

    # print 'x y z', x, y, z
    #
    x /= (a-1)
    y /= (b-1)
    z /= (c-1)
    # print 'x y z', x, y, z

    pi2 = np.pi*np.pi

    sval = (1 - (-1)**m) * (1-(-1)**p)



    amp_val = 16.0*float(A0) * sval / (m*(2*n-1)*p)/(pi2*np.pi)

    # print 'amp val == ', np.sum(np.isnan(amp_val))



    # print 'sum amp_val = ', np.sum(amp_val)
    # lmn = pi2 *( (m/a)**2 + ((n-0.5)/b)**2 + (p/c)**2 )
    lmn = ( m**2 + ((n-0.5)/(float(b)/a) )**2 + (p/(float(c)/a))**2 )/t_star

    # print 'lm nan == ', np.sum(np.isnan(lmn))
    #
    # print 't_star ', t_star
    #
    # t_val= -lmn
    #
    # print 't_val', t_val
    #
    # print 't_val isnan ', np.sum(np.isnan(t_val))

    val = amp_val * np.sin((m*np.pi)*x) * np.sin( ((n-0.5)*np.pi)*y) * np.sin((p*np.pi)*z) * np.exp(-t * lmn)

    # print 'val == ', val

    return np.sum(val)


def create_mnp(L=10):
    mm,nn,pp = np.meshgrid(np.linspace(1,L,L), np.linspace(1,L,L), np.linspace(1,L,L))

    m = np.array(np.reshape(mm, (1,-1) ), dtype = np.float)
    n = np.array(np.reshape(nn, (1,-1) ), dtype = np.float)
    p = np.array(np.reshape(pp, (1,-1) ), dtype = np.float)
    return  m,n,p


def update_temp_3d(u, m=[], n=[], p=[], a=10, b=10, c=10, data_unit =1, D = 1,  A0= 10, t = 0):

    xx0 = np.linspace(0, a-1, a)
    yy0 = np.linspace(0, b-1, b)
    zz0 = np.linspace(0, c-1, c)

    z2, y2, x2 = np.meshgrid(zz0, yy0, xx0)

    # x2 = xx
    # y2 = yy
    # z2 = zz

    # a = float(a) * data_unit
    # b = float(b) * data_unit
    # c = float(c) * data_unit

    t_star = ((a*data_unit/np.pi)**2)/D

    # print 't_star = %e sec  t/t_star == %e' % (t_star, t/t_star,)


    for z in zz0:
        for y in yy0:
            for x in xx0:
                u[z, y, x] = u_total_v3d(x, y, z, a=a,b=b,c=c, m=m, n=n, p=p, A0 = A0, t_star=t_star, t = t)
                # u[z, y, x] = u_total_v3d(x, y, z, m, n, p, a = a, b = b, c = c, A0 = A0,  t = t)

    print 'u update done at %.3e sec ' % (t,)

    xx, yy = None, None
    return x2, y2, z2, t_star




L = 40
a = int(20)
b = a
c = 4
distance_unit = 1e-3 # mm

D_val = 5e-5 # m
A0 = 10

x_star = float(a)/1000.0
y_star = x_star

u = np.zeros((c, b, a), dtype=np.float)

t = 30   # sec
# fig = plt.figure(1001)
# fig.show()

data_t = np.array([])

m,n,p = create_mnp(L)

x2,y2,z2, ts = update_temp_3d(u, a=a, b=b, c=c, m=m, n=n, p=p, data_unit=distance_unit, D = D_val,  A0 = A0, t = 0)

zsurface = 2
mm.plot2DSurfaceAndContour(x2[:,zsurface, :], y2[:,zsurface, :], u[zsurface,:,:], figure_id=1111, cstride=2, rstride=2, \
                           xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)', cbarLabel='T (C)')


ysurface = 10
mm.plot2DSurfaceAndContour(x2[:,:, ysurface].transpose(), y2[:,:, ysurface].transpose, u[:,ysurface,:], figure_id=1111)


data_t = np.array([])

mid_point = [c/2, b/2, a/2]

NN_points = np.array([[0,0,5], [0,0,-5], [0, 5, 0], [0,-5, 0]])

n_point = NN_points + mid_point

tp = tuple(map(tuple, n_point.transpose()))

data_t = np.array([])

for t in np.arange(0, 0.3, 0.01):
    x2,y2,z2, ts = update_temp_3d(u, a=a, b=b, c=c, m=m, n=n, p=p, data_unit=distance_unit, D = D_val,  A0 = A0, t = t)

    # zsurface = 2
    # mm.plot2DSurfaceAndContour(x2[:,zsurface, :], y2[:,zsurface, :], u[zsurface,:,:], figure_id=1111, cstride=2, rstride=2, xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)')
    # plt.show()
    #
    print 't = ', t

    v = u #[1:-1,1:,1:-1]
    mid_temp = u[mid_point[0],mid_point[1],mid_point[2]]

    side_temp = u[tp]

    st_xp = u[n_point[0, 0], n_point[0, 1],n_point[0, 2]]
    st_xm = u[n_point[1, 0], n_point[1, 1],n_point[1, 2]]
    st_yp = u[n_point[2, 0], n_point[2, 1],n_point[2, 2]]
    st_ym = u[n_point[3, 0], n_point[3, 1],n_point[3, 2]]

    print 'side_temp ', st_xp, st_xm, st_yp, st_ym
    data_t = np.append(data_t, [t , np.mean(v), mid_temp, st_xp, st_xm, st_yp, st_ym])
    # data_t2 = np.append(data_t2, [t , st_xp, st_xm, st_yp, st_ym])
    # data_t = np.append(data_t, side_temp)


data_t = np.reshape(data_t, (-1,3+4))
data_t2 = np.reshape(data_t2, (-1,5))

# np.savez_compressed('thermal_sim3D_analytic_sm_data_300ms_10ms_step_with_midpoint_data.npz', data_analytic1 = data_t)
dd  =np.load('thermal_sim3D_analytic_sm_data_300ms_10ms_step_with_midpoint_data.npz')
# dd = np.load('thermal_sim3D_analytic_sm_data.npz')

data_t = dd['data_analytic1']


plt.figure(5551)
plt.clf()
plt.cla()

h1, = plt.plot(data_t[:,0], data_t[:,1], color = 'red', marker = 'o', label = 'Mean Temprature')
h2, = plt.plot(data_t[:,0], data_t[:,2], color = 'blue', marker = 'o', label = 'MidPoint (2,10,10) Temprature')
plt.legend(handles = [h1, h2])
plt.xlabel('Time (sec)')
plt.ylabel('Temperature (C)')
plt.xlim(0,0.3)
title_str = 'Thermal diffusion of in 3D material with dimension (z,y,x) = %s mm^3\nTemperature measured at the ' % (u.shape,)
plt.title(title_str)
plt.xlim(0,0.3)
print 'Done!!'



plt.figure(5552)
plt.clf()
plt.cla()

h1, = plt.plot(data_t[:,0], data_t[:,1], color = 'red', marker = 'o', label = 'Mean Temprature')
h2, = plt.plot(data_t[:,0], data_t[:,2], color = 'blue', marker = 'o', label = 'MidPoint (2,10,10) Temprature')
hxp, = plt.plot(data_t[:,0], data_t[:,3], color = 'green', marker = 'o', label = 'MidPoint (2,10,15) Temprature')
hxm, = plt.plot(data_t[:,0], data_t[:,4], color = 'cyan', marker = 'o', label = 'MidPoint (2,10,5) Temprature')
hyp, = plt.plot(data_t[:,0], data_t[:,5], color = 'orange', marker = 'o', label = 'MidPoint (2,15,10) Temprature')
hym, = plt.plot(data_t[:,0], data_t[:,6], color = 'purple', marker = 'o', label = 'MidPoint (2,5,10) Temprature')
plt.legend(handles = [h1, h2,hxp,hxm, hyp, hym])
plt.xlabel('Time (sec)')
plt.ylabel('Temperature (C)')
plt.xlim(0,0.3)
title_str = 'Thermal diffusion of in 3D material with dimension (z,y,x) = %s mm^3\nTemperature measured at the ' % (u.shape,)
plt.title(title_str)
plt.xlim(0,0.3)
print 'Done!!'





x2,y2,z2, ts = update_temp_3d(u, a=a, b=b, c=c, m=m, n=n, p=p, data_unit=distance_unit, D = D_val,  A0 = A0, t = 0)
title_str = 'Thermal distribution at z = 2 at t = %.3f' % 0,
zsurface = 2
mm.plot2DSurfaceAndContour(x2[:,zsurface, :], y2[:,zsurface, :], u[zsurface,:,:], figure_id=881, cstride=2, rstride=2, xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)',zlim=(-5,10), title_str=title_str)
plt.show()


x2,y2,z2, ts = update_temp_3d(u, a=a, b=b, c=c, m=m, n=n, p=p, data_unit=distance_unit, D = D_val,  A0 = A0, t = 0.02)
title_str = 'Thermal distribution at z = 2 at t = %.3f sec' % 0.02,
zsurface = 2
mm.plot2DSurfaceAndContour(x2[:,zsurface, :], y2[:,zsurface, :], u[zsurface,:,:], figure_id=882, cstride=2, rstride=2, xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)',zlim=(-5,10), title_str=title_str)
plt.show()


x2,y2,z2, ts = update_temp_3d(u, a=a, b=b, c=c, m=m, n=n, p=p, data_unit=distance_unit, D = D_val,  A0 = A0, t = 0.05)
title_str = 'Thermal distribution at z = 2 at t = %.3f sec' % 0.05,
zsurface = 2
mm.plot2DSurfaceAndContour(x2[:,zsurface, :], y2[:,zsurface, :], u[zsurface,:,:], figure_id=883, cstride=2, rstride=2, xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)',zlim=(-5,10), title_str=title_str)
plt.show()


x2,y2,z2, ts = update_temp_3d(u, a=a, b=b, c=c, m=m, n=n, p=p, data_unit=distance_unit, D = D_val,  A0 = A0, t = 0.1)
title_str = 'Thermal distribution at z = 2 at t = %.3f sec' % 0.1,
zsurface = 2
mm.plot2DSurfaceAndContour(x2[:,zsurface, :], y2[:,zsurface, :], u[zsurface,:,:], figure_id=884, cstride=2, rstride=2, xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)',zlim=(-5,10), title_str=title_str)
plt.show()




