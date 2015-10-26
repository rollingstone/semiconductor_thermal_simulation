__author__ = 'Marzuk Kamal'


import numpy as np
import scipy as sp
from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import itertools
# import Generator as G

from matplotlib.tri import triangulation
from scipy.spatial import ConvexHull

import MKMultiMaterialSystem as mm
reload(mm)

def u_total_0dim(x, y, m, n, A0 = 10,  a = 1, b = 1, D = 1 , t =0):
    sval = 1 - (1 - (np.array(m, dtype = np.int) & 1))
    amp_val = 8*A0 * sval / (m*(2*n-1))/(np.pi*np.pi)

    lmn = np.pi*np.pi *( (m*m) + (n*n) )
    val = amp_val * np.sin(m*np.pi*x) * np.sin((n-0.5)*np.pi*y) * np.exp(-lmn * t)
    # val = np.sum(u_mn(x,y, m, n, A0=A0, a=a,b=b,D=D,t=t))
    return np.sum(val)

def update_temp(u, a=10, a_unit =1, L = 10, D =1,  A0= 10, t = 0):

    mm,nn = np.meshgrid(np.linspace(1,L+1,L), np.linspace(1,L+1,L))

    m = np.array(np.resize(mm, (1, L*L))[0], dtype = np.float)
    n = np.array(np.resize(nn, (1, L*L))[0], dtype = np.float)

    xx, yy = np.meshgrid(np.linspace(0,a-1,a), np.linspace(0, a-1,a))

    x2 = xx.copy()
    y2 = yy.copy()

    xx = np.array(np.reshape(xx, (1,-1)), dtype=np.int)[0]
    yy = np.array(np.reshape(yy, (1,-1)), dtype=np.int)[0]

    x_star = np.array(sp.copy(xx), dtype = np.float)/(a-1)
    y_star = sp.copy(x_star)

    a = float(a) * a_unit

    t_star = a*a/D

    print 't_star = %e sec  t/t_star == %e' % (t_star, t/t_star,)

    for y in yy:
        for x in xx:
            u[y, x] = u_total_0dim(x_star[x], y_star[y], m, n, A0 = A0,  t = t/t_star)





    print 'u update done at %e sec ' % (t,)

    xx, yy = None, None
    return x2, y2, t_star

L = 40
a = int(20)
b = a
distance_unit = 1e-3

D_val = 5e-5 # m

A0 = 10


x_star = float(a)/1000.0
y_star = x_star

u = np.zeros((a,a), dtype=np.float)

t = 30   # sec
fig = plt.figure(1001)
fig.show()

data_t = np.array([])

for t in np.arange(0, 3, 0.5):

    x2, y2,ts = update_temp(u, a=a, a_unit=distance_unit, L= L, D = D_val,  A0 = A0, t = t)

    plt.clf()
    plt.cla()

    # import MKMultiMaterialSystem as mm
    # reload(mm)

    mm.plot2DSurfaceAndContour(x2, y2, u,fig = fig, cbarLabel='Temperature (C)', cstride=2, rstride=2, xlabel= 'X (mm)', ylabel= 'Y (mm)', zlabel='T (C)')
    plt.show()
    print 't = ', t

    data_t = np.append(data_t, [t , np.mean(u)])


x2, y2,ts = update_temp(u, a=a, a_unit=distance_unit, L= L, D = D_val,  A0 = A0, t = 0)
mm.plot2DSurfaceAndContour(x2, y2, u,figure_id=9901, cbarLabel='Temperature (C)', cstride=2, rstride=2, xlabel= 'X (mm)', ylabel= 'Y (mm)', zlabel='T (C)', zlim=[0,8])

x2, y2,ts = update_temp(u, a=a, a_unit=distance_unit, L= L, D = D_val,  A0 = A0, t = 1)
mm.plot2DSurfaceAndContour(x2, y2, u,figure_id = 9902, cbarLabel='Temperature (C)', cstride=2, rstride=2, xlabel= 'X (mm)', ylabel= 'Y (mm)', zlabel='T (C)', zlim=[0,8])

x2, y2,ts = update_temp(u, a=a, a_unit=distance_unit, L= L, D = D_val,  A0 = A0, t = 2)
mm.plot2DSurfaceAndContour(x2, y2, u,figure_id = 9903, cbarLabel='Temperature (C)', cstride=2, rstride=2, xlabel= 'X (mm)', ylabel= 'Y (mm)', zlabel='T (C)', zlim=[0,8])


data_t  = np.reshape(data_t, (-1,2))


plt.figure(555)
plt.clf()
plt.cla()

plt.plot(data_t[:,0], data_t[:,1], color = 'red')
plt.xlabel('Time (sec)')
plt.ylabel('Temperature (C)')
plt.xlim(0,2)
plt.title('thermal diffusion of in 2d square over time')

print 'Done!!'




####

import MKMultiMaterialSystem as mm
reload(mm)


mfactor = 2

model1 = mm.MKModel(model_id= 1, model_size=(3,a*mfactor,a*mfactor), D = D_val, initial_condition = A0)
main_system = mm.MKMultiMaterialSystem(base_model = model1, delta_pos=(0.1/mfactor, 0.1/mfactor, 0.1/mfactor), time_step_scaling = 1, padding=(1, 0.2, 0.2))
main_system.finalizeSystemSettings()


# apply Dirichlet BC
# main_system.imposeDirichletBC(mat_id = 1, wall = 'zbottom', bc_value= 0)
# main_system.imposeDirichletBC(mat_id = 1, wall = 'ztop', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'xbottom', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'xtop', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'ybottom', bc_value= 0)

#  apply Neumann BC
main_system.imposeNeumannBC(mat_id = 1, wall = 'ztop', bc_value = 0)
main_system.imposeNeumannBC(mat_id = 1, wall = 'zbottom', bc_value = 0)
neu_pos = main_system.imposeNeumannBC(mat_id = 1, wall = 'ytop', bc_value = 0)

neu_pos = main_system.neumanBCCoordinates[0].coords
# neu_pos = main_system.imposeNeumannBCExcludinDirichlet(bc_value = 0)

u_init  = sp.copy(main_system.ui)
mm.plot2DImage(2, main_system.main_model[:,15,:], xAxisLabel= 'X (mm)', yAxisLabel= 'Y (mm)' )
mm.plot2DImage(21, np.nan_to_num(u_init[:,15,:]), xAxisLabel= 'X (mm)', yAxisLabel= 'Y (mm)')
mm.plot2DImage(25, np.nan_to_num(main_system.D[:,15,:]), xAxisLabel= 'X (mm)', yAxisLabel= 'Y (mm)')

dpos = main_system.dirichletBCCoordinates[0].coords
# xpos = main_system.dirichletBCCoordinates[2].coords

# npos = main_system.neumanBCCoordinates[0].coords

# pos = main_system.getAllAirContactSurfacePoints(retainData = 0.4)
pos = main_system.getAllMaterialContactSurfacePoints(retainData=0.4)


f3 = plt.figure(1001)
plt.clf()
plt.cla()
ax = f3.add_subplot(111, projection = '3d')
# ax.invert_zaxis()
# ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = 'yellow')
ax.scatter(dpos[:,2], dpos[:,1], dpos[:,0], marker = '.', color = 'blue')
# ax.scatter(xpos[:,2], xpos[:,1], xpos[:,0], marker = '.', color = 'green')
# ax.scatter(xpos[:,2], xpos[:,1], xpos[:,0], marker = '.', color = 'yellow')
ax.scatter(neu_pos[:,2], neu_pos[:,1], neu_pos[:,0], marker = '.', color = 'red')
ax.set_xlabel('X')
ax.set_ylabel('Y')
# ax.invert_zaxis()
ax.set_zlabel('Z')
plt.show()



f3 = plt.figure(2001)
plt.clf()
plt.cla()
ax = f3.add_subplot(111, projection = '3d')
# ax.invert_zaxis()
ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = 'red')
# ax.scatter(dpos[:,2], dpos[:,1], dpos[:,0], marker = '.', color = 'blue')
# ax.scatter(xpos[:,2], xpos[:,1], xpos[:,0], marker = '.', color = 'green')
# ax.scatter(xpos[:,2], xpos[:,1], xpos[:,0], marker = '.', color = 'yellow')
ax.scatter(neu_pos[:,2], neu_pos[:,1], neu_pos[:,0], marker = '.', color = 'yellow')
ax.set_xlabel('X')
ax.set_ylabel('Y')
# ax.invert_zaxis()
ax.set_zlabel('Z')
plt.show()




len = main_system.neumanBCCoordinates.__len__()

colors = cm.rainbow(np.linspace(0,1,len))

cc = itertools.cycle(colors)

f3 = plt.figure(19001)
plt.clf()
plt.cla()
ax = f3.add_subplot(111, projection = '3d')
# ax.invert_zaxis()

for d in main_system.neumanBCCoordinates:
    #
    # if d.mat_id == 1:
    #     continue

    # if d.fluxDirection != mm.MKDirection.XMAX and d.fluxDirection != mm.MKDirection.XMIN:
    #     next(cc)
    #     continue

    # if d.fluxDirection != mm.MKDirection.ZMAX and d.fluxDirection != mm.MKDirection.ZMIN and d.fluxDirection != mm.MKDirection.YMAX and d.fluxDirection != mm.MKDirection.YMIN:
    #     next(cc)
    #     continue

    pos = d.coords
    print 'material id  %d    direction %d' % (d.mat_id, d.fluxDirection,)
    # print 'pos ', pos
    # if d.mat_id == 1:
    #     ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = 'red')
    # else:

    # cc = (0.1*float(d.mat_id), 0.1*float(d.mat_id), 0)
    # print 'color = ', cc
    ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = next(cc))
    #
    # raw_input('Next Wall')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('Neumann BC')

plt.show()



len = main_system.dirichletBCCoordinates.__len__()
colors = cm.rainbow(np.linspace(0,1,len))
cc = itertools.cycle(colors)

f3 = plt.figure(55001)
plt.clf()
plt.cla()
ax = f3.add_subplot(111, projection = '3d')
# ax.invert_zaxis()

for d in main_system.dirichletBCCoordinates:
    pos = d.coords
    print 'material id  %d    direction %d' % (d.mat_id, d.fluxDirection,)
    # print 'pos ', pos
    # if d.mat_id == 1:
    #     ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = 'red')
    # else:

    # cc = (0.1*float(d.mat_id), 0.1*float(d.mat_id), 0)
    # print 'color = ', cc
    ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = next(cc))
    #
    # raw_input('Next Wall')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('Dirichlet BC')

plt.show()


################

u_init = sp.copy(main_system.u)

####


data = main_system.runSimulation_v2(iterationTime=20000, initU=main_system.u, figureID=401, ySection= 4, applyDiri=True, applyNeumann=True, deltaTimeFactor= 1)

plt.title('Thermal distribution')

plt.show()
zs, ys, xs = main_system.main_model.shape
xx, yy = np.meshgrid( np.linspace(0, xs-1, xs), np.linspace(0,ys-1,ys))

str = 'D  = %.3e, dimension %s \n dt = %.3e sec ds  = %.3e mm' % (main_system.D_dict[1], main_system.main_model.shape, \
                                                           main_system.delta_time, main_system.delta_pos[0]*1000.0,)

# mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[4,:,:]+main_system.u[3,:,:]+main_system.u[4,:,:]+main_system.u[2,:,:])/3.0), \
#                            figure_id= 701, rstride=2, cstride=2, title_str= str,\
#                            xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)', cbarLabel='T (C)')  #, xlim = (-10,60), zlim = (-30,70))2

zpos = int(main_system.main_model.shape[0]/2)

print 'zpos == ', zpos

mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[zpos,:,:])/1.0), \
                           figure_id= 701, rstride=8, cstride=8, title_str= str,\
                           xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)', cbarLabel='T (C)', zlim = [-5,10])  #, xlim = (-10,60), zlim = (-30,70))2



# mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[:,:,:])/1.0), \
#                            figure_id= 702, rstride=8, cstride=8, title_str= str,\
#                            xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)', cbarLabel='T (C)')  #, xlim = (-10,60), zlim = (-30,70))2

# plt.suptitle(str)


fig2 = plt.figure(5001)
plt.clf()
plt.cla()

plt.plot(data[0][:]/1, np.nan_to_num(data[1][1][:]), marker='o')
plt.xlabel('Time (sec)')
plt.ylabel('T (C)')
# plt.xlim(0, 2)
plt.show()

