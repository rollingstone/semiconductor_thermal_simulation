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


def A(m,n, A0):

    lv = (np.array(m, dtype = np.int) & 1)
    sval = 1- (1 - (lv))

    val = 8*A0 * sval / (m*(2*n-1))
    # print 'A = ', val
    return val

def l_mn(m,n, a =1, b=1, D=1):
    val = -np.sqrt(D) * np.pi**2 *( (m/a)**2 + ((n - 0.5)/b)**2 )
    # print 'l_mn = ', val
    return val

def u_mn(x,y, m, n, A0 = 10, a = 1, b = 1, D = 1, t = 0):
    val = A(m,n, A0) * np.sin(m*np.pi*x/a) * np.sin((n-0.5)*np.pi*y/b) * np.exp(l_mn(m, n, a=a, b=b, D=D) * t)
    return val

def u_total(x, y, m, n, A0 = 10,  a = 1, b = 1, D = 1 , t =0):
    sval = 1 - (1 - (np.array(m, dtype = np.int) & 1))
    # sval = 1- (1 - (lv))
    amp_val = 8*A0 * sval / (m*(2*n-1))

    lmn = -np.sqrt(D) * np.pi**2 *( (m/a)**2 + ((n - 0.5)/b)**2 )

    val = amp_val * np.sin(m*np.pi*x/a) * np.sin((n-0.5)*np.pi*y/b) * np.exp(lmn * t)
    # val = np.sum(u_mn(x,y, m, n, A0=A0, a=a,b=b,D=D,t=t))
    return np.sum(val)


L = 20
a = 20
b = a
D_val = 5e-5 # m

A0 = 10

# for t in xrange(15):
t = 3.0   # sec


mm,nn = np.meshgrid(np.linspace(1,L+1,L), np.linspace(1,L+1,L))

m = np.array(np.resize(mm, (1, L*L))[0], dtype = np.float)
n = np.array(np.resize(nn, (1, L*L))[0], dtype = np.float)

xx, yy = np.meshgrid(np.linspace(0,a-1,a), np.linspace(0, b-1,b))

x2 = xx.copy()
y2 = yy.copy()

xx = np.array(np.reshape(xx, (1,-1)), dtype=np.int)[0]
yy = np.array(np.reshape(yy, (1,-1)), dtype=np.int)[0]

u = np.zeros((a,b), dtype=np.float)

# u_total(1, 1, m, n, a = a, b = b, D = 1, A0 = A0,  t = 5)


# # plt.ioff()
# plt.ion()

u_dict = {}

fig = plt.figure(1001)
fig.show()
plt.clf()
plt.cla()

count  = 0
for y in yy:
    for x in xx:
        u[y, x] = u_total(x, y, m, n, a = float(a)/1000.0, b = float(b)/1000.0, D = D_val, A0 = A0,  t = t)


    # print 'time t = ', txx


    # u_dict[t] = u


import MKMultiMaterialSystem as mm
reload(mm)


mm.plot2DSurfaceAndContour(x2, y2, u,fig = fig, cbarLabel='Temperature (C)', cstride=4, rstride=4)

    # plt.show()

plt.show()
print 'Done!!'




####


model1 = mm.MKModel(model_id= 1, model_size=(4,60,60), D = 5e-1, initial_condition = 35)
main_system = mm.MKMultiMaterialSystem(base_model = model1, delta_pos=(0.1, 0.1, 0.1), time_step_scaling = 1/8.0, padding=(1, 0.2, 0.2))
main_system.finalizeSystemSettings()


# apply Dirichlet BC
main_system.imposeDirichletBC(mat_id = 1, wall = 'zbottom', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'ztop', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'xbottom', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'xtop', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'ytop', bc_value= 50)


#  apply Neumann BC
neu_pos = main_system.imposeNeumannBCExcludinDirichlet(bc_value = 0)

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
    plt.title('Dirichlet BC')

plt.show()


################

u_init = sp.copy(main_system.u)

####

data = main_system.runSimulation(iterationTime=100000, initU=u_init, figureID=401, ySection= 5, applyDiri=True, applyNeumann=True)

plt.title('Thermal distribution')

plt.show()
zs, ys, xs = main_system.main_model.shape
xx, yy = np.meshgrid( np.linspace(0, xs-1, xs), np.linspace(0,ys-1,ys))

str = 'D  = %e, dimension %s \n dt = %e sec ds  = %e mm' % (main_system.D_dict[1], main_system.main_model.shape, main_system.delta_time, main_system.delta_pos[0],)
mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[4,:,:]+main_system.u[3,:,:]+main_system.u[4,:,:]+main_system.u[2,:,:])/3.0), figure_id= 701, rstride=2, cstride=2, title_str= str)  #, xlim = (-10,60), zlim = (-30,70))2


