__author__ = 'Marzuk Kamal'


import matplotlib
# matplotlib.use('TkAgg')
from matplotlib import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import itertools
import numpy as np
import scipy as sp
# import Generator as G

from matplotlib.tri import triangulation
from scipy.spatial import ConvexHull

import MKMultiMaterialSystem as mm
reload(mm)


# plt.ioff()

# model1 = mm.MKModel(model_id= 1, model_size=(16,16), D = 0.5, initial_condition = 4)
# model2 = mm.MKModel(model_id= 2, model_size=(4,4), D = 0.5, initial_condition = 3)
# model3 = mm.MKModel(model_id= 3, model_size=(3,3), D = 0.5, initial_condition = 2)
# main_system = mm.MKMultiMaterialSystem(base_model = model1, delta_pos=(0.1,0.1), time_step_scaling = 1, padding=(0.2,0.2))

# D is in m^2/sec
# model_size is unitless will be multiplied by delta_pos

model1 = mm.MKModel(model_id= 1, model_size=(8,100,100), D = 0.5e-0, initial_condition = 10)
# model2 = mm.MKModel(model_id= 2, model_size=(5, 20, 20), D = 0.5e-0, initial_condition = 30)
# model3 = mm.MKModel(model_id= 3, model_size=(5, 10, 10), D = 0.5e-0, initial_condition = 25)

# model1.plotModelSurface(101,0)

main_system = mm.MKMultiMaterialSystem(base_model = model1, delta_pos=(0.1, 0.1, 0.1), time_step_scaling = 1, padding=(0.2, 0.2, 0.2))

# main_system.addModel(model2, placement_position=(15+25,10,20))
# main_system.addModel(model3, placement_position=(15+30,10,25))

main_system.finalizeSystemSettings()


main_system.imposeDirichletBC(mat_id = 1, wall = 'zbottom', bc_value= 50)
# main_system.imposeDirichletBC(mat_id = 1, wall = 'ztop', bc_value= 0.01)
# main_system.imposeDirichletBC(mat_id = 1, wall = 'ybottom', bc_value= 0.01)
# main_system.imposeDirichletBC(mat_id = 1, wall = 'ytop', bc_value= 50)
# main_system.imposeDirichletBC(mat_id = 1, wall = 'xbottom', bc_value= 0.01)
# main_system.imposeDirichletBC(mat_id = 1, wall = 'xtop', bc_value= 0.01)

#
# len = main_system.dirichletBCCoordinates.__len__()
#
# colors = cm.rainbow(np.linspace(0, 1,len))
# cc = itertools.cycle(colors)
# f3 = plt.figure(50001)
# plt.clf()
# plt.cla()
# ax = f3.add_subplot(111, projection = '3d')
# # ax.invert_zaxis()
#
# for d in main_system.dirichletBCCoordinates:
#     #
#     # if d.mat_id == 1:
#     #     continue
#
#     # if d.fluxDirection != mm.MKDirection.XMAX and d.fluxDirection != mm.MKDirection.XMIN:
#     #     next(cc)
#     #     continue
#
#     # if d.fluxDirection != mm.MKDirection.ZMAX and d.fluxDirection != mm.MKDirection.ZMIN and d.fluxDirection != mm.MKDirection.YMAX and d.fluxDirection != mm.MKDirection.YMIN:
#     #     next(cc)
#     #     continue
#
#     pos = d.coords
#     # print 'material id  %d    direction %d' % (d.mat_id, d.fluxDirection,)
#     # print 'pos ', pos
#     # if d.mat_id == 1:
#     #     ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = 'red')
#     # else:
#
#     # cc = (0.1*float(d.mat_id), 0.1*float(d.mat_id), 0)
#     # print 'color = ', cc
#     ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = next(cc))
#     #
#     # raw_input('Next Wall')
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#
# plt.show()

u_init  = sp.copy(main_system.ui)

################

data = main_system.runSimulation(iterationTime=1000, initU=u_init, figureID=777001, ySection= 15, applyNeumann=False)
plt.title('Thermal distribution')

zs,ys,xs = main_system.main_model.shape

yy, xx = np.meshgrid(np.linspace(0,ys-1,ys), np.linspace(0, xs-1, xs))

zsection = 4

# plt.ioff()
mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num(main_system.u[zsection,:,:]),figure_id=999001, cbarLabel='Temperature (C)' )

plt.show()


mm.plot2DImage(fig_id= 2001, image0 = main_system.u[4,:,:])




# main_system.imposeDirichletBC(mat_id = 1, wall = 'xtop', bc_value= 50)

# neu_pos = main_system.imposeNeumannBCExcludinDirichlet(bc_value = 0)

u_init  = sp.copy(main_system.ui)

mm.plot2DImage(2, main_system.main_model[:,15,:], xAxisLabel= 'X (mm)', yAxisLabel= 'Y (mm)' )

mm.plot2DImage(21, np.nan_to_num(u_init[:,15,:]), xAxisLabel= 'X (mm)', yAxisLabel= 'Y (mm)')

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
ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = 'red')
ax.scatter(dpos[:,2], dpos[:,1], dpos[:,0], marker = '.', color = 'blue')
# ax.scatter(xpos[:,2], xpos[:,1], xpos[:,0], marker = '.', color = 'green')
# ax.scatter(xpos[:,2], xpos[:,1], xpos[:,0], marker = '.', color = 'yellow')
ax.scatter(neu_pos[:,2], neu_pos[:,1], neu_pos[:,0], marker = '.', color = 'yellow')
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




# len = main_system.neumanBCCoordinates.__len__()
#
# colors = cm.rainbow(np.linspace(0,1,len))
#
# cc = itertools.cycle(colors)
#
# f3 = plt.figure(19001)
# plt.clf()
# plt.cla()
# ax = f3.add_subplot(111, projection = '3d')
# # ax.invert_zaxis()
#
# for d in main_system.dirichletBCCoordinates:
#     #
#     # if d.mat_id == 1:
#     #     continue
#
#     # if d.fluxDirection != mm.MKDirection.XMAX and d.fluxDirection != mm.MKDirection.XMIN:
#     #     next(cc)
#     #     continue
#
#     # if d.fluxDirection != mm.MKDirection.ZMAX and d.fluxDirection != mm.MKDirection.ZMIN and d.fluxDirection != mm.MKDirection.YMAX and d.fluxDirection != mm.MKDirection.YMIN:
#     #     next(cc)
#     #     continue
#
#     pos = d.coords
#     print 'material id  %d    direction %d' % (d.mat_id, d.fluxDirection,)
#     # print 'pos ', pos
#     # if d.mat_id == 1:
#     #     ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = 'red')
#     # else:
#
#     # cc = (0.1*float(d.mat_id), 0.1*float(d.mat_id), 0)
#     # print 'color = ', cc
#     ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = next(cc))
#     #
#     # raw_input('Next Wall')
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#
# plt.show()
#

################

data = main_system.runSimulation(iterationTime=1000, initU=u_init, figureID=777001, ySection= 15, applyNeumann=True)

plt.title('Thermal distribution')

########################


import numpy as np
def analylicSolutionPointSource(u,  pos, t, D = 1, init_value = 10, pos0 = (10,10,10), t0 = 0):

    # print 't = ', t
    # print 'D = ', D

    u[pos[:,0], pos[:,1], pos[:,2] ] = 1/(4.0 * np.pi * (t - t0) * D)**(3.0/2.0) *\
                                       np.exp( ((pos[:,2] - pos0[2])**2  +  (pos[:,1] - pos0[1])**2 +\
                                                (pos[:,0] - pos0[0])**2)/(4*t*D) \
                                             )


L = 20
init_value = 10e8

cube = np.zeros((L,L,L), dtype = np.float)
init_center = np.array([L/2,L/2,L/2], dtype=np.int)
cube[init_center[2], init_center[1], init_center[0]] = init_value

zz,yy,xx = np.meshgrid(np.linspace(0, L-1, L), np.linspace(0, L-1, L), np.linspace(0, L-1, L))

zz = np.reshape(zz, (1,-1))
yy = np.reshape(yy, (1,-1))
xx = np.reshape(xx, (1,-1))

pos = np.array(np.reshape(np.append(np.append(zz,yy), xx), (-1,3)), dtype = np.int)


delta_t = 0.01

tval = 100000

for i in xrange(1,tval):
    analylicSolutionPointSource(cube, pos = pos, t = float(i)*delta_t, D = 2, t0 = 0)



mm.plot2DImage(1000, cube[:,:,10])

z2, x2 = np.meshgrid(np.linspace(0,L-1,L), np.linspace(0,L-1,L))


plt.figure(8001)
plt.clf()
plt.cla()
ax = plt.subplot(111,projection='3d')
data = cube[:,10,:]*1e8

mx = np.max(data)
minx  =np.min(data)

ax.plot_surface(x2, z2, data, rstride=4, cstride=4, alpha=0.3)
cset = ax.contour(x2, z2, data, zdir='z', offset=-10, cmap=cm.coolwarm)
cset = ax.contour(x2, z2, data, zdir='x', offset=0, cmap=cm.coolwarm)
cset = ax.contour(x2, z2, data, zdir='y', offset=20, cmap=cm.coolwarm)

ax.set_xlabel('X')
ax.set_xlim(0, 20)
ax.set_ylabel('Y')
ax.set_ylim(0, 20)
ax.set_zlabel('Z')
# ax.set_zlim(0, 10)
ax.set_zlim(-10, 30)

plt.show()





def cube_marginals(cube, normalize=False):
    c_fcn = np.mean if normalize else np.sum
    xy = c_fcn(cube, axis=0)
    xz = c_fcn(cube, axis=1)
    yz = c_fcn(cube, axis=2)
    return(xy,xz,yz)

def plotcube(cube,x=None,y=None,z=None,normalize=False,plot_front=False):
    """Use contourf to plot cube marginals"""
    (Z,Y,X) = cube.shape
    (xy,xz,yz) = cube_marginals(cube,normalize=normalize)
    if x == None: x = np.arange(X)
    if y == None: y = np.arange(Y)
    if z == None: z = np.arange(Z)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # draw edge marginal surfaces
    offsets = (Z-1,0,X-1) if plot_front else (0, Y-1, 0)
    cset = ax.contourf(x[None,:].repeat(Y,axis=0), y[:,None].repeat(X,axis=1), xy, zdir='z', offset=offsets[0], cmap=plt.cm.coolwarm, alpha=0.75)
    cset = ax.contourf(x[None,:].repeat(Z,axis=0), xz, z[:,None].repeat(X,axis=1), zdir='y', offset=offsets[1], cmap=plt.cm.coolwarm, alpha=0.75)
    cset = ax.contourf(yz, y[None,:].repeat(Z,axis=0), z[:,None].repeat(Y,axis=1), zdir='x', offset=offsets[2], cmap=plt.cm.coolwarm, alpha=0.75)

    # draw wire cube to aid visualization
    ax.plot([0,X-1,X-1,0,0],[0,0,Y-1,Y-1,0],[0,0,0,0,0],'k-')
    ax.plot([0,X-1,X-1,0,0],[0,0,Y-1,Y-1,0],[Z-1,Z-1,Z-1,Z-1,Z-1],'k-')
    ax.plot([0,0],[0,0],[0,Z-1],'k-')
    ax.plot([X-1,X-1],[0,0],[0,Z-1],'k-')
    ax.plot([X-1,X-1],[Y-1,Y-1],[0,Z-1],'k-')
    ax.plot([0,0],[Y-1,Y-1],[0,Z-1],'k-')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
























######################
main_system.ui[:] = u_init[:]
time_steps = 100000
mm.tic()
for i in xrange(time_steps):
    main_system.evolveTempWithTime_v4(enableDiri=True, enableNeumann=False)

mm.plot2DImage(666001, np.nan_to_num(main_system.u[:,15,:]))

ax = plt.gca()
ax.set_xlabel('X (mm)')
ax.set_ylabel('Y (mm)')

print 'Total simulation time = %.4e seconds with time steps of %f ms' % (time_steps * main_system.delta_time, main_system.delta_time*1000.0,)
print 'Total time taken ', mm.toc(), 'seconds'





fig = plt.figure(1001)
plt.clf()
ax = plt.cla()

t = data[0]*1000000
m1 = data[1][1]
m2 = data[1][2]
m3 = data[1][3]

h1, = plt.plot(t, m1, color = 'red', label = 'Base')
h2, = plt.plot(t, m2, color = 'green', label = 'Middle')
h3, = plt.plot(t, m3, color = 'blue', label = 'Top')
ax = plt.gca()
# ax.xais.get_majot_formatter().set_scientific(True)

plt.legend(handles=[h1, h2, h3])
ax.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))

str = r"temperature (C) vs Time (us)\n of various materials in contact \n dpos = " + main_system.delta_pos.tostring() + "  dt = "

plt.title(str)
# plt.title(r'temperature (C) vs Time (us)\n of various materials in contact ')
plt.xlabel('t (us)')
plt.ylabel('T (C)')
plt.show()



####

f = plt.figure(2001)
plt.clf()
plt.cla()
im = plt.imshow(main_system.u[:,18,:])
# plt.colorbar(im)
plt.colorbar(im)
plt.show()



ff = plt.figure(3)
plt.clf()
plt.cla()
ax = ff.add_subplot(111, projecttion = '3d')
im = plt.imshow(main_system.u[:,15,:])
plt.imshow(main_system.ui[:,30,:])
ax.scatter3d(pos[:,2], pos[:,1], pos[:,0], marker = 'o' )
plt.title('')
plt.colorbar()






####################################

fig = plt.figure(2)
plt.clf()
plt.cla()
ax = fig.add_subplot(111, projection='3d')

plt.imshow(main_system.ui[:,30,:])
ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = 'o' )
plt.colorbar()


# c2 = model2.surfaceCoordinates
#
# fig9 = plt.figure(9)
# plt.clf()
# plt.cla()
# ax = fig9.add_subplot(111, projection='3d')
# Axes3D.plot_surface(c2[:,2], c2[:,1], c2[:,0])
#


plt.figure(2)
plt.clf()
plt.cla()
plt.imshow(u_init[:,30,:])
plt.colorbar()


plt.figure(3)
plt.clf()
plt.cla()
im = plt.imshow(main_system.u[:,30,:])
plt.colorbar(im)