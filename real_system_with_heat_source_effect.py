__author__ = 'Marzuk Kamal'


import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
# from matplotlib.tri import triangulation
# from matplotlib import cm
# from mpl_toolkits.mplot3d import axes3d
# import itertools
#
# from matplotlib.tri import triangulation
# from scipy.spatial import ConvexHull
# from scipy.optimize import curve_fit
#
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.cm as cm
# import itertools
# # import Generator as G
#
# from matplotlib.tri import triangulation
# from scipy.spatial import ConvexHull
#
# from tempfile import NamedTemporaryFile, TemporaryFile


import MKMultiMaterialSystem as mm
reload(mm)

####
# data reference http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/thermal.html

L = 30
a = int(20)
b = a
c = 40

voxel_L = 10*1e-6 * 100 # 10 um to cm

voxel_volume_cm3 = voxel_L**3.0  #  cm3

I = 70e-3   # amp
V = 2.5     # volt

C_joul_per_gram = 0.33 # J g-1°C -1

rho_GaAs = 5.317  # g/cm^3
print rho_GaAs

rho_c = C_joul_per_gram*rho_GaAs   # joul/cm^3/C

Q = I * V/voxel_volume_cm3  # Joul/s
print 'Q = %.4e  J/cm3/s' % (Q,)


heat_fraction = 1
heat_source = Q/(rho_GaAs * C_joul_per_gram) * heat_fraction


D_val = 3.1e-5 # m^2/s
A0 = 20

x_star = float(a)/1000.0
y_star = x_star

boundary_temp = 25


mfactor = 1
# plt.ioff()
# model1 = mm.MKModel(model_id= 1, model_size=(16,16), D = 0.5, initial_condition = 4)
# model2 = mm.MKModel(model_id= 2, model_size=(4,4), D = 0.5, initial_condition = 3)
# model3 = mm.MKModel(model_id= 3, model_size=(3,3), D = 0.5, initial_condition = 2)
# main_system = mm.MKMultiMaterialSystem(base_model = model1, delta_pos=(0.1,0.1), time_step_scaling = 1, padding=(0.2,0.2))

# D is in m^2/sec
# model_size is unitless will be multiplied by delta_pos

model1 = mm.MKModel(model_id= 1, model_size=(20,50,50), D = 5e-5, initial_condition = 35)
model2 = mm.MKModel(model_id= 2, model_size=(5, 20, 20), D = 5e-5, initial_condition = 30)
model3 = mm.MKModel(model_id= 3, model_size=(5, 10, 10), D = 5e-5, initial_condition = 25)

# model1.plotModelSurface(101,0)

main_system = mm.MKMultiMaterialSystem(base_model = model1, delta_pos=0.1, time_step_scaling = 1/8.0, padding=(1, 0.2, 0.2))
main_system.addModel(model2, placement_position=(15+25,10,20))
main_system.addModel(model3, placement_position=(15+30,10,25))
main_system.finalizeSystemSettings()


main_system.imposeDirichletBC(mat_id = 1, wall = 'zbottom', bc_value= 50)
main_system.imposeDirichletBC(mat_id = 1, wall = 'xbottom', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'xtop', bc_value= 0)
# main_system.imposeDirichletBC(mat_id = 1, wall = 'xtop', bc_value= 50)


# main_system.ui = main_system.ui * 0.0 + boundary_temp
u_init = np.array(sp.copy(main_system.ui), dtype=np.float)


mm.plot2DImage(2, main_system.main_model[:,20,:], xAxisLabel= 'X (um)', yAxisLabel= 'Z (um)', showContour=False )
mm.plot2DImage(21, np.nan_to_num(u_init[:,20,:]), xAxisLabel= 'X (um)', yAxisLabel= 'Z (um)', showContour=False)



# apply Dirichlet BC
main_system.imposeDirichletBC(mat_id = 1, wall = 'zbottom', bc_value= boundary_temp)
main_system.imposeDirichletBC(mat_id = 1, wall = 'ztop', bc_value= boundary_temp)
main_system.imposeDirichletBC(mat_id = 1, wall = 'xbottom', bc_value= boundary_temp)
main_system.imposeDirichletBC(mat_id = 1, wall = 'xtop', bc_value= boundary_temp)
main_system.imposeDirichletBC(mat_id = 1, wall = 'ybottom', bc_value= boundary_temp)
main_system.imposeDirichletBC(mat_id = 1, wall = 'ytop', bc_value= boundary_temp)

# main_system.imposeNeumannBC(mat_id = 1, wall = 'ztop', bc_value = 0)
# main_system.imposeNeumannBC(mat_id = 1, wall = 'zbottom', bc_value = 0)
# main_system.imposeNeumannBC(mat_id = 1, wall = 'ytop', bc_value = 0)
# main_system.imposeNeumannBC(mat_id = 1, wall = 'ybottom', bc_value = 0)
# main_system.imposeNeumannBC(mat_id = 1, wall = 'xtop', bc_value = 0)
# main_system.imposeNeumannBC(mat_id = 1, wall = 'xbottom', bc_value = 0)

main_system.imposeHeatSource(mat_id= 1, heat_source_id= 1, Q = Q, c = C_joul_per_gram, \
                             rho = rho_GaAs, positionOffset= [20,0,0], sourceType='zplane',\
                             heatSourceRange = [100,2,2])


mm.plot2DImage(2, main_system.main_model[:,10,:], xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)' )
mm.plot2DImage(21, np.nan_to_num(u_init[:,10,:]), xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)')
# mm.plot2DImage(25, np.nan_to_num(main_system.D[:,10,:]), xAxisLabel= 'X (mm)', yAxisLabel= 'Z (mm)')

main_system.showAllMaterialContactWalls()
main_system.showAllNeumannBCWalls()
main_system.showAllDirichletBCWalls()

main_system.showAllHeatSources(figure_id= 10)


main_system.showAllWalls(figure_id=1001)

u_init = np.array(sp.copy(main_system.ui), dtype=np.float)

################

data, mid_temp = main_system.runSimulation_v2(iterationTime=20000, initU=u_init, figureID=401,\
                                    ySection= 4, applyDiri=True, applyNeumann=True, \
                                    deltaTimeFactor= 1, shutDownExternalHeatingAtTime=1*1e443, shutDownDuration=1e100) # , getZSectionData =


plt.title('Thermal distribution')
plt.show()

# np.savez_compressed('thermal_sumation_2D_numerical_solution_for_analytical_comparison_win_ymax_neumann_with_mid_point_temp_data_Neumann=True_afterFixingNeumann.npz', data = data, mid_temp = mid_temp)
# dt = np.load('thermal_sumation_2D_numerical_solution_for_analytical_comparison.npz')
# data = dt['data']


zs, ys, xs = main_system.main_model.shape
xx, yy, zz = np.meshgrid( np.linspace(0, xs-1, xs), np.linspace(0,ys-1,ys), np.linspace(0,zs-1,zs))
xx, yy = np.meshgrid( np.linspace(0, xs-1, xs), np.linspace(0,ys-1,ys))

str = 'D  = %.1e, dimension %s \n dt = %.1e sec ds  = %.2f mm' % (main_system.D_dict[1], main_system.main_model.shape, \
                                                           main_system.delta_time, main_system.delta_pos*1000.0,)

# mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[4,:,:]+main_system.u[3,:,:]+main_system.u[4,:,:]+main_system.u[2,:,:])/3.0), \
#                            figure_id= 701, rstride=2, cstride=2, title_str= str,\
#                            xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)', cbarLabel='T (C)')  #, xlim = (-10,60), zlim = (-30,70))2

zpos = int(main_system.main_model.shape[0]/2)

print 'zpos == ', zpos

# mm.plot2DSurfaceAndContour(zz,yy, np.nan_to_num((main_system.u[:,10,:])/1.0)
mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[10,:,:])/1.0), \
                           figure_id= 701, rstride=4, cstride=4, title_str= str,\
                           xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)', cbarLabel='T (C)', zlim = [-5,10])  #, xlim = (-10,60), zlim = (-30,70))2


# mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[:,:,:])/1.0), \
#                            figure_id= 702, rstride=8, cstride=8, title_str= str,\
#                            xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)', cbarLabel='T (C)')  #, xlim = (-10,60), zlim = (-30,70))2
# plt.suptitle(str)


# np.savez_compressed('I=70mA_V=2p5V_data.npz', data = data)
data_I50mAV2p5V =  np.load('I=50mA_V=2p5V_data.npz')
data_I70mAV2p5V =  np.load('I=70mA_V=2p5V_data.npz')


fig2 = plt.figure(17001)
plt.clf()
plt.cla()

plt.plot(data[0][:]*1000.0, data[1][1])
# plt.plot(data[0][:], mid_temp, color = 'green', marker = '.')
# plt.plot(data_t[:,0], data_t[:,1], marker='o', color='red')
# plt.plot(data_t[:,0], data_t[:,2], marker='o', color='cyan')
plt.xlabel('Time (ms)')
plt.ylabel('T (C)')
plt.title('T vs time ')
# plt.ylim(24.99,25.002)
# plt.xlim(0, 1)
plt.show()



data1 = data_I50mAV2p5V['data']
data2 = data_I70mAV2p5V['data']

fig2 = plt.figure(17002)
plt.clf()
plt.cla()

plt.plot(data1[0][:]*1000.0, data1[1][1])
plt.plot(data2[0][:]*1000.0, data2[1][1])
# plt.plot(data[0][:], mid_temp, color = 'green', marker = '.')
# plt.plot(data_t[:,0], data_t[:,1], marker='o', color='red')
# plt.plot(data_t[:,0], data_t[:,2], marker='o', color='cyan')
plt.xlabel('Time (ms)')
plt.ylabel('T (C)')
plt.title('T vs time ')
# plt.ylim(24.99,25.002)
# plt.xlim(0, 1)
plt.show()




fig2 = plt.figure(17003)
plt.clf()
plt.cla()

plt.plot(data1[0][:]*1000.0, data2[1][1] - data1[1][1])
# plt.plot(data[0][:], mid_temp, color = 'green', marker = '.')
# plt.plot(data_t[:,0], data_t[:,1], marker='o', color='red')
# plt.plot(data_t[:,0], data_t[:,2], marker='o', color='cyan')
plt.xlabel('Time (ms)')
plt.ylabel('T (C)')
plt.title('T vs time ')
# plt.ylim(24.99,25.002)
# plt.xlim(0, 1)
plt.show()


# data_analytic = np.load('thermal_sim3D_analytic_sm_data_300ms_10ms_step.npz')
# data_analytic = data_analytic['data_analytic1']


pnum, f, h = mm.TestExponentialFit(x = data[0][:]/1, y = np.nan_to_num(data[1][1][:]))
panalytic, f, h = mm.TestExponentialFit(x = data_analytic[:,0], y = data_analytic[:,1])


fig2 = plt.figure(5002)
plt.clf()
plt.cla()

hn, = plt.plot(data[0][:], np.nan_to_num(data[1][1][:]), label = 'Numerical', marker = '.', color = 'blue')
hnfit, = plt.plot(data[0][:], f(data[0][:], *pnum), color = 'cyan', label = 'Numerical data Exp fit')
ha, = plt.plot(data_analytic[:,0], data_analytic[:,1], marker='o', color='green', label = 'Analytical')
hafit, = plt.plot(data_analytic[:,0], f(data_analytic[:,0], *panalytic), color='red', label = 'Analycal data Exp fit')

title_str = 'Comparison of Numerical and Analytical solution of thermal diffusion'
plt.title(title_str)

rate_factor = pnum[1]/panalytic[1]
fit_text = """In fit equation y = a * exp(-b * x)
              Numerical data fit params %s
              Analytical data fit params %s
              Numerical temp falls %.3f time faster than analytic data""" % (pnum, panalytic, rate_factor)

plt.text(0.05,3, fit_text)
plt.xlabel('Time (sec)')
plt.ylabel('T (C)')
plt.legend(handles = [hn, hnfit, ha, hafit])
plt.ylim(0,5)
plt.xlim(0, 0.3)
plt.show()

pnum
panalytic


data_num = f(data[0][:],*pnum)
data_analytic = f(data[0][:], *panalytic)

data_diff = data_analytic - data_num


fig3 = plt.figure(5003)
plt.clf()
plt.cla()
hdiff, = plt.plot(data[0][:], data_diff, color = 'red', label = '$T_A - T_N$')
title_str = 'Difference between Analytical ($T_A$) and Numerical ($T_N$)\n temperature decay with time'
plt.title(title_str)
#
# rate_factor = pnum[1]/panalytic[1]
# fit_text = """In fit equation y = a * exp(-b * x)
#               Numerical data fit params %s
#               Analytical data fit params %s
#               Numerical temp falls %.3f time faster than analytic data""" % (pnum, panalytic, rate_factor)

# plt.text(0.05,3, fit_text)
plt.xlabel('Time (sec)')
plt.ylabel('$T_A - T_N$ (C)')
plt.legend(handles = [hdiff])
plt.ylim(0,1)
# plt.xlim(0, 1)
plt.show()
