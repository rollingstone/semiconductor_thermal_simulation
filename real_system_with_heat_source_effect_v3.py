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

import cProfile, pstats, StringIO
import re


import cPickle as pk

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


led_volume = 4e-6 # cm3

I = 16e-3   # amp
V = 2.5     # volt

C_joul_per_gram = 0.33 # J g-1°C -1

rho_GaAs = 5.317  # g/cm^3
print rho_GaAs

rho_c = C_joul_per_gram*rho_GaAs   # joul/cm^3/C

heat_fraction = 1


Q = I * V * heat_fraction#/voxel_volume_cm3   # Joul/s
print 'Q = %.4e  J/cm3/s' % (Q,)


# heat_source = Q/(rho_GaAs * C_joul_per_gram) * heat_fraction


D_val = 3.1e-5 # m^2/s
A0 = 20

x_star = float(a)/1000.0
y_star = x_star

boundary_temp = 25


mfactor = 1

# D are in cm2/s

# Aluminium
model1 = mm.MKModel(model_id= 1, model_size=(4, 50, 50), D = 0.8418, rho = 2.712, Cp =  0.902, initial_condition = 25) # IMS-PCB Alu
#  Plastic
# model2 = mm.MKModel(model_id= 2, model_size=(10, 40, 40), D = 3.4e-3, rho = 1.2, Cp =  1.5, initial_condition = 25) # Lead Frame... epoxy
model2 = mm.MKModel(model_id= 2, model_size=(10, 40, 40), D = 0.1, rho = 1.2, Cp =  1.5, initial_condition = 25) # Lead Frame... epoxy
# GaAs
model3 = mm.MKModel(model_id= 3, model_size=(10, 20, 20), D = 0.31, rho = 5.317, Cp =  0.33, initial_condition = 25) # GaAs P-N junction semiconductor
#  Gold
model4 = mm.MKModel(model_id= 4, model_size=(18, 20, 50), D = 1.27, rho = 19.320, Cp =  0.129, initial_condition = 25) # GaAs P-N junction semiconductor


# model1 = mm.MKModel(model_id= 1, model_size=(10, 40, 40), D = 1.2e-5, initial_condition = 25) # IMS-PCB Alu
# model2 = mm.MKModel(model_id= 2, model_size=(5* 1, 30, 30), D = 3.4e-7, initial_condition = 25) # Lead Frame... epoxy
# model3 = mm.MKModel(model_id= 3, model_size=(10, 20, 20), D = 3.1e-5, initial_condition = 25) # GaAs P-N junction semiconductor
#
# model4 = mm.MKModel(model_id= 4, model_size=(10, 20, 20), D = 3.1e-5, initial_condition = 25) # GaAs P-N junction semiconductor
# main_system = mm.MKMultiMaterialSystem(base_model = model1, delta_pos=voxel_L, delta_time= 1e-5, time_step_scaling = 1,\
#                                        padding=(2, 2, 2))

main_system = mm.MKMultiMaterialSystem(base_model = model1, padding=(4, 4, 4))
main_system.addModel(model2, onTopOfObject=model1)
main_system.addModel(model3, onTopOfObject=model2, placement_position=(0,0,0))
main_system.addModel(model4, onTopOfObject=model1, insertionPriority=-10)

# main_system.addModel(model4, onTopOfObject=model2, placement_position=(0,0,30))

main_system.finalizeSystemSettings()



main_system.imposeDirichletBC(mat_id = 1, wall = 'zbottom', bc_value= 25)
main_system.imposeNeumannBCExcludinDirichlet(bc_value=0)

# main_system.imposeHeatSource(heat_source_id=303, mat_id=1, Q=Q, c=C_joul_per_gram,rho=rho_GaAs,\
#                               positionOffset=(5,0,0))

main_system.imposeHeatSource(heat_source_id=303, mat_id=3, Q=Q, c=C_joul_per_gram,rho=rho_GaAs,\
                              positionOffset=(5,0,0))

# main_system.imposeHeatSource(heat_source_id=305, mat_id=4, Q=Q,c=C_joul_per_gram,rho=rho_GaAs,\
#                               positionOffset=(5,0,0))

u_init  = np.array(sp.copy(main_system.ui), dtype=np.float)

xz_frame = 25

mm.plot2DImage(2, main_system.main_model[:,xz_frame,:], xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)', showContour=False )
mm.plot2DImage(21, np.nan_to_num(u_init[:,xz_frame,:]), xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)' , showContour= False)

mm.plot2DImage(221, np.nan_to_num(main_system.D[:,xz_frame,:]), xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)' , showContour= True)
# mm.plot2DImage(25, main_system.D[:,10,:]*1e4, xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)', showContour= True)
# mm.plot2DImage(26, main_system.D[:,20,:]*1e4, xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)', showContour= True)
# mm.plot2DImage(27, main_system.D[:,xz_frame,:]*1e4, xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)', showContour= True)
# mm.plot2DImage(28, main_system.D[:,35,:]*1e4, xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)', showContour= True)


main_system.showAllWalls()
main_system.showAllNeumannBCWalls()
main_system.showAllDirichletBCWalls()
main_system.showAllHeatSources()

# u_init = np.array(sp.copy(main_system.ui), dtype=np.float)

################

# import profile
# pr = cProfile.Profile()
# pr.enable()
# full_data = []

data, mid_temp = main_system.runSimulation_v2(\
                                iterationTime=100000, initU=u_init, figureID=401,\
                                deltaTime= 1e-7,\
                                heatSourceVolume=20*20*10,\
                                voxelL = voxel_L,\
                                ySection= 25, applyDiri=True, \
                                applyNeumann=True, \
                                applyHeatDecay = False,\
                                shutDownExternalHeatingAtTime=1*1e33, shutDownDuration=1e100) # , getZSectionData =

# plt.title('Thermal distribution')
plt.show()

mm.plot2DImage(402, main_system.u[:,25,:], xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)', \
                    interpolationMethod='none')

plt.savefig('')


# np.savez_compressed('Neumann_0p1_data_v1.npz', data = data)
    # full_data.append(data)

# np.savez_compressed('data_10ms_simulation_100x200x200LED_500x400x400_Lead_Frame.npz', data = data)


data_combined = {}
time_data = np.array([])
for ix in xrange(full_data.__len__()):
    if ix == 0:
        time_data = full_data[ix][0]
    else:
        time_data = np.append(time_data, full_data[ix][0] + full_data[ix][0][-1] )




mm.saveObject('save_ms.pk1', main_system)
vv = mm.loadObjects('save_ms.pk1')
# pr.disable()
# s = StringIO.StringIO()
# sortby = 'cumulative'
# ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
# ps.print_stats()
# print s.getvalue()
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


fig2 = plt.figure(99001)
plt.clf()
plt.cla()
ax = plt.subplot(111)
h1, = plt.plot(data[0][:]*1e6, data[1][1][:,0], color = 'red', label = 'Mean temp', linewidth=3)
h2, = plt.plot(data[0][:]*1e6, data[1][1][:,1], color = 'green', label = 'Top 1/3rd mean temp')
h3, = plt.plot(data[0][:]*1e6, data[1][1][:,2], color = 'blue', label = 'Middle 1/3rd mean temp')
h4, = plt.plot(data[0][:]*1e6, data[1][1][:,3], color = 'cyan', label = 'Bottom 1/3rd mean temp')
plt.legend(handles = [h1,h2,h3,h4], loc = 2)
plt.xlabel('Time (us)')
plt.ylabel('T (C)')
locs,labels = plt.xticks()
plt.xticks(locs, map(lambda x: "%2.0e" % x, locs))
# ax.set_xscale('log')
plt.title('IMS Temp (C) vs Time (us) ')
plt.show()




fig2 = plt.figure(99002)
plt.clf()
plt.cla()
ax = plt.subplot(111)
h1, = plt.plot(data[0][:]*1e6, data[1][2][:,0], color = 'red', label = 'Mean temp', linewidth=3)
h2, = plt.plot(data[0][:]*1e6, data[1][2][:,1], color = 'green', label = 'Top 1/3rd mean temp')
h3, = plt.plot(data[0][:]*1e6, data[1][2][:,2], color = 'blue', label = 'Middle 1/3rd mean temp')
h4, = plt.plot(data[0][:]*1e6, data[1][2][:,3], color = 'cyan', label = 'Bottom 1/3rd mean temp')
plt.legend(handles = [h1,h2,h3,h4], loc=2)

plt.ylabel('T (C)')
plt.xlabel('Time (us)')
locs,labels = plt.xticks()
plt.xticks(locs, map(lambda x: "%2.0e" % x, locs))
# ax.set_xscale('log')
plt.title('Lead Frame Temp (C) vs Time (us)')
plt.show()


fig2 = plt.figure(99003)
plt.clf()
plt.cla()
ax = plt.subplot(111)
h1, = plt.plot(data[0][:]*1e6, data[1][3][:,0], color = 'red', label = 'Mean temp', linewidth=3)
h2, = plt.plot(data[0][:]*1e6, data[1][3][:,1], color = 'green', label = 'Top 1/3rd mean temp')
h3, = plt.plot(data[0][:]*1e6, data[1][3][:,2], color = 'blue', label = 'Middle 1/3rd mean temp')
h4, = plt.plot(data[0][:]*1e6, data[1][3][:,3], color = 'cyan', label = 'Bottom 1/3rd mean temp')
plt.xlabel('Time (us)')
plt.ylabel('T (C)')
locs,labels = plt.xticks()
plt.xticks(locs, map(lambda x: "%2.0e" % x, locs))
# ax.set_xscale('log')
plt.title('LED Temp (C) vs Time (us)')
plt.legend(handles = [h1,h2,h3,h4], loc=2)
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
plt.xlim(0, 1)
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
