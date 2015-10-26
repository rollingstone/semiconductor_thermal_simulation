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

# plt.ioff()

####

L = 30
a = int(20)
b = a
c = 4
distance_unit = 1e-3 # mm

D_val = 5e-5 # m
A0 = 10

x_star = float(a)/1000.0
y_star = x_star


mfactor = 1

model1 = mm.MKModel(model_id= 1, model_size=(c*mfactor,a*mfactor,a*mfactor), D = D_val*1.0, initial_condition = A0)
main_system = mm.MKMultiMaterialSystem(base_model = model1, delta_time = 1e-4, delta_pos = 0.1/mfactor, time_step_scaling = 1, padding=(1, 0.2, 0.2))
main_system.finalizeSystemSettings()

u_init = np.array(sp.copy(main_system.ui), dtype=np.float)


# apply Dirichlet BC
main_system.imposeDirichletBC(mat_id = 1, wall = 'zbottom', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'ztop', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'xbottom', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'xtop', bc_value= 0)
main_system.imposeDirichletBC(mat_id = 1, wall = 'ybottom', bc_value= 0)

#  apply Neumann BC
main_system.imposeNeumannBC(mat_id = 1, wall = 'ytop', bc_value = 0)

# main_system.imposeNeumannBC(mat_id = 1, wall = 'ztop', bc_value = 0)
# main_system.imposeNeumannBC(mat_id = 1, wall = 'zbottom', bc_value = 0)
# neu_pos = main_system.neumanBCCoordinates[0].coords
# neu_pos = main_system.imposeNeumannBCExcludinDirichlet(bc_value = 0)

# main_system.imposeHeatSource(mat_id= 1, heat_source_id= 1, Q = 100, c = 2, rho = 2, positionOffset= [1,0,0], sourceType='zplane')


mm.plot2DImage(2, main_system.main_model[:,15,:], xAxisLabel= 'X (mm)', yAxisLabel= 'Z (mm)' )
mm.plot2DImage(21, np.nan_to_num(u_init[:,15,:]), xAxisLabel= 'X (mm)', yAxisLabel= 'Z (mm)')
mm.plot2DImage(25, np.nan_to_num(main_system.D[:,15,:]), xAxisLabel= 'X (mm)', yAxisLabel= 'Y (mm)')

main_system.showAllMaterialContactWalls()
main_system.showAllNeumannBCWalls(ylim=[0,23])
main_system.showAllDirichletBCWalls()

main_system.showAllHeatSources()
################

data, mid_temp = main_system.runSimulation_v2(iterationTime=3000, initU=u_init, figureID=401,\
                                    ySection= 4, applyDiri=True, applyNeumann=True, \
                                    deltaTimeFactor= 1) # , getZSectionData =


plt.title('Thermal distribution')
plt.show()

# np.savez_compressed('thermal_sumation_2D_numerical_solution_for_analytical_comparison_win_ymax_neumann_with_mid_point_temp_data_Neumann=True_afterFixingNeumann.npz', data = data, mid_temp = mid_temp)
# dt = np.load('thermal_sumation_2D_numerical_solution_for_analytical_comparison.npz')
# data = dt['data']


zs, ys, xs = main_system.main_model.shape
xx, yy = np.meshgrid( np.linspace(0, xs-1, xs), np.linspace(0,ys-1,ys))

str = 'D  = %.1e, dimension %s \n dt = %.1e sec ds  = %.2f mm' % (main_system.D_dict[1], main_system.main_model.shape, \
                                                           main_system.delta_time, main_system.delta_pos*1000.0,)

# mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[4,:,:]+main_system.u[3,:,:]+main_system.u[4,:,:]+main_system.u[2,:,:])/3.0), \
#                            figure_id= 701, rstride=2, cstride=2, title_str= str,\
#                            xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)', cbarLabel='T (C)')  #, xlim = (-10,60), zlim = (-30,70))2

zpos = int(main_system.main_model.shape[0]/2)

print 'zpos == ', zpos

# mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[zpos,:,:])/1.0), \
mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[zpos,:,:])/1.0), \
                           figure_id= 701, rstride=4, cstride=4, title_str= str,\
                           xlabel='X (mm)', ylabel='Z (mm)', zlabel='T (C)', cbarLabel='T (C)', zlim = [-5,10])  #, xlim = (-10,60), zlim = (-30,70))2


# mm.plot2DSurfaceAndContour(xx,yy, np.nan_to_num((main_system.u[:,:,:])/1.0), \
#                            figure_id= 702, rstride=8, cstride=8, title_str= str,\
#                            xlabel='X (mm)', ylabel='Y (mm)', zlabel='T (C)', cbarLabel='T (C)')  #, xlim = (-10,60), zlim = (-30,70))2
# plt.suptitle(str)


fig2 = plt.figure(5001)
plt.clf()
plt.cla()

plt.plot(data[0][:], data[1][1][:])
# plt.plot(data[0][:], mid_temp, color = 'green', marker = '.')
plt.plot(data_t[:,0], data_t[:,1], marker='o', color='red')
# plt.plot(data_t[:,0], data_t[:,2], marker='o', color='cyan')
plt.xlabel('Time (sec)')
plt.ylabel('T (C)')

plt.xlim(0, 0.4)
plt.show()



data_analytic = np.load('thermal_sim3D_analytic_sm_data_300ms_10ms_step.npz')
data_analytic = data_analytic['data_analytic1']


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
