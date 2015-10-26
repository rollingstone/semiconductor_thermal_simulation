
__author__ = 'Marzuk Kamal'


import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib.tri import triangulation
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import itertools

from matplotlib.tri import triangulation
from scipy.spatial import ConvexHull
from scipy.optimize import curve_fit

from joblib import Parallel, delayed
# import dill

import gc

import pickle



import time as tm
import multiprocessing.pool as pool
import numexpr as nx
from numba import  vectorize, jit, float32, float64, int32, int64

from multiprocessing import Pool


# import numpy as np
# import scipy as sp
# from matplotlib import pylab as plt
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.cm as cm
# import Generator as G


# from tempfile import NamedTemporaryFile, TemporaryFile
#
#
# from copy import copy


top = 1 << 3
bottom = -top

left = 1 << 2
right = -left

up = 1 << 1
down = -up

center = 1

class MKBoundaryCondition(object):
    def __init__(self, maximum_materials_count = 10):
        self.boundary_value_lookup = np.ones((maximum_materials_count, 1)) * -99
        self.boundary_type_lookup = np.ones((maximum_materials_count, 1)) * -99
        self.materials_in_contact = np.zeros((maximum_materials_count, 1))


    def addBoundaryCondition(self, mat1_id, mat2_id, boundaryConditionValue = 0, boundaryConditionType = 0):

        self.boundary_value_lookup[mat2_id] = boundaryConditionValue
        self.boundary_type_lookup[mat2_id] = boundaryConditionType
        self.materials_in_contact[mat1_id] = mat2_id

        # print

    def getBoundaryCondition(self, m1, m2):
        bc_value = self.boundary_value_lookup[self.materials_in_contact[m1]]
        bc_type = self.boundary_type_lookup[self.materials_in_contact[m1]]

        if bc_value == -99:
            print 'ERROR: boundary contition cannot be ', bc_value
            return -99

        return [bc_value, bc_type]


class MKSurfaceContact(object):
    def __init__(self):
        self.id = 0
        self.bcType = 0 # 0 for Dirichlet, 1 for Neumann
        self.bcValue = 5
        self.coords = np.array([], dtype=np.int)
        self.bcBirection = np.array([-1,0,0]) # + z by default


class MKDirection:
    ZMAX = 4
    ZMIN = -4

    YMAX = 3
    YMIN = -3

    XMAX = 2
    XMIN = -2

    Uniform = 0
    pass


class MKModeOjectData(object):
    def __init__(self, model = [], position = (0,0,0)):
        self.model = model
        self.position = np.array(position, dtype = np.int)

        x = int(position[2])
        y = int(position[1])
        z = int(position[0])

        zs,ys,xs = model.model.shape
        self.coords = np.array([], dtype=np.int)

        for zi in xrange(z, int(z+zs)):
            for yi in xrange(y, int(y+ys)):
                for xi in xrange(x, int(x+xs)):
                    self.coords = np.append(self.coords, [zi, yi, xi])

        self.coords =np.reshape(self.coords, (-1,3))

        print 'model_id %d self.coords shape == %s' % (model.model_id, self.coords.shape,)


class MKBoundaryConditionData(object):
    def __init__(self):
        self.mat_id = 0
        self.bcValue = 0
        self.bcType = 0 # 0 Dirichelet 1 Neumann
        self.coords = np.array([])
        self.fluxDirection = MKDirection.ZMAX # +z direction

class MKRuntimeData(object):
    def __init__(self):
        self.temperature = {}
        self.flux = {}


class MKHeatSourceData(object):
    def __init__(self, heat_id = None, mat_id = None, Q = 1, c = 1, rho = 1, fraction = 1):
        self.heat_id =  heat_id
        self.mat_id = mat_id
        self.Q = Q          #   J/cm3/s  Jin onevoxel
        self.c = c          #   J/g/C
        self.rho = rho      #   g/cm3
        self.fraction = fraction
        self.surface_area = 1
        self.total_voxel_volume = 1

        self.H = self.Q * self.fraction /(self.c * self.rho )
        self.H_per_pixel = self.H
        self.coords = np.array([])

        print 'heat_source H =Q/(rho*c) = %.4e C/s' % self.H,



    def changeHeatSourceParameters(self,Q = None, c = None, rho = None, fraction = None, \
                                   total_pixel = 1, voxel_volume = 1e-3):
        if Q is not None:
            self.Q =Q

        if c is not None:
            self.c = c

        if rho is not None:
            self.rho = rho

        if fraction is not None:
            self.fraction = fraction


        self.H = (self.Q * self.fraction)/(self.c*self.rho)
        self.H_per_pixel = self.H / self.surface_area # in pixels



class MKModel(object):
    def __init__(self, model_id = -1,  model_size = (10,10), D = 0.5, rho = 1, Cp = 1, initial_condition = 1, dirichlet_bc = 1, neumann_bc = 1):

        len = np.prod(model_size)

        if len == 0:
            print 'Error: model_size cannot have any zero values,  model_size == ', model_size
            return

        self.model_id = model_id
        self.model_ndim = model_size.__len__()
        self.model_size = np.array(model_size)
        self.model = np.zeros(self.model_size, dtype=np.byte)
        self.model[:] = self.model_id
        #  self.model = np.reshape(self.model, model_size)
        self.D = D
        self.rho = rho
        self.Cp = Cp

        self.global_offset = 0

        self.ui = np.zeros(self.model_size, dtype=np.float)
        self.ui[:] = initial_condition
        self.u = sp.copy(self.ui)
        self.dirichlet_boundary_condition_data = np.zeros(self.model_size, dtype = np.float) # with respect to air
        self.neumann_boundary_condition_data = np.zeros(self.model_size, dtype = np.float) # with respect to air
        self.dirichlet_bc_value =  dirichlet_bc
        self.neumann_bc_value = neumann_bc

        self.dirichlet_vector = []
        self.neumann_vector = []



        # self.surfaceCoordinates = self.getSurfaceCoordinates()
        # self.setBoundaryConditions()

    def isSurface(self, pos):

        if self.model_ndim == 3:
            s = self.model.shape
            xmax = s[2]
            ymax = s[1]
            zmax = s[0]

            x = pos[2]
            y = pos[1]
            z = pos[0]

            if x <= 0 or x+1 >= xmax:
                return True
            elif y <= 0 or y+1 >= ymax:
                return True
            elif z <= 0 or z+1 >= zmax:
                return  True

            if np.prod(self.model[z-1:z+2, y-1:y+2, x-1:x+2]) != 0:
                return  True

            return  False

        elif self.model_ndim == 2:
            s = self.model.shape
            ymax = s[0]
            xmax = s[1]

            y = pos[0]
            x = pos[1]


            if x < 0 or x >= xmax:
                return True
            elif y < 0 or y >= ymax:
                return True

            if np.prod(self.model[y-1:y+2, x-1:x+2]) != 0:
                return  True

            return  False
        else:
            raise Exception("ERROR: model dimension must be 2 or 3")

        return False

    def setBoundaryConditions(self):
        coords = self.getSurfaceCoordinates()

        for p in coords:
            if self.model_ndim == 3:
                self.dirichlet_boundary_condition_data[p[0], p[1], p[2]] = self.dirichlet_bc_value
                self.neumann_boundary_condition_data[p[0], p[1], p[2]] = self.neumann_bc_value

                self.dirichlet_vector.append(self.dirichlet_bc_value)
                self.neumann_vector.append(self.neumann_bc_value)

            elif self.model_ndim ==2:
                self.dirichlet_boundary_condition_data[p[0], p[1]] = self.dirichlet_bc_value
                self.neumann_boundary_condition_data[p[0], p[1]] = self.neumann_bc_value

                self.dirichlet_vector.append(self.dirichlet_bc_value)
                self.neumann_vector.append(self.neumann_bc_value)

            else:
                raise Exception("ERROR: worng dimension, 2 or 3 required!")


    def getSurfacePoints(self):

        zs,ys,xs = self.model.shape

        x = np.linspace(0, xs-1, xs)
        y = np.linspace(0, ys-1, ys)
        z = np.linspace(0, zs-1, zs)

        # xplane = [[z]]


    def getModel(self):
        return self.model

    def plotModelSurface(self, fig_id, filter_out = 0.95):
        surface_pos = self.getSurfaceCoordinates()

        r,c = surface_pos.shape

        print "r , c ", r,c
        count = 0
        data = np.array([])

        data_len = r * (1-filter_out)

        skip_val = int(np.ceil(r / data_len))

        print "skip_val = ", skip_val

        rval = np.array(np.round(np.random.rand(1, data_len) * (r-1)), dtype=np.int)

        print "rval ==", rval

        count = 0
        while count < r:
            data = np.append(data, surface_pos[count,:])
            count += skip_val

        #        data = data.reshape((-1,3))

        print "data = ", data.shape

        print data
        print "data.shape ", data.shape

        data = surface_pos

        fig = plt.figure(fig_id)
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(data[:,2], data[:,1], data[:,0])
        ax.set_aspect('equal')
        plt.xlabel('X data')
        plt.ylabel('Y data')
        # plt.zlabel('Z data')
        # forceAspect(ax, aspect=1)
        # plt.set_x_label
        fig.show()


class MKMultiMaterialSystem(object):

    NN_offset = np.array( [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])

    ztop = 1
    zbottom = -1

    ytop = 2
    ybottom = -2

    xtop = 3
    xbottom = -3



    # def __init__(self, base_model,  max_materials = 10, delta_time = 1e-3, delta_pos = 0.01, time_step_scaling = 1, padding = (0.2, 0.2, 0.2)):
    def __init__(self, base_model,  max_materials = 10, padding = (0.2, 0.2, 0.2)):

        gc.enable()

        np.set_printoptions(precision=8)

        # self.pool = Pool(processes = 6)

        self.base_model = base_model
        self.system_dimension = base_model.model_ndim

        self.main_model = [] #copy.deepcopy(base_model)
        self.model_object_base_coord = []
        self.boundaryCondition = MKBoundaryCondition()

        self.model_objects = {}

        self.neumanBCValue = {}
        self.dirichletBCValue = {}

        self.neumanBCCoordinates = None
        self.dirichletBCCoordinates = None

        self.heatSource = {}

        self.voxelL  = 10e-6*10 # 10 um in cm
        self.voxelArea  = self.voxelL**2 # 10 um in cm
        self.voxelVolume  = self.voxelL**3 # 10 um in cm

        self.simulation_result = []

        # self.model_objects.append(base_model)
        self.model_position = []

        self.all_bc_coords = []
        self.all_dirichlet_value = []
        self.all_neumann_value = []

        # self.dpos = delta_pos / 100.0   # convert from cm to m

        self.delta_pos = 1  # self.dpos
        self.delta_pos2 = 1 # self.dpos**2    #

        self.thermal_update = """self.ui[1:-1,1:-1,1:-1] + self.D[1:-1,1:-1,1:-1] * self.delta_time * \
                    ( (self.ui[2:, 1:-1, 1:-1] + self.ui[:-2, 1:-1, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1]) /self.delta_pos + \
                      (self.ui[1:-1, 2:, 1:-1] + self.ui[1:-1, :-2, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1]) /self.delta_pos + \
                      (self.ui[1:-1, 1:-1, 2:] + self.ui[1:-1, 1:-1, :-2] - 2 * self.ui[1:-1, 1:-1, 1:-1]) /self.delta_pos \
                    )"""

        new_size = np.round(self.base_model.model_size  * (1.0 + np.array(padding)))
        mx_val = np.max(new_size)
        new_size = np.array([mx_val, mx_val, mx_val])

        self.heatSourceVolume = None

        print 'self.base_model.model_size == ', self.base_model.model_size
        print 'new size', new_size
        self.main_model = np.zeros(new_size, dtype=np.int);
        self.ui = np.zeros(new_size, dtype=np.float)
        self.ui[:] = np.nan
        self.u = np.array(sp.copy(self.ui), dtype = np.float)

        print 'self.ui  type and self u type', self.ui.dtype

        self.zfilter = [] #np.logical_not(np.isnan(self.ui[pos[:,0]+1,pos[:,1],pos[:,2]] + self.ui[pos[:,0]-1,pos[:,1],pos[:,2]] - 2 * self.ui[pos[:,0],pos[:,1],pos[:,2]]))
        self.yfilter = [] #np.logical_not(np.isnan(self.ui[pos[:,0],pos[:,1]+1,pos[:,2]] + self.ui[pos[:,0],pos[:,1]-1,pos[:,2]] - 2 * self.ui[pos[:,0],pos[:,1],pos[:,2]]))
        self.xfilter = [] #np.logical_not(np.isnan(self.ui[pos[:,0],pos[:,1],pos[:,2]+1] + self.ui[pos[:,0],pos[:,1],pos[:,2]-1] - 2 * self.ui[pos[:,0],pos[:,1],pos[:,2]]))

        self.zfmatrix = []
        self.yfmatrix = []
        self.xfmatrix = []



        self.u_neumann = np.zeros(new_size, dtype=np.float)
        self.D = np.zeros(new_size, dtype=np.float) + np.nan
        self.RhoC = np.zeros(new_size, dtype=np.float)
        # self.D = np.ones(new_size, dtype=np.float)
        self.D_dict = {}

        self.mainModelObjectCoords = self.getMainModelObjectPositions() ## save all the object coords

        self.airContactSurface = {}
        self.materialContactSurface = {}
        self.materialTemperature = {}

        # self.modelDirichletBoundaryConditions = np.zeros(new_size, dtype= float)
        # self.modelNeumannBoundaryConditions = np.zeros(new_size, dtype= float)

        self.sxx = 1.0
        self.syy = 1.0
        self.szz = 1.0

        self.max_D = -9999999
        # self.min_deltaPosition = np.min(delta_pos)

        self.externalHeatSourceActive = True

        center_pos = new_size / 2;
        # self.base_model.model_size  * (1.0 + np.array(padding))
        offset_pos = center_pos - self.base_model.model_size/2;
        self.model_position.append(offset_pos);

        self.mat_index_dict = {}


        print "self.delta_pos == ", self.delta_pos


        self.delta_time = 1e-5
        # self.delta_time_scale = time_step_scaling

        self.inv_delta_pos = 1/self.delta_pos
        self.inv_delta_pos2 = 1/self.delta_pos2


        tt = 8*self.delta_pos2/self.D[1]

        # self.delta_time = delta_time

        self.addModel(self.base_model, placement_position=offset_pos)


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.pool.close()
        gc.collect()
        gc.disable()



    def addModel(self, new_model = [], onTopOfObject = None, offsetFromCenter = (0,0,0), placement_position = (0,0,0), place_next_to_object = -1, placement_mode = center, \
                 initial_condition = None, boundary_condition = 0, bc_type = 0, dirichlet_bc_value = 0, \
                 neumann_bc_value = 0, insertionPriority = 0):

        s = new_model.model.shape

        print 'model %d shape %s' % (new_model.model_id, s, )

        zoffset = 0
        yoffset = 0
        xoffset = 0

        if self.main_model.ndim == 3:

            if onTopOfObject is not None:
                pos = self.getMaterialCoords(onTopOfObject.model_id, recalculateAll= True)

                print 'pos ', pos.shape

                zmax = np.max(pos[:,0])

                print 'zmax ', zmax

                sbase = self.model_objects[onTopOfObject.model_id].model.model.shape

                zoffset = int(zmax+1)
                xoffset = int(np.min(pos[:,2]) + sbase[2]/2 - s[2]/2)
                yoffset = int(np.min(pos[:,1]) + 1*( sbase[1]/2 - s[1]/2))


                print 'zoffset ', zoffset
                print 'yoffset ', yoffset
                print 'xoffset ', xoffset

                z = placement_position[0] + zoffset + np.array([0, s[0]])
                y = placement_position[1] + yoffset + np.array([0, s[1]])
                x = placement_position[2] + xoffset + np.array([0, s[2]])

                print 'x  ==', x
                print 'y  ==', y
                print 'z  ==', z
            else:
                z = placement_position[0] + zoffset + np.array([0, s[0]])
                y = placement_position[1] + yoffset + np.array([0, s[1]])
                x = placement_position[2] + xoffset + np.array([0, s[2]])

            if insertionPriority < 0:
                idx = self.main_model[:] != 0
                mmatrix = sp.copy(self.main_model)
                self.main_model[z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.model_id
                self.main_model[idx] = mmatrix[idx]
                mmatrix = None
                idx = None
            else:
                self.main_model[z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.model_id


            # self.model_objects.append([new_model, np.array(placement_position, dtype=np.int)])
            # self.model_objects[new_model.model_id] = [new_model, np.array(placement_position, dtype = np.int)]
            self.model_objects[new_model.model_id] = MKModeOjectData(new_model, placement_position)

            # self.all_bc_coords.append(new_model.)

            # self.D[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.D

            print 'matid = %d D %.4e' % (new_model.model_id, new_model.D,)

            self.D[self.main_model[:] == new_model.model_id] = new_model.D
            self.RhoC[self.main_model[:] == new_model.model_id] = new_model.rho*new_model.Cp

            self.D_dict[new_model.model_id] = new_model.D

            if new_model.D >= self.max_D:
                self.max_D = new_model.D

            # self.modelDirichletBoundaryConditions[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.dirichlet_boundary_condition_data
            # self.modelNeumannBoundaryConditions[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.neumann_boundary_condition_data

            if initial_condition == None:
                print 'x == ', x
                print 'y == ', y
                print 'z == ', z

                print 'ui shape', self.ui.shape
                self.ui[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.ui
                # self.ui[ z[0]:z[1], y[0]:y[1], 120:150] = new_model.ui
            else:
                self.ui[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = initial_condition

            self.model_position.append(np.array(placement_position))
            self.boundaryCondition.addBoundaryCondition(0, new_model.model_id, boundary_condition, bc_type)

            # self.mainModelObjectCoords = self.getMainModelObjectPositions()

        elif self.main_model.ndim == 2:
            x = placement_position[1] + np.array([0, s[1]])
            y = placement_position[0] + np.array([0, s[0]])

            self.main_model[ y[0]:y[1], x[0]:x[1]] = new_model.model

            # self.model_objects.append([new_model, np.array(placement_position)])

            self.model_objects[new_model.model_id] = [new_model, np.array(placement_position)]

            self.neumanBCValue[new_model.model_id] = neumann_bc_value
            self.dirichletBCValue[new_model.model_id] = dirichlet_bc_value

            self.D[ y[0]:y[1], x[0]:x[1]] = new_model.D


            print "self.D === ", self.D

            #            self.modelDirichletBoundaryConditions[ y[0]:y[1], x[0]:x[1]] = new_model.

            if initial_condition == -9999:
                self.ui[ y[0]:y[1], x[0]:x[1]] = new_model.ui
            else:
                self.ui[ y[0]:y[1], x[0]:x[1]] = float(initial_condition)

            self.boundaryCondition.addBoundaryCondition(0, new_model.model_id, boundary_condition, bc_type)

            # self.D[new_model.model_id] = new_model.D

        else:
            print 'ERROR: Model must be 2 or 3 dimensional, given ', self.main_model.shape
            return



    def shrinkMainModel(self, pad = 8):

        idx = self.main_model[:] != 0

        pos = np.array(self.getMainModelObjectPositions(), dtype = np.int)

        xmin = np.min(pos[:,2])
        ymin = np.min(pos[:,1])
        zmin = np.min(pos[:,0])

        xsize = np.max(pos[:,2]) - xmin + 1
        ysize = np.max(pos[:,1]) - ymin + 1
        zsize = np.max(pos[:,0]) - zmin + 1

        offset = pad/2
        self.global_offset = offset

        new_xsize = xsize + pad
        new_ysize = ysize + pad
        new_zsize = zsize + pad

        new_model = np.zeros((new_zsize, new_ysize, new_xsize), dtype= np.int)

        new_model[offset:(offset+zsize), offset:(offset+ysize), offset:(offset+xsize)] = \
            self.main_model[zmin:(zmin+zsize), ymin:(ymin+ysize), xmin:(xmin+xsize)]

        self.main_model = np.array(sp.copy(new_model), dtype = np.int)

        new_model = np.array(sp.copy(new_model), dtype = np.float)

        new_model[offset:(offset+zsize), offset:(offset+ysize), offset:(offset+xsize)] = \
            self.ui[zmin:(zmin+zsize), ymin:(ymin+ysize), xmin:(xmin+xsize)]

        self.ui = np.array(sp.copy(new_model), dtype =np.float)
        self.ui[self.main_model[:] == 0] = np.nan

        self.u = sp.copy(self.ui)

        new_model[:] = 0.0

        new_model[offset:(offset+zsize), offset:(offset+ysize), offset:(offset+xsize)] = \
            self.D[zmin:(zmin+zsize), ymin:(ymin+ysize), xmin:(xmin+xsize)]

        self.D = np.array(sp.copy(new_model), dtype = np.float)

        new_model[:] = 0.0

        new_model[offset:(offset+zsize), offset:(offset+ysize), offset:(offset+xsize)] = \
            self.RhoC[zmin:(zmin+zsize), ymin:(ymin+ysize), xmin:(xmin+xsize)]

        self.RhoC = new_model.copy()
        # self.D = np.array(sp.copy(self.u), dtype = np.float)
        # self.D[:] = 0.0

        for mat_id, D in self.D_dict.iteritems():

            print 'Updating D   ',
            print 'mat_id == %d  D = %.3e ' % (mat_id, D)
            idx = self.main_model[:] == int(mat_id)

            print 'idx true count ', np.sum(idx.shape)
            # self.D[idx] = D

            old_pos = self.model_objects[mat_id].position
            self.model_objects[mat_id].position = [old_pos[0] - zmin + offset , old_pos[1] - ymin +  offset, old_pos[2] - xmin +  offset]

            self.mat_index_dict[mat_id] = idx


    def finalizeSystemSettings(self, pad = 4):

        if self.base_model.model_ndim == 3:
            self.shrinkMainModel(pad=pad)

            self.mainModelObjectCoords = np.array(self.getMainModelObjectPositions(), dtype = np.int)



            self.airContactSurface = {}
            self.materialContactSurface = {}

            self.calculateContactSurfaces()

            for mat_id in self.materialContactSurface.keys():
                coords = np.array(self.materialContactSurface[mat_id].coords, dtype=np.int)
                self.u_neumann[coords[:,0], coords[:,1], coords[:,2]] = self.materialContactSurface[mat_id].bcValue


            self.sxx = self.delta_time/self.delta_pos2
            self.syy = self.delta_time/self.delta_pos2
            self.szz = self.delta_time/self.delta_pos2

            print 'Data finalized'
            print 'Model shape ', self.main_model.shape


            for mat_id, d in self.D_dict.iteritems():
                print 'Material %d: D = %.3e m^2/s' % (mat_id, d,)
                self.model_objects[mat_id].position = self.getMaterialBasePosition(mat_id)


            self.d1 = sp.copy(self.main_model) *0.0
            self.d2 = sp.copy(self.main_model) *0.0
            self.d3 = sp.copy(self.main_model) *0.0



        # elif self.base_model.model_ndim == 2:
        #     last_err_setting = np.seterr(divide= 'ignore')
        #
        #     self.delta_time = np.prod(self.delta_pos2)/(2*2 * self.D * np.sum(self.delta_pos2)) / self.delta_time_scale
        #     self.delta_time[np.isinf(self.delta_time) | np.isnan(self.delta_time)] = 0
        #
        #     np.seterr(**last_err_setting)
        #
        # else:
        #     raise Exception("ERROR: Must be 2 or 3 dimensional model")

        # self.delta_time = np.min(self.delta_time[self.delta_time > 0])
        #
        # print "delta time is %.4e %s" % (self.delta_time,'  seconds')

        # self.airContactSurface, self.materialContactSurface =  self.calculateSurfaceContactsAndCoordinates()
        # self.calculateAirContactSurfaces(air_contacT_surface=self.airContactSurface)

        print 'self.airContactSurface  ', self.airContactSurface
        print 'self.materialContactSurface  ', self.materialContactSurface

        # self.airSurfaceCoords = self.getAirContactSurfaceCoordinates()


    def placeModelIntoMainTemplate(self, position = (0,0)):
        pass


    def isInLimit(self, pos):
        zs, ys, xs = self.main_model.shape

        xmin = min(abs(pos[2] - 1), 0)
        ymin = min(abs(pos[1] - 1), 0)
        zmin = min(abs(pos[0] - 1), 0)

        xmax = max(pos[2] + 1, xs-1)
        ymax = max(pos[1] + 1, ys-1)
        zmax = max(pos[0] + 1, zs-1)

        dx = xmax - xmin + 1
        dy = ymax - ymin + 1
        dz = zmax - zmin + 1


    def getMainModelObjectPositions(self):

        zs, ys, xs = self.main_model.shape

        t0 = tm.time()
        print 'model shape', self.main_model.shape

        # idx = self.main_model[:,:,:] != 0
        #
        # zz,yy,xx = np.meshgrid(np.linspace(0, zs-1,zs), np.linspace(0,ys-1,ys),\
        #                         np.linspace(0, xs-1, xs))
        #
        # pos = np.array([zz.flatten(), yy.flatten(), xx.flatten()], dtype=np.int).flatten().reshape((3,-1)).transpose()
        #
        # idx = self.main_model[pos[:,0], pos[:,1], pos[:,2]] != 0
        #
        # print '********************idx shape ', idx.shape
        #
        # pos0 = pos[idx,:]
        #
        # print 'pos0 shape ', pos0.shape
        #
        # pos = pos0
        #
        # zz = None
        # yy = None
        # xx = None
        #
        # return pos
        #
        # # print 'idx shape', idx.shape
        # xx = xx.flatten()
        # yy = yy.flatten()
        # zz = zz.flatten()
        #
        # pos = np.array([zz,yy,xx], dtype =  np.int).transpose()
        #
        # idx = self.main_model[pos[:,0], pos[:,1], pos[:,2]] != 0
        #
        # pos = pos[idx,:]

        #
        # pos = np.array(np.reshape( np.append(zp,np.append(yp,xp)), (3,-1)).transpose(), dtype = np.int)
        #
        # print 'pos shape', pos.shape
        # v1 = np.linspace(0,3,4)
        # v2 = np.linspace(0,5,6)
        # v3 = np.linspace(0,4,5)
        #
        # zz,yy,xx = np.meshgrid(v3,v2,v1)


        # xx = np.reshape(xx, (1,-1))
        # yy = np.reshape(yy, (1,-1))
        # zz = np.reshape(zz, (1,-1))
        #
        # ppp = np.append(zz,np.append(yy,xx))
        #
        # ppp = np.array(np.reshape(ppp, (3,-1)).transpose(), dtype = np.int)
        #
        # idx = np.logical_not(self.main_model[ppp[:,0], ppp[:,1], ppp[:,2]] == 0)
        #
        # pos = ppp[idx,:]

        xx = None
        yy = None
        zz = None




        # xx,yy,zz = np.meshgrid(np.linspace(0, xs-1,xs), np.linspace(0,ys-1,ys),\
        #                        np.linspace(0, zs-1, zs))
        #
        # L = zs*ys*xs
        #
        # zz = np.array(np.reshape(zz,(1,-1)), dtype=np.int)
        # yy = np.array(np.reshape(yy,(1,-1)), dtype=np.int)
        # xx = np.array(np.reshape(xx,(1,-1)), dtype=np.int)
        #
        # print 'xx yy zz ',
        # print xx.shape,
        # print yy.shape,
        # print zz.shape
        #
        # pos = np.reshape(np.append(zz,np.append(yy,xx)),(3,-1)).transpose()
        #
        # idx = self.main_model[pos[:,0],pos[:,1],pos[:,2]] != 0
        #
        # print 'idx sum', np.sum(idx)
        #
        # pos = pos[idx,:]
        # # pos = np.array([zz[idx], yy[idx], xx[idx]], dtype=np.int)

        print 'Calculating object coordinates....'

        # pos = np.zeros((zs*ys*xs,3), dtype= np.int)
        pos = np.array([], dtype= np.int)

        count = 0

        for z in xrange(0, zs):
            for y in xrange(0, ys):
                for x in xrange(0, xs):
                    if self.main_model[z,y,x] != 0:
                        # pos[count,:] = [z,y,x]
                        # count += 1
                        pos = np.append(pos, [z,y,x])

        pos = np.reshape(pos,(-1,3))

        t1 = tm.time()
        print 'Time taken %f seconds' % (t1-t0,)
        print 'getMainModelObjectPositions .... pos shape = ', pos.shape
        print 'Object coordinate calculation done!'
        return pos


    def getMainModelObjectPositions_under_development(self):

        zs, ys, xs = self.main_model.shape

        # t0 = tm.time()
        tic()
        print '****************model shape', self.main_model.shape

        print 'Calculating object coordinates....'


        # xx,yy,zz = np.meshgrid(np.linspace(0, xs-1,xs), np.linspace(0,ys-1,ys),\
        #                        np.linspace(0, zs-1, zs))
        #
        # L = zs*ys*xs
        #
        # zz = np.array(np.reshape(zz,(1,-1)), dtype=np.int)
        # yy = np.array(np.reshape(yy,(1,-1)), dtype=np.int)
        # xx = np.array(np.reshape(xx,(1,-1)), dtype=np.int)
        #
        # print 'xx yy zz ',
        # print xx.shape,
        # print yy.shape,
        # print zz.shape
        #
        # pos = np.reshape(np.array([zz,yy,xx]), (-1,3))
        #
        # idx = self.main_model[pos[:,0],pos[:,1],pos[:,2]] != 0
        #
        # print 'idx sum', np.sum(idx)
        #
        # pos = pos[idx,:]
        # pos = np.array([zz[idx], yy[idx], xx[idx]], dtype=np.int)

        # print 'pos shape  ', pos.shape
        # pos = np.reshape(pos, (-1,3))
        # print 'pos reshape  ', pos.shape


        pos = np.array([])

        for z in xrange(0, zs):
            for y in xrange(0, ys):
                for x in xrange(0, xs):
                    if self.main_model[z,y,x] != 0:
                        pos = np.append(pos, [z,y,x])

        pos = np.reshape(pos,(-1,3))

        t1 = toc()
        print 'Time taken %f seconds' % (t1,)
        print 'getMainModelObjectPositions .... pos shape = ', pos.shape
        print 'Object coordinate calculation done!'
        return pos


    def calculateContactSurfaces(self):

        print 'Calculating contact surface...'

        for pos in self.mainModelObjectCoords:
            self.updateContactNN(pos)

        print 'updateContactNN done!'

        print 'aircontactSurface === ', self.airContactSurface.keys()

        for k in self.airContactSurface.keys():
            self.airContactSurface[k].coords = np.reshape(self.airContactSurface[k].coords, (-1,3))

        print 'airContact done!'

        for k in self.materialContactSurface.keys():
            self.materialContactSurface[k].coords = np.reshape(self.materialContactSurface[k].coords, (-1,3))


        print 'materialContact done!'


    def updateContactNN(self, pos = (0,0,0)):

        zs, ys, xs = self.main_model.shape

        air_contact = self.airContactSurface
        mat_contact = self.materialContactSurface

        NN = np.array(pos + MKMultiMaterialSystem.NN_offset, dtype=np.int)
        # print 'NN == ', NN

        c_id = self.main_model[pos[0], pos[1], pos[2]]

        if np.prod(self.main_model[ NN[:,0], NN[:,1], NN[:,2]]) == 0:
            if c_id in air_contact.keys():
                air_contact[c_id].coords = np.append(air_contact[c_id].coords, pos)
            else:
                air_contact[c_id] = MKSurfaceContact()
                air_contact[c_id].id = 0
                air_contact[c_id].coords = np.array(pos, dtype= np.int)

        for p in NN:
            m_id = self.main_model[p[0], p[1], p[2]]
            if m_id == 0:
                continue

            if c_id != m_id:
                # if m_id in mat_contact.keys():
                #     if mat_contact[m_id] == c_id:
                #         continue

                if c_id in mat_contact.keys():
                    mat_contact[c_id].id = m_id
                    mat_contact[c_id].coords = np.append(mat_contact[c_id].coords, p)
                else:
                    mat_contact[c_id] = MKSurfaceContact()
                    mat_contact[c_id].id = m_id
                    mat_contact[c_id].coords = np.array(p, dtype=np.int)

    def getSurfaceTemperatureData(self):
        pass

    def setBoundaryCondition(self, mat1, mat2):
        pass


    def getD(self, pos):
        # d = self.
        pass

    def getAllAirContactSurfacePoints(self, retainData = 1.0):

        pos = np.array([])
        for k in self.airContactSurface.keys():
            pos = np.append(pos, self.airContactSurface[k].coords.reshape((1,-1)))

        pos = pos.reshape((-1,3))

        if retainData == 1:
            return pos

        r,c = pos.shape
        skip_val = int(1.0/retainData)

        pos_new = np.array([])

        for  idx in xrange(0, r):
            if np.random.ranf() <= retainData:
                pos_new = np.append(pos_new, pos[idx,:])

        return pos_new.reshape((-1,3))

    def getAllMaterialContactSurfacePoints(self, retainData = 1.0):

        pos = np.array([])
        for k in self.materialContactSurface.keys():
            pos = np.append(pos, self.materialContactSurface[k].coords.reshape((1,-1)))

        pos = pos.reshape((-1,3))

        if retainData == 1:
            return pos

        r,c = pos.shape
        skip_val = int(1.0/retainData)

        pos_new = np.array([])

        for  idx in xrange(0, r):
            if np.random.ranf() <= retainData:
                pos_new = np.append(pos_new, pos[idx,:])

        return pos_new.reshape((-1,3))


    def imposeDirichletBC(self, mat_id, wall = 0, bc_value = 0):

        coords = np.array(self.getMainModelObjectPositions(), dtype=np.int)
        idx = self.main_model[coords[:,0],coords[:,1],coords[:,2]] == mat_id

        if wall == 'zbottom' or wall == 'ztop':
            if wall == 'zbottom':
                zb = np.min(coords[idx, 0])
            else:
                zb = np.max(coords[idx, 0])

            idz = coords[idx, 0] == zb

        elif wall == 'ytop' or wall == 'ybottom':
            if wall == 'ybottom':
                zb = np.min(coords[idx, 1])
            else:
                zb = np.max(coords[idx, 1])

            idz = coords[idx, 1] == zb

        elif wall == 'xtop' or wall == 'xbottom':
            if wall == 'xbottom':
                zb = np.min(coords[idx, 2])
            else:
                zb = np.max(coords[idx, 2])

            idz = coords[idx, 2] == zb


        plane = np.array(coords[idz, :].reshape((1,-1)), dtype=np.int)

        mdata = MKBoundaryConditionData()
        mdata.coords = np.array(plane, dtype = np.int)
        mdata.mat_id = mat_id
        mdata.bcType = 0
        mdata.bcValue = bc_value

        if self.dirichletBCCoordinates is None:
            self.dirichletBCCoordinates =  [mdata]  #np.array(plane, dtype=np.int)
        else:
            self.dirichletBCCoordinates.append(mdata)
            # self.dirichletBCCoordinates[mat_id].coords = np.append(self.dirichletBCCoordinates[mat_id].coords, plane)

        for m_dt in self.dirichletBCCoordinates:
            m_dt.coords = m_dt.coords.reshape((-1,3))
            print 'self.dirichletBCCoordinates mat %d shape %s' % (m_dt.mat_id, m_dt.coords.shape, )


    def imposeNeumannBC(self, mat_id, wall = 0, bc_value = 0):

        coords = self.mainModelObjectCoords
        idx = self.main_model[coords[:,0],coords[:,1],coords[:,2]] == mat_id
        direction = 0

        if wall == 'zbottom' or wall == 'ztop':
            if wall == 'zbottom':
                zb = np.min(coords[idx, 0])
                direction = MKDirection.ZMIN
            else:
                zb = np.max(coords[idx, 0])
                direction = MKDirection.ZMAX

            idz = coords[idx, 0] == zb

        elif wall == 'ytop' or wall == 'ybottom':
            if wall == 'ybottom':
                zb = np.min(coords[idx, 1])
                direction = MKDirection.YMIN
            else:
                zb = np.max(coords[idx, 1])
                direction = MKDirection.YMAX

            idz = coords[idx, 1] == zb

        elif wall == 'xtop' or wall == 'xbottom':
            if wall == 'xbottom':
                zb = np.min(coords[idx, 2])
                direction = MKDirection.XMIN

            else:
                zb = np.max(coords[idx, 2])
                direction = MKDirection.XMAX

            idz = coords[idx, 2] == zb

        plane = np.array(coords[idz, :].reshape((1,-1)), dtype=np.int)

        # plane = plane.reshape((-1,3))

        mdata = MKBoundaryConditionData()
        mdata.coords = np.array(plane, dtype = np.int)
        mdata.mat_id = mat_id
        mdata.bcType = 1 # 1 for Neumann
        mdata.bcValue = bc_value
        mdata.fluxDirection = direction

        if self.neumanBCCoordinates is None:
            self.neumanBCCoordinates = [mdata]
        else:
            self.neumanBCCoordinates.append(mdata)

        for m_dt in self.neumanBCCoordinates:
            m_dt.coords = m_dt.coords.reshape((-1,3))
            print 'self.neumannBCCoordinates mat %d shape %s' % (m_dt.mat_id, m_dt.coords.shape, )


    def getMaterialCoords(self, mat_id = None, recalculateAll = False):

        if mat_id is not None:
            if recalculateAll == True:
                self.mainModelObjectCoords = self.getMainModelObjectPositions()

            p = np.array(self.mainModelObjectCoords, dtype=np.int)
            idx = self.main_model[p[:,0], p[:,1], p[:,2]] == mat_id

            return p[idx,:]

    def getMaterialBasePosition(self, mat_id = None):

        coords = self.getMaterialCoords(mat_id=mat_id)

        zmin = np.min(coords[:,0])

        idx = coords[:,0] == zmin

        xmin = np.min(coords[idx,2])
        ymin = np.min(coords[idx,1])

        return [zmin, ymin, xmin]


    def imposeHeatSource(self, heat_source_id = None, Q = 1, c = 1, rho = 1, fraction = 1, \
                         mat_id = None, positionOffset= (0,0,0), sourceType = 'zplane', \
                         radius = 1000000, heatSourceRange = [10000,10000,10000]):

        hs = MKHeatSourceData(heat_id = heat_source_id, mat_id=mat_id, Q = Q, c = c, rho = rho, fraction = fraction)


        base_pos = self.getMaterialBasePosition(mat_id=mat_id)

        # base_pos = np.array(self.model_objects[mat_id].position, dtype=np.int)

        print 'base_pos == %s', base_pos

        idx =  positionOffset != 0

        source_position = np.array(base_pos + np.array(positionOffset), dtype=np.int) #+ idx * self.global_offset

        print 'base_pos', base_pos
        print 'source_pos', source_position

        model_coords = self.getMaterialCoords(mat_id= mat_id)

        print 'model ID === %d coords shape === %s' % (mat_id, model_coords.shape, )

        if sourceType == 'zplane':
            idx = model_coords[:,0] == source_position[0]
            # print ''
        elif sourceType == 'yplane':
            idx = model_coords[:,1] == source_position[1]
        elif sourceType == 'xplane':
            idx = model_coords[:,2] == source_position[2]

        coords_vals =  model_coords[idx,:]


        cval = np.mean(coords_vals, axis=0)

        print 'coords_vals shape', coords_vals.shape
        print 'cvals', cval

        idx = np.abs(cval[0] - coords_vals[:,0]) <  heatSourceRange[0]
        coords_vals =  coords_vals[idx,:]

        idx = np.abs(cval[1] - coords_vals[:,1]) <  heatSourceRange[1]
        coords_vals =  coords_vals[idx,:]

        idx = np.abs(cval[2] - coords_vals[:,2]) <  heatSourceRange[2]
        coords_vals =  coords_vals[idx,:]

        zs,ys,xs = self.model_objects[mat_id].model.model.shape

        pixel_volume = zs*ys*xs

        hs.coords = coords_vals # model_coords[idx,:]

        # hs.coords = model_coords
        hs.surface_area_of_heat_source_plane_in_pixels = hs.coords.__len__()
        hs.volume_of_material_in_pixel = pixel_volume
        # hs.total_voxel_volume = hs.surface_area * (self.delta_pos)**3.0  # in cm3

        hs.H = Q/(rho*c)

        hs.H_per_pixel = (hs.H * 1.0) #/ hs.surface_area
        self.heatSource[heat_source_id] = hs

        print 'total_voxel_volume == ', hs.total_voxel_volume
        print 'surface area === ', hs.surface_area
        print 'H   == ', hs.H
        print 'H_per_pixel == ', hs.H_per_pixel
        print 'hs.coords.shape ', hs.coords.shape



    def updateNeumannBCArray(self, mat_id, bcType, bcValue, fluxDirection, coords):

        mdata = MKBoundaryConditionData()
        mdata.bcType = bcType
        mdata.mat_id = mat_id
        mdata.bcValue = bcValue
        mdata.fluxDirection = fluxDirection
        mdata.coords = coords.reshape((-1,3))

        print 'mat_id == ', mat_id
        print 'coords shape ', coords.shape

        # vv = raw_input('Press Enter to continue...')

        if self.neumanBCCoordinates is None:

            self.neumanBCCoordinates = [mdata]
        else:
            self.neumanBCCoordinates.append(mdata)

    def updateParameters(self):
        self.sxx = self.delta_time/self.delta_pos2
        self.syy = self.sxx
        self.szz = self.syy


    def imposeNeumannBCExcludinDirichlet(self, bc_value = 0):
        # for k, coords

        self.neumanBCCoordinates = None

        zs,ys,xs = self.main_model.shape
        zz, yy, xx = np.meshgrid(np.linspace(0, zs-1, zs), np.linspace(0, ys-1, ys), np.linspace(0, xs-1, xs))

        pos = np.array([])

        for k, data in self.airContactSurface.iteritems():
            pos = np.append(pos, data.coords.reshape((1,-1)))

        pos = pos.reshape((-1,3))

        print 'pos before filter ', pos.shape

        if self.dirichletBCCoordinates is not None:
            for diri_data in self.dirichletBCCoordinates:
                for c in diri_data.coords:
                    idx = np.logical_not(np.all(pos == c, axis = 1))
                    pos = pos[idx,:]

            # for k, data in self.materialContactSurface.iteritems():
            #     for c in data.coords:
            #         if c.__len__() != 3:
            #             continue
            #
            #         idx = np.logical_not(np.all(pos == c, axis = 1))
            #         pos = pos[idx,:]

            ### return pos

        pos = np.array(pos, dtype = np.int)
        # coords = pos


        for m_id, model_data in self.model_objects.iteritems():

            print 'model id  == ', m_id,
            print '  model_data  === ', model_data

            model = model_data.model
            # m_id = model.model_id

            print 'Calculating coords for model ', m_id

            idx = self.main_model[pos[:,0], pos[:,1],pos[:,2]] == m_id

            coords = pos[idx,:]

            print '************* mat_id', m_id,
            print 'coords shape', coords.shape

            # print 'z idx shape', idx.shape
            # raw_input('Enter to continue')

            dmax = np.max(coords[:, 0])
            dmin = np.min(coords[:, 0])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.ZMAX, \
                                      coords = coords[coords[:, 0] == dmax, :])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.ZMIN, \
                                      coords = coords[coords[:, 0] == dmin, :])


            dmax = np.max(coords[:, 1])
            dmin = np.min(coords[:, 1])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.YMAX, \
                                      coords = coords[coords[:, 1] == dmax, :])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.YMIN, \
                                      coords = coords[coords[:, 1] == dmin, :])


            dmax = np.max(coords[:, 2])
            dmin = np.min(coords[:, 2])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.XMAX, \
                                      coords = coords[coords[:, 2] == dmax, :])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.XMIN, \
                                      coords = coords[coords[:, 2] == dmin, :])



        print 'pos after filter ', pos.shape

        return pos


    def calculatesTemperatures(self, tempData):
        for model_id, mdata in self.model_objects.iteritems():
            tempData[model_id] = np.mean(self.u[self.main_model[:] == model_id])


    def nan_add(self, v1, v2):
        idx1 = np.isnan(v1)
        idx2 = np.isnan(v2)

        idx_surface = np.logical_xor(np.isnan(v1) , np.isnan(v2))
        idx_inside = np.logical_and(np.isnan(v1) , np.isnan(v2))

        v3 = v1 + v2

        # v3[idx_surface] = v1[idx_surface]
        v3[idx] = 0
        return  v3


    # @vectorize([float64(float64)])
    # @jit(nogil = True)
    @jit(nogil = True, cache =False)
    def sum_2nd_derivative_z(ui, d, nan_flag):
        d[:,:,:] = ui[2:, 1:-1, 1:-1] + ui[:-2, 1:-1, 1:-1] - 2 * ui[1:-1, 1:-1, 1:-1]
        d[nan_flag] = 0

    # @vectorize([float64(float64)])
    # @jit(nogil = True)
    @jit(nogil = True, cache =False)
    def sum_2nd_derivative_y(ui,d, nan_flag):
        d[:,:,:] = ui[1:-1, 2:, 1:-1] + ui[1:-1, :-2, 1:-1] - 2 * ui[1:-1, 1:-1, 1:-1]
        d[nan_flag] = 0

    # @vectorize([float64(float64)])
    # @jit(nogil = True)
    @jit(nogil = True, cache = False)
    def sum_2nd_derivative_x(ui, d, nan_flag):
        d[:,:,:] = ui[1:-1, 1:-1, 2:] + ui[1:-1, 1:-1, :-2] - 2 * ui[1:-1, 1:-1, 1:-1]
        d[nan_flag] = 0

    @jit(nogil = False, cache = False)
    def sum_2nd_derivative_zyx(ui):
        return (ui[2:, 1:-1, 1:-1] + ui[:-2, 1:-1, 1:-1] - 2 * ui[1:-1, 1:-1, 1:-1], \
                ui[1:-1, 2:, 1:-1] + ui[1:-1, :-2, 1:-1] - 2 * ui[1:-1, 1:-1, 1:-1], \
                ui[1:-1, 1:-1, 2:] + ui[1:-1, 1:-1, :-2] - 2 * ui[1:-1, 1:-1, 1:-1]
                )


    def evolveTempWithTime_v5p3(self, enableDiri = True, enableNeumann = True, data = None, deltaTimeFactor = 1):

        self.u[:] = self.ui[:]

        pos = self.mainModelObjectCoords


        def mt_func(ui,pos,ax):
            if ax == 0:
                return (ui[pos[:,0]+1, pos[:,1], pos[:,2]] + ui[pos[:,0]-1, pos[:,1], pos[:,2]] - 2 * ui[pos[:,0], pos[:,1], pos[:,2]])
            if ax == 1:
                return (ui[pos[:,0], pos[:,1]+1, pos[:,2]] + ui[pos[:,0], pos[:,1]-1, pos[:,2]] - 2 * ui[pos[:,0], pos[:,1], pos[:,2]])
            if ax == 2:
                return (ui[pos[:,0], pos[:,1], pos[:,2]+1] + ui[pos[:,0], pos[:,1], pos[:,2]-1] - 2 * ui[pos[:,0], pos[:,1], pos[:,2]])

        # self.u[pos[:,0], pos[:,1], pos[:,2]] += self.D[pos]

        d1 = self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] + self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]
        d1[self.zfilter] = 0

        d2 =  self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] + self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]
        d2[self.yfilter] = 0

        d3 = self.ui[pos[:,0], pos[:,1], pos[:,2]+1] + self.ui[pos[:,0], pos[:,1], pos[:,2]-1] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]
        d3[self.xfilter] = 0

        # d1 = MKMultiMaterialSystem.sum_2nd_derivative_z(self.ui)
        # d2 = MKMultiMaterialSystem.sum_2nd_derivative_y(self.ui)
        # d3 = MKMultiMaterialSystem.sum_2nd_derivative_x(self.ui)

        self.u[pos[:,0], pos[:,1], pos[:,2]] += self.D[pos[:,0], pos[:,1], pos[:,2]] * self.delta_time / self.delta_pos2 * \
                                                ( d1 + d2 +  d3)


        d1 = None
        d2 = None
        d3 = None

        if self.externalHeatSourceActive == True:
            for heat_id, heat_source in self.heatSource.iteritems():
                pos = heat_source.coords
                # delta_heat = heat_source.H_per_pixel * self.delta_time
                delta_heat = heat_source.H * self.delta_time
                self.u[pos[:,0], pos[:,1], pos[:,2]] = self.u[pos[:,0], pos[:,1], pos[:,2]] + delta_heat

                # print 'H  %f   H * dt %.3e  H_per_pixel = %.4e  H_per_pixel * dt = %.4e  ' % (heat_source.H, heat_source.H * self.delta_time, heat_source.H_per_pixel, heat_source.H_per_pixel*self.delta_time, )
                # print 'delta_heat == ', delta_heat
                # print ' heat_source.H_per_pixel ', heat_source.H_per_pixel


        # print 'enableNeumann ', enableNeumann
        # print 'enableDiri', enableDiri
        # self.u[1:-1,1:-1,1:-1] = nx.evaluate(self.thermal_update)

        # Dirichlet bundary condition with air contact surface

        if enableDiri == True and self.dirichletBCCoordinates is not None:
            for d_data in self.dirichletBCCoordinates:
                pos = d_data.coords
                self.u[pos[:,0], pos[:,1], pos[:, 2]] = d_data.bcValue

        # Neumann boundary condition along z-axis

        if enableNeumann == True and self.neumanBCCoordinates is not None:
            # if self.neumanBCCoordinates is None:

            for ndata in self.neumanBCCoordinates:
                pos = ndata.coords
                mat_id = ndata.mat_id
                bc_val = ndata.bcValue

                # print 'entered neumann loop'
                # print 'mat_id == ', mat_id

                dfactor = 2* self.D_dict[mat_id] * self.delta_time / self.delta_pos2

                # print 'Neumann bc value ', ndata.bcValue
                # print 'mat_id %d    bcValue %d  flux   ' % (ndata.mat_id, ndata.bcValue, ndata.fluxDirection,)


                if ndata.fluxDirection == MKDirection.ZMAX:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] + self.delta_pos * bc_val) )
                    # print 'ZMAX'

                elif ndata.fluxDirection == MKDirection.ZMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - self.delta_pos * bc_val) )
                    # print 'ZMIN'

                elif ndata.fluxDirection == MKDirection.YMAX:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] + self.delta_pos * bc_val) )
                    # print 'YMAX'
                elif ndata.fluxDirection == MKDirection.YMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] - self.delta_pos * bc_val) )
                    # print 'YMIN'

                elif ndata.fluxDirection == MKDirection.XMAX:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]-1] + self.delta_pos * bc_val) )
                    # print 'XMAX'

                elif ndata.fluxDirection == MKDirection.XMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] - self.delta_pos * bc_val) )
                    # print 'XMIN'


        self.ui[:] = self.u[:] #sp.copy(self.u)


    def evolveTempWithTime_v5p2(self, enableDiri = True, enableNeumann = True, data = None, deltaTimeFactor = 1, enableHeatDecay = True):

        # MKMultiMaterialSystem.sum_2nd_derivative_z(self.ui, self.d1, self.zfmatrix)
        # MKMultiMaterialSystem.sum_2nd_derivative_y(self.ui, self.d2, self.yfmatrix)
        # MKMultiMaterialSystem.sum_2nd_derivative_x(self.ui, self.d3, self.xfmatrix)

        self.d1[:,:,:] = self.ui[2:, 1:-1, 1:-1] + self.ui[:-2, 1:-1, 1:-1] - 2.0 * self.ui[1:-1, 1:-1, 1:-1]
        self.d2[:,:,:] = self.ui[1:-1, 2:, 1:-1] + self.ui[1:-1, :-2, 1:-1] - 2.0 * self.ui[1:-1, 1:-1, 1:-1]
        self.d3[:,:,:] = self.ui[1:-1, 1:-1, 2:] + self.ui[1:-1, 1:-1, :-2] - 2.0 * self.ui[1:-1, 1:-1, 1:-1]

        self.d1[self.zfmatrix] = 0.0
        self.d2[self.yfmatrix] = 0.0
        self.d3[self.xfmatrix] = 0.0

        self.u[1:-1, 1:-1, 1:-1] = self.ui[1:-1, 1:-1, 1:-1] + self.D[1:-1,1:-1,1:-1] * self.delta_time / self.delta_pos2 * \
                                  ( self.d1 + self.d2 +  self.d3 )

        if self.externalHeatSourceActive == True:
            for i in xrange(0,self.heatSource.__len__()):

                heat_id = self.heatSource.keys()[i]
                heat_source = self.heatSource.values()[i]

                pos = heat_source.coords
                # pos = np.array(self.model_objects[3].position, dtype=np.int)

                delta_heat = float( heat_source.H / self.heatSourceVolume )

                # print 'delta_heat === ', delta_heat
                # print 'delta_heat*dt === ', delta_heat*self.delta_time


                self.u[pos[:,0], pos[:,1], pos[:,2]] = self.ui[pos[:,0], pos[:,1], pos[:,2]] + delta_heat * self.delta_time

                # print 'H  %f   H * dt %.3e  H_per_pixel = %.4e  H_per_pixel * dt = %.4e  ' % (heat_source.H, heat_source.H * self.delta_time, heat_source.H_per_pixel, heat_source.H_per_pixel*self.delta_time, )
                # print 'delta_heat == ', delta_heat
                # print ' heat_source.H_per_pixel ', heat_source.H_per_pixel


        # print 'enableNeumann ', enableNeumann
        # print 'enableDiri', enableDiri
        # self.u[1:-1,1:-1,1:-1] = nx.evaluate(self.thermal_update)

        # Neumann boundary condition along z-axis

        if self.neumanBCCoordinates is not None:
            if enableNeumann == True:
                # if self.neumanBCCoordinates is None:

                for ndata in self.neumanBCCoordinates:
                    pos = ndata.coords
                    mat_id = ndata.mat_id
                    bc_val = ndata.bcValue

                    # print 'entered neumann loop'
                    # print 'mat_id == ', mat_id

                    dfactor = 2* self.D_dict[mat_id] * self.delta_time / self.delta_pos2

                    # print 'Neumann bc value ', ndata.bcValue
                    # print 'mat_id %d    bcValue %d  flux   ' % (ndata.mat_id, ndata.bcValue, ndata.fluxDirection,)


                    if ndata.fluxDirection == MKDirection.ZMAX:
                        self.u[pos[:,0], pos[:,1], pos[:,2]] = self.ui[pos[:,0], pos[:,1], pos[:,2]] + -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] + self.delta_pos * -bc_val) )
                        # print 'ZMAX'

                    elif ndata.fluxDirection == MKDirection.ZMIN:
                        self.u[pos[:,0], pos[:,1], pos[:,2]] = self.ui[pos[:,0], pos[:,1], pos[:,2]] + -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - self.delta_pos * bc_val) )
                        # print 'ZMIN'

                    elif ndata.fluxDirection == MKDirection.YMAX:
                        self.u[pos[:,0], pos[:,1], pos[:,2]] = self.ui[pos[:,0], pos[:,1], pos[:,2]] + -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] + self.delta_pos * -bc_val) )
                        # print 'YMAX'
                    elif ndata.fluxDirection == MKDirection.YMIN:
                        self.u[pos[:,0], pos[:,1], pos[:,2]] = self.ui[pos[:,0], pos[:,1], pos[:,2]] + -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] - self.delta_pos * bc_val) )
                        # print 'YMIN'

                    elif ndata.fluxDirection == MKDirection.XMAX:
                        self.u[pos[:,0], pos[:,1], pos[:,2]] = self.ui[pos[:,0], pos[:,1], pos[:,2]] + -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]-1] + self.delta_pos * -bc_val) )
                        # print 'XMAX'

                    elif ndata.fluxDirection == MKDirection.XMIN:
                        self.u[pos[:,0], pos[:,1], pos[:,2]] = self.ui[pos[:,0], pos[:,1], pos[:,2]] + -dfactor * (self.ui[pos[:,0], pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] - self.delta_pos * bc_val) )
                        # print 'XMIN'

            if enableHeatDecay == True:
                for ndata in self.neumanBCCoordinates:
                    pos = ndata.coords
                    mat_id = ndata.mat_id
                    bc_val = ndata.bcValue

                    RhoCp = self.model_objects[mat_id].model.rho * self.model_objects[mat_id].model.Cp

                    v1 =  50e-4 * self.voxelArea * (self.ui[pos[:,0], pos[:,1], pos[:,2]]- 25)
                    v2 =  0.5 * 5.67e-12 *self.voxelArea * ((self.ui[pos[:,0], pos[:,1], pos[:,2]]+273)**4 - 298**4)

                    # val =  0.5 * 5.67e-12 *self.voxelArea * ( (self.ui[pos[:,0], pos[:,1], pos[:,2]]+273)**4 - 298**4)

                    # print 'RhoCp', RhoCp
                    # print '*********************v1 max ', np.max(v1), 'temp maxx == ', np.max(self.ui[pos[:,0], pos[:,1], pos[:,2]])
                    # print '************v2 max ', np.max(v2)

                    # print '********************************************* dT === ', val

                    self.u[pos[:,0], pos[:,1], pos[:,2]] -= (v1+v2) / RhoCp * self.delta_time
                    v1 = None
                    v2 = None



        # Dirichlet bundary condition with air contact surface

        if enableDiri == True and self.dirichletBCCoordinates is not None:
            for d_data in self.dirichletBCCoordinates:
                pos = d_data.coords
                self.u[pos[:,0], pos[:,1], pos[:, 2]] = d_data.bcValue


        self.ui[:] = self.u[:] #sp.copy(self.u)


    def calculateNANIndicesAndTempStorage(self):

        pos = np.array(self.getMainModelObjectPositions(), dtype = np.int)

        d =  self.ui[pos[:,0]+1,pos[:,1],pos[:,2]] + self.ui[pos[:,0]-1,pos[:,1],pos[:,2]]
        self.zfilter = np.isnan( d )


        d =  self.ui[pos[:,0],pos[:,1]+1,pos[:,2]] + self.ui[pos[:,0],pos[:,1]+1,pos[:,2]]
        self.yfilter = np.isnan( d )

        d = self.ui[pos[:,0],pos[:,1],pos[:,2]+1] + self.ui[pos[:,0],pos[:,1],pos[:,2]-1]
        self.xfilter = np.isnan( d )

        print 'self.zfilter   ', self.zfilter.shape, '  sum ', np.sum(self.zfilter)
        print 'self.yfilter   ', self.yfilter.shape, '  sum ', np.sum(self.yfilter)
        print 'self.xfilter   ', self.xfilter.shape, '  sum ', np.sum(self.xfilter)


        self.zfmatrix = np.isnan(self.ui[2:,1:-1,1:-1] + self.ui[:-2,1:-1,1:-1])
        self.yfmatrix = np.isnan(self.ui[1:-1,2:,1:-1] + self.ui[1:-1,:-2,1:-1])
        self.xfmatrix = np.isnan(self.ui[1:-1,1:-1,2:] + self.ui[1:-1,1:-1,:-2])

        self.d1 = np.zeros(self.zfmatrix.shape, dtype = np.float)
        self.d2 = np.zeros(self.yfmatrix.shape, dtype = np.float)
        self.d3 = np.zeros(self.xfmatrix.shape, dtype = np.float)

        self.mat_index_dict = {}

        for mat_id, d_data in self.D_dict.iteritems():
            self.mat_index_dict[mat_id] = self.main_model[:] == mat_id


        pass



    # def startSimulation(self, maxTimeIteration = 1000):
    #
    #     for i in xrange(maxTimeIteration):
    #         self.evolveT()


    def plot2DSectionData(self, fig_no, data2D):

        f = plt.figure(fig_id)
        plt.clf()
        plt.cla()
        im = plt.imshow(data2D)
        plt.colorbar(im)


    def plotSurface3D(self):
        pass


    def analylicSolutionPointSource(u_init, t, D = 1, pos0 = (0,0,0), t0 = 0):
        c = 1/(4.0 * np.pi * (t - t0) * D)**(3.0/2.0) * np.exp( ((pos[:,2] - pos0[2])**2 +  (pos[:,1] - pos0[1])**2 + (pos[:,0] - pos0[0])**2)/(4*t*D)      )
        return c




    def runSimulation(self, iterationTime = 1000, initU = None, figureID = 666, ySection = 15, \
                      applyDiri = True, applyNeumann = True, deltaTimeFactor = 1, toggleHeating=False, toggleDuration=1000):
        ######################
        self.ui[:] = initU[:]

        data = []

        time_steps = iterationTime
        count = 0
        tic()

        temp_dict = {}
        for mid, d_val in self.D_dict.iteritems():
            temp_dict[mid] = np.array([])

        time_data = np.array([])

        for i in xrange(time_steps):
            self.evolveTempWithTime_v5p1(enableDiri=applyDiri, enableNeumann=applyNeumann)

            for mid, d_val in self.D_dict.iteritems():
                p = self.model_objects[mid].position
                zs,ys,xs = self.model_objects[mid].model.model.shape

                pz = p[0]
                py = p[1]
                px = p[2]

                pzm = pz + zs
                pym = py + ys
                pxm = px + xs

                temp_dict[mid] = np.append(temp_dict[mid], np.nanmean(self.u[pz:pzm, py:pym, px:pxm]))

            time_data = np.append(time_data, i * self.delta_time)

            progress =  np.ceil(float(i)/time_steps*100)

            if progress >= count*10:
                count += 1
                print 'Simulation completed %d' % (progress,),
                print '%'

        plot2DImage(figureID, self.u[:,ySection,:], xAxisLabel= 'X (mm)', yAxisLabel= 'Y (mm)')

        print 'Total simulation time = %.4e seconds with time steps of %f ms' % (time_steps * self.delta_time, self.delta_time*1000.0,)
        t2 =  toc()

        print 'Total time taken ',

        t2 = toc()
        if t2 < 60:
            print t2, 'seconds'
        else:
            print t2/60.0, 'minutes'

        return [time_data, temp_dict, t2]



    def runSimulation_v2(self, iterationTime = 1000, initU = None, figureID = 666, ySection = 15, \
                         applyDiri = True, applyNeumann = True, applyHeatDecay = True,\
                         deltaTime = None, \
                         heatSourceVolume = 1,\
                         voxelL = 10e-6*100.0,\
                         shutDownExternalHeatingAtTime = None, shutDownDuration = None):
        ######################


        if deltaTime is None:
            print 'deltaTime not given!'

        # if deltaPos is None:
        #     print 'deltaPos (in cm) not given'

        self.voxelL = voxelL                # cm
        self.voxelArea = self.voxelL**2     # cm2
        self.voxelVolume = self.voxelL**3   # cm3

        self.delta_time = deltaTime

        self.heatSourceVolume = heatSourceVolume * self.voxelVolume# in cm3

        self.delta_pos = self.voxelL     # cm to m
        self.delta_pos2 = self.delta_pos**2.0

        self.updateParameters()
        self.calculateNANIndicesAndTempStorage()

        self.u[:] = initU[:]
        self.ui[:] = initU[:]
        # data = []


        zs,ys,xs = self.main_model.shape
        mid_point = np.array([zs, ys, xs])/2 -1

        self.simulation_result = []

        print 'mid_point ', mid_point

        mid_temp_data = np.array([])

        # time_steps = iterationTime
        count = 0
        tic()

        temp_dict = {}
        for mid, d_val in self.D_dict.iteritems():
            temp_dict[mid] = np.array([])

        time_data = np.array([], dtype=np.float)


        shutdown_time = 1e100

        if shutDownExternalHeatingAtTime is not None:
            shutdown_time = shutDownExternalHeatingAtTime


        for mid, data in self.D_dict.iteritems():
            temp_dict[mid] = np.zeros((iterationTime, 4), dtype = np.float)

        mid_temp_data = np.zeros((iterationTime,))
        time_data = np.zeros((iterationTime,))

        it_time = tm.time()


        status_update_interval = 10

        if iterationTime > 10000000L:
            status_update_interval = 1
        elif iterationTime > 1000000L:
            status_update_interval = 2
        elif iterationTime > 10000L:
            status_update_interval = 5


        for i in xrange(iterationTime):
            ####################

            if self.externalHeatSourceActive == True:
                if self.delta_time *i > shutdown_time:
                    self.externalHeatSourceActive = False




            # self.evolveTempWithTime_v5p2(enableDiri=applyDiri, enableNeumann=applyNeumann)
            self.evolveTempWithTime_v5p2(enableDiri=applyDiri, enableNeumann=applyNeumann, enableHeatDecay = applyHeatDecay)





            mid_temp = self.u[mid_point[0], mid_point[1], mid_point[2]]
            mid_temp_data[i] =  mid_temp
            ###################

            for jdx  in xrange(self.D_dict.__len__()):
                mid = self.D_dict.keys()[jdx]
                d_val = self.D_dict.values()[jdx]

                p = self.model_objects[mid].position
                zs,ys,xs = self.model_objects[mid].model.model.shape

                pz = p[0]
                py = p[1]
                px = p[2]

                pz_section = zs/3

                pzm = pz + zs
                pym = py + ys
                pxm = px + xs

                # v = self.u[pz:pzm, py:pym, px:pxm]
                # print 'zs ys xs p pmax ', zs,ys, xs, p, (pzm, pym, pxm)

                # temp = np.nanmean(v)
                # temp = np.nanmean(v)

                temp = np.mean(self.u[self.mat_index_dict[mid]])
                temp_bottom = np.mean(self.u[pz:(pz+pz_section), py:pym, px:pxm])
                temp_middle = np.mean(self.u[(pz+pz_section):(pz+2*pz_section), py:pym, px:pxm])
                temp_top =    np.mean(self.u[(pz+2*pz_section):pzm, py:pym, px:pxm])

                # print 'time %.3e mean_temp %.4f ' % ( i * self.delta_time, temp, )
                # print 'temp == ', temp

                temp_dict[mid][i,:] = [temp, temp_top, temp_middle, temp_bottom ]

            time_data[i] = i * self.delta_time

            progress =  np.ceil(float(i)/iterationTime*100.0)

            if progress >= count:
                gt = tm.time()
                gc.collect()
                print 'Time (sec) taken for gc ==  ', tm.time() - gt

                print 'Simulation completed %d' % (progress,),
                print '%'
                print 'calculation took time (sec) ', tm.time() - it_time
                it_time = tm.time()

                count += status_update_interval
                # print 'count == ', count


        t2 =  toc()


        if t2 < 60:
            t2_total = t2
            tstr = '%.1e sec' % t2_total,
        else:
            t2_total =  t2/60.0
            tstr = '%.2f minutes' % t2_total,


        result_str = '''Total simulation time = %4.0e sec (dt %1.0e sec)\nComputation Time %s''' % (iterationTime * self.delta_time, self.delta_time, tstr,)
        print result_str


        plot2DImage(figureID, self.u[:,:,ySection], xAxisLabel= 'X (10x um)', yAxisLabel= 'Z (10x um)', \
                    title_str=result_str, interpolationMethod='none')




        self.externalHeatSourceActive = True

        return [time_data, temp_dict, t2], mid_temp_data



    def showAllMaterialTemperature(self):

        # if self.sim

        pass


    def showAllNeumannBCWalls(self, figure_id = 19001, xlabel = 'X (mm)', ylabel = 'Y (mm)', zlabel = 'Z (mm)', \
                              title_str = 'Neumann BC', ax = None, \
                              xlim = None, ylim = None, zlim = None):

        if self.neumanBCCoordinates is None:
            return

        len = self.neumanBCCoordinates.__len__()


        colors = cm.rainbow(np.linspace(0,1,len))

        cc = itertools.cycle(colors)

        if ax is None:
            f3 = plt.figure(figure_id)
            plt.clf()
            plt.cla()
            ax = f3.add_subplot(111, projection = '3d')

        for d in self.neumanBCCoordinates:
            pos = d.coords
            print 'material id  %d    direction %d' % (d.mat_id, d.fluxDirection,)
            # ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = next(cc))

            x = pos[:,2].flatten()
            y = pos[:,1].flatten()
            z = pos[:,0].flatten()

            # ax.scatter(x, y, z, marker = '.', color = next(cc))

            if x.__len__() <=6 or x.__len__() <=6 or x.__len__() <=6:
                ax.scatter(x, y, z, marker = '.', color = next(cc))
            else:
                # try:
                #     ax.plot_trisurf(x, y, z, color = next(cc), alpha = 0.6)
                # except:
                ax.scatter(x, y, z, marker = '.', color = next(cc))

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_zlabel(zlabel)
            # plt.title(title_str)

        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)
        if zlim is not None:
            plt.zlim(zlim)

        plt.title(title_str)
        plt.draw()

    def showAllDirichletBCWalls(self, figure_id = 19002, xlabel = 'X (mm)', ylabel = 'Y (mm)', zlabel = 'Z (mm)', \
                                title_str = 'Dirichlet BC', ax = None, \
                                xlim = None, ylim = None, zlim = None):

        if self.dirichletBCCoordinates is None:
            return

        len = self.dirichletBCCoordinates.__len__()
        colors = cm.rainbow(np.linspace(0,1,len))
        cc = itertools.cycle(colors)

        if ax is None:
            f3 = plt.figure(figure_id)
            plt.clf()
            plt.cla()
            ax = f3.add_subplot(111, projection = '3d')

        # ax.invert_zaxis()

        for d in self.dirichletBCCoordinates:
            pos = d.coords
            print 'material id  %d    direction %d' % (d.mat_id, d.fluxDirection,)
            # ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = next(cc))
            x = pos[:,2].flatten()
            y = pos[:,1].flatten()
            z = pos[:,0].flatten()

            # ax.scatter(x, y, z, marker = '.', color = next(cc))

            if x.__len__() <=3 or x.__len__() <=3 or x.__len__() <=3:
                ax.scatter(x, y, z, marker = '.', color = next(cc))
            else:
                # try:
                #     ax.plot_trisurf(x, y, z, color = next(cc), alpha = 0.6, linewidth=0.01)
                # except:
                ax.scatter(x, y, z, marker = '.', color = next(cc))

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_zlabel(zlabel)
            # plt.title(title_str)

        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)
        if zlim is not None:
            plt.zlim(zlim)

        plt.title(title_str)
        plt.draw()


    def showAllHeatSources(self, figure_id = 88002, xlabel = 'X (mm)', ylabel = 'Y (mm)', zlabel = 'Z (mm)', \
                           title_str='Heat Sources', ax = None, \
                           xlim=None, ylim=None,zlim=None):

        if self.heatSource is None:
            return

        len = self.heatSource.__len__()
        colors = cm.rainbow(np.linspace(0, 1, len))
        cc = itertools.cycle(colors)

        if ax is None:
            f3 = plt.figure(figure_id)
            plt.clf()
            plt.cla()
            ax = f3.add_subplot(111, projection = '3d')

        for h_id, h in self.heatSource.iteritems():
            pos = h.coords
            print 'heat source id  %d    direction %.4f' % (h.heat_id, h.H,)

            # ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = next(cc))
            x = pos[:,2].flatten()
            y = pos[:,1].flatten()
            z = pos[:,0].flatten()

            # ax.scatter(x, y, z, marker = '.', color = next(cc))

            # if x.__len__() <=3 or x.__len__() <=3 or x.__len__() <=3:
            #     ax.scatter(x, y, z, marker = '.', color = next(cc))
            # else:
            ax.plot_trisurf(x, y, z, color = next(cc), alpha = 0.6)

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_zlabel(zlabel)
            # plt.title(title_str)


        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)
        if zlim is not None:
            plt.zlim(zlim)

        plt.title(title_str)
        plt.draw()


    def showAllMaterialContactWalls(self, figure_id = 19000, xlabel = 'X (mm)', ylabel = 'Y (mm)', zlabel = 'Z (mm)', title_str = '', \
                                    xlim = None, ylim = None, zlim = None, ax = None):

        pos = self.getAllMaterialContactSurfacePoints(retainData=0.4)

        if ax is None:
            f3 = plt.figure(figure_id)
            plt.clf()
            plt.cla()
            ax = f3.add_subplot(111, projection = '3d')

        # ax.scatter(pos[:,2], pos[:,1], pos[:,0], marker = '.', color = 'red')

        x = pos[:,2].flatten()
        y = pos[:,1].flatten()
        z = pos[:,0].flatten()

        ax.scatter(x, y, z, marker = '.', color = 'red')
        # if x.__len__() <=3 or x.__len__() <=3 or x.__len__() <=3:
        #     ax.scatter(x, y, z, marker = '.', color = 'red')
        # else:
        #     ax.plot_trisurf(x, y, z, color = 'red', alpha = 0.6)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)




        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)
        if zlim is not None:
            plt.zlim(zlim)

        plt.title(title_str)
        plt.draw()



    def showAllWalls(self, figure_id = 19000,  title_str = 'Diri, Neumann , Heat And MAterial Contact Walls'):

        len = 0

        if self.dirichletBCCoordinates is not None:
            len += self.dirichletBCCoordinates.__len__()

        if self.neumanBCCoordinates is not None:
            len += self.neumanBCCoordinates.__len__()

        if self.heatSource is not None:
            len += self.heatSource.__len__()

        if self.materialContactSurface is not None:
            len += self.materialContactSurface.__len__()

        colors = cm.rainbow(np.linspace(0,1,len))
        cc = itertools.cycle(colors)


        plt.figure(figure_id)
        plt.clf()
        plt.cla()
        ax = plt.gca(projection='3d')


        self.showAllDirichletBCWalls(ax =ax)
        self.showAllNeumannBCWalls(ax=ax)
        self.showAllHeatSources(ax=ax, title_str='')

        plt.suptitle(title_str)
        plt.show()



        pass


def plot2DImage(fig_id, image0, doInvertY = True, xAxisLabel = 'X (mm)', yAxisLabel= 'Y (mm)', title_str = '',\
                showContour = True, interpolationMethod = 'none'):
    ff = plt.figure(fig_id)
    plt.clf()
    plt.cla()
    ax = ff.add_subplot(111)
    img = np.round(image0, 7)

    im = ax.imshow(img, cmap=cm.coolwarm, alpha = 0.5, interpolation = interpolationMethod)

    if doInvertY:
        ax.invert_yaxis()


    if showContour:
        ys,xs = img.shape

        xx,yy = np.meshgrid(np.linspace(0, xs-1,xs)*1, np.linspace(0,ys-1,ys)*1)

        cs = plt.contour(xx,yy,img, linewidths=1)
        plt.clabel(cs, inline = 1, fontsize = 10)


    plt.colorbar(im, cmap = cm.coolwarm)
    # colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
    # cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])# horizontal colorbar

    # vmin = np.min(image0)
    # vmax = np.max(image0)
    # vmean = np.mean(image0)
    #
    # ts = [vmin, vmean, vmax]
    #
    # vals = '%.4f '
    #
    # s1 = '%.4f' % vmin
    # s2 = '%.4f' % vmean
    # s3 = '%.4f' % vmax
    #
    #
    # cb = ff.colorbar(im, orientation = 'vertical')
    # cb.ax.set_yticklabels([s1,s2,s3])

    ax.set_xlabel(xAxisLabel)
    ax.set_ylabel(yAxisLabel)
    plt.title(title_str)



def plot2DSurfaceAndContour(X, Y, Z, fig = None, figure_id = 1001, xlim = None, ylim = None, zlim = None, \
                            xlabel = 'X', ylabel = 'Y', zlabel = 'Z', cmap = cm.coolwarm , \
                            cticks = None, cbarLabel = '', rstride = 4, cstride = 4, title_str = ''):

    if fig is None:
        fig = plt.figure(figure_id)
        # plt.ioff()
        # plt.clf()
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


    ph =ax.plot_surface(X, Y, Z, rstride=rstride, cstride=cstride, alpha=0.5, cmap= cmap)
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
    plt.title(title_str)

    # plt.colorbar(ph, ticks = ct)

    fig.canvas.draw()
    tm.sleep((0.01))

    return  None



tm_val = 0

def tic():
    global tm_val

    tm_val = tm.time()

def toc():
    global tm_val
    return tm.time() - tm_val


def TestExponentialFit(figure_id = 99, x= [], y = [], doplot = False):

    def func1(x, a, b):
        return a * np.exp(-b * x)

    popt, pcov = curve_fit(func1, x, y, p0=(1, 1e-6))

    h = 0

    if doplot == True:
        plt.figure(figure_id)
        plt.clf()
        plt.cla()
        hdata, = plt.plot(x,y, color='red', marker = 'o', label = 'Data')
        hfit, = plt.plot(x, func1(x, *popt), color='blue', label = 'Exp Fit')
        plt.legend(handles = [hdata, hfit])
        title_str  = 'y = a*exp(-b * x) + c, (a,b,c) = ', popt
        plt.title(title_str)
        h = plt.show()

    return popt, func1, h


def saveObject(fname, obj_to_save):

    with open(fname, 'wb') as fptr:
        pp = pickle.Pickler(fptr, pickle.HIGHEST_PROTOCOL)
        pp.dump(obj_to_save)



def loadObjects(fname):

    with open(fname, 'rb') as fptr:
        pp = pickle.Pickler(fptr, pickle.HIGHEST_PROTOCOL)
        return pp.load()