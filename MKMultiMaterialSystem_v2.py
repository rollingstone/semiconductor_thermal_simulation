
__author__ = 'kamal'


import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from matplotlib.tri import triangulation
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
import time as tm
import multiprocessing.pool as pool
import numexpr as nx
from numba import  vectorize, jit, float32, float64, int32, int64



from copy import copy


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

        x = position[2]
        y = position[1]
        z = position[0]

        zs,ys,xs = model.model.shape
        zz,yy,xx = np.meshgrid(np.linspace(zs, z+zs-1, zs), np.linspace(ys, y+ys-1, ys), np.linspace(xs, x+xs-1, xs))

        zz = zz.reshape((1,-1))
        yy = yy.reshape((1,-1))
        xx = xx.reshape((1,-1))

        self.coords = np.array(np.append(np.append(zz,yy), xx).reshape((-1,3)), dtype= np.int)

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

class MKModel(object):
    def __init__(self, model_id = -1,  model_size = (10,10), D = 0.5, initial_condition = 1, dirichlet_bc = 1, neumann_bc = 1):

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

    def plotModelSurface(self, fig_id, filter_out = 0.8):
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

    def __init__(self, base_model,  max_materials = 10, delta_pos = (0.01,0.01, 0.01), time_step_scaling = 1, padding = (0.2, 0.2, 0.2)):

        np.set_printoptions(precision=8)

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

        # self.model_objects.append(base_model)
        self.model_position = []

        self.all_bc_coords = []
        self.all_dirichlet_value = []
        self.all_neumann_value = []


        new_size = np.round(self.base_model.model_size  * (1.0 + np.array(padding)))
        mx_val = np.max(new_size)
        new_size = np.array([mx_val, mx_val, mx_val])

        self.main_model = np.zeros(new_size, dtype=np.byte);
        self.ui = np.zeros(new_size, dtype=np.float64) + np.nan
        self.u = np.zeros(new_size, dtype=np.float64) + np.nan

        self.u_neumann = np.zeros(new_size, dtype=np.float)
        self.D = np.zeros(new_size, dtype=np.float) + np.nan
        # self.D = np.ones(new_size, dtype=np.float)
        self.D_dict = {}

        self.mainModelObjectCoords = [] ## save all the object coords

        self.airContactSurface = {}
        self.materialContactSurface = {}
        self.materialTemperature = {}

        self.modelDirichletBoundaryConditions = np.zeros(new_size, dtype= float)
        self.modelNeumannBoundaryConditions = np.zeros(new_size, dtype= float)

        self.sxx = 1.0
        self.syy = 1.0
        self.szz = 1.0

        self.max_D = -9999999
        self.min_deltaPosition = np.min(delta_pos)


        center_pos = new_size / 2;
        # self.base_model.model_size  * (1.0 + np.array(padding))
        offset_pos = center_pos - self.base_model.model_size/2;
        self.model_position.append(offset_pos);

        self.delta_pos = np.array(delta_pos, dtype = np.float) /100.0     # convert from cm to m

        print "self.delta_pos == ", self.delta_pos


        self.delta_time = 0
        self.delta_time_scale = time_step_scaling

        self.delta_pos2 = np.array(delta_pos, dtype=np.float)**2
        self.inv_delta_pos = np.floor(1.0 / np.array(self.delta_pos, dtype=np.float))

        self.addModel(self.base_model, placement_position=offset_pos)




    def addModel(self, new_model, onTopOfObject = [], offsetFromCenter = (0,0,0), placement_position = (0,0,0), place_next_to_object = -1, placement_mode = center, \
                 initial_condition = -9999, boundary_condition = 0, bc_type = 0, dirichlet_bc_value = 0,\
                 neumann_bc_value = 0):

        s = new_model.model.shape

        if self.main_model.ndim == 3:

            z = placement_position[0] + np.array([0, s[0]])
            y = placement_position[1] + np.array([0, s[1]])
            x = placement_position[2] + np.array([0, s[2]])

            self.main_model[z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.model
            # self.model_objects.append([new_model, np.array(placement_position, dtype=np.int)])
            # self.model_objects[new_model.model_id] = [new_model, np.array(placement_position, dtype = np.int)]
            self.model_objects[new_model.model_id] = MKModeOjectData(new_model, placement_position)

            # self.all_bc_coords.append(new_model.)

            self.D[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.D

            self.D_dict[new_model.model_id] = new_model.D

            if new_model.D >= self.max_D:
                self.max_D = new_model.D


            self.delta_time = (self.min_deltaPosition/100.0)**2/8/self.max_D/self.delta_time_scale

            self.modelDirichletBoundaryConditions[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.dirichlet_boundary_condition_data
            self.modelNeumannBoundaryConditions[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.neumann_boundary_condition_data

            if initial_condition == -9999:
                self.ui[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.ui
            else:
                self.main_model[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = initial_condition

            self.model_position.append(np.array(placement_position))
            self.boundaryCondition.addBoundaryCondition(0, new_model.model_id, boundary_condition, bc_type)

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
                self.ui[ y[0]:y[1], x[0]:x[1]] = initial_condition

            self.boundaryCondition.addBoundaryCondition(0, new_model.model_id, boundary_condition, bc_type)

            # self.D[new_model.model_id] = new_model.D

        else:
            print 'ERROR: Model must be 2 or 3 dimensional, given ', self.main_model.shape
            return



    def shrinkMainModel(self, pad = 4):

        idx = self.main_model[:] != 0

        pos = np.array(self.getMainModelObjectPositions(), dtype = np.int)

        xmin = np.min(pos[:,2])
        ymin = np.min(pos[:,1])
        zmin = np.min(pos[:,0])

        xsize = np.max(pos[:,2]) - xmin + 1
        ysize = np.max(pos[:,1]) - ymin + 1
        zsize = np.max(pos[:,0]) - zmin + 1

        offset = pad/2

        new_xsize = xsize + 4
        new_ysize = ysize + 4
        new_zsize = zsize + 4

        new_model = np.zeros((new_zsize, new_ysize, new_xsize), dtype= np.byte)

        new_model[offset:(offset+zsize), offset:(offset+ysize), offset:(offset+xsize)] = \
                                        self.main_model[zmin:(zmin+zsize), ymin:(ymin+ysize), xmin:(xmin+xsize)]

        self.main_model = sp.copy(new_model)

        new_model[:] = 0

        new_model[offset:(offset+zsize), offset:(offset+ysize), offset:(offset+xsize)] = \
                                        self.ui[zmin:(zmin+zsize), ymin:(ymin+ysize), xmin:(xmin+xsize)]

        self.ui = np.array(new_model.copy(), dtype = np.float)
        self.ui[self.main_model == 0] = np.nan

        self.u = sp.copy(self.ui)


        new_model[offset:(offset+zsize), offset:(offset+ysize), offset:(offset+xsize)] = \
                                        self.D[zmin:(zmin+zsize), ymin:(ymin+ysize), xmin:(xmin+xsize)]

        self.D = new_model.copy()

        # for









    def finalizeSystemSettings(self):

        if self.base_model.model_ndim == 3:
            self.shrinkMainModel()

            self.mainModelObjectCoords = np.array(self.getMainModelObjectPositions(), dtype = np.int)

            self.airContactSurface = {}
            self.materialContactSurface = {}

            self.calculateContactSurfaces()

            for mat_id in self.materialContactSurface.keys():
                coords = np.array(self.materialContactSurface[mat_id].coords, dtype=np.int)
                self.u_neumann[coords[:,0], coords[:,1], coords[:,2]] = self.materialContactSurface[mat_id].bcValue


            ds = self.delta_pos2[0]*self.delta_pos2[1] + self.delta_pos2[1]*self.delta_pos2[2] + self.delta_pos2[2]*self.delta_pos2[0]

            # print "ds == ", ds
            # print "self.D  == ", self.D
            # print "self.delta_time_scale == ", self.delta_time_scale

            # last_err_setting = np.seterr(divide= 'ignore')
            #
            # self.delta_time = (self.delta_pos2[0] * self.delta_pos2[1] * self.delta_pos2[2]) / (2*2*2 * self.D * ds) / self.delta_time_scale
            #
            # self.delta_time[np.isinf(self.delta_time) | np.isnan(self.delta_time)] = 0
            #
            #
            # np.seterr(**last_err_setting)


            self.sxx = self.delta_time/self.delta_pos2[2]
            self.syy = self.delta_time/self.delta_pos2[1]
            self.szz = self.delta_time/self.delta_pos2[0]

        elif self.base_model.model_ndim == 2:
            last_err_setting = np.seterr(divide= 'ignore')

            self.delta_time = np.prod(self.delta_pos2)/(2*2 * self.D * np.sum(self.delta_pos2)) / self.delta_time_scale
            self.delta_time[np.isinf(self.delta_time) | np.isnan(self.delta_time)] = 0

            np.seterr(**last_err_setting)

        else:
            raise Exception("ERROR: Must be 2 or 3 dimensional model")

        self.delta_time = np.min(self.delta_time[self.delta_time > 0])

        print "delta time is %.4e %s" % (self.delta_time,'  seconds')

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

        print 'Calculating object coordinates....'

        pos = np.array([])
        for z in xrange(0, zs):
            for y in xrange(0, ys):
                for x in xrange(0, xs):
                    if self.main_model[z,y,x] != 0:
                        pos = np.append(pos, [z,y,x])

        pos = np.reshape(pos,(-1,3))

        t1 = tm.time()
        print 'Time taken %f seconds' % (t1-t0,)
        print 'getMainModelObjectPositions .... pos shape = ', pos.shape
        print 'Object coordinate calculation done!'
        return pos


    def calculateContactSurfaces(self):

        print 'Calculating contact surface...'

        for pos in self.mainModelObjectCoords:
            self.updateContactNN(pos)

        for k in self.airContactSurface.keys():
            self.airContactSurface[k].coords = np.reshape(self.airContactSurface[k].coords, (-1,3))

        for k in self.materialContactSurface.keys():
            self.materialContactSurface[k].coords = np.reshape(self.materialContactSurface[k].coords, (-1,3))


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

        coords = self.mainModelObjectCoords
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

        plane = plane.reshape((-1,3))

        mdata = MKBoundaryConditionData()
        mdata.coords = np.array(plane, dtype = np.int)
        mdata.mat_id = mat_id
        mdata.bcType = 0
        mdata.bcValue = bc_value
        mdata.fluxDirection = direction

        if self.neumanBCCoordinates is None:
            self.neumanBCCoordinates = [mdata]
        else:
            self.neumanBCCoordinates.append(mdata)

        for m_dt in self.neumanBCCoordinates:
            m_dt.coords = m_dt.coords.reshape((-1,3))
            print 'self.neumannBCCoordinates mat %d shape %s' % (m_dt.mat_id, m_dt.coords.shape, )


    def setAllOuterSurfaceNeumannBC(self):
        nc = ['zbottom', 'ztop' , 'ybottom', 'ytop', 'xbottom', 'xtop']
        pass


    def cleanNeumannBCCoordinates(self):
        pass


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

            # print 'z idx shape', idx.shape
            # raw_input('Enter to continue')

            dmax = np.max(coords[:, 0])
            dmin = np.min(coords[:, 0])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.ZMAX,\
                                      coords = coords[coords[:, 0] == dmax, :])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.ZMIN,\
                                      coords = coords[coords[:, 0] == dmin, :])


            dmax = np.max(coords[:, 1])
            dmin = np.min(coords[:, 1])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.YMAX,\
                                      coords = coords[coords[:, 1] == dmax, :])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.YMIN,\
                                      coords = coords[coords[:, 1] == dmin, :])


            dmax = np.max(coords[:, 2])
            dmin = np.min(coords[:, 2])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.XMAX,\
                                      coords = coords[coords[:, 2] == dmax, :])

            self.updateNeumannBCArray(mat_id = m_id, bcType = 1, bcValue = bc_value, fluxDirection= MKDirection.XMIN,\
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
    @jit(nogil = True, cache = True)
    def sum_2nd_derivative_z(ui):
        # ui = np.array(ui, dtype=np.float64)
        new_u = ui[2:, 1:-1, 1:-1] + ui[:-2, 1:-1, 1:-1] - 2 * ui[1:-1, 1:-1, 1:-1]
        new_u[np.isnan(new_u)] = 0
        return new_u

    # @vectorize([float64(float64)])
    @jit(nogil = True, cache = True)
    def sum_2nd_derivative_y(ui):
        # ui = np.array(ui, dtype=np.float64)
        new_u = ui[1:-1, 2:, 1:-1] + ui[1:-1, :-2, 1:-1] - 2 * ui[1:-1, 1:-1, 1:-1]
        new_u[np.isnan(new_u)] = 0
        return new_u

    # @vectorize([float64(float64)])
    @jit(nogil = True, cache = True)
    def sum_2nd_derivative_x(ui):
        # ui = np.array(ui, dtype=np.float64)
        new_u = ui[1:-1, 1:-1, 2:] + ui[1:-1, 1:-1, :-2] - 2 * ui[1:-1, 1:-1, 1:-1]
        new_u[np.isnan(new_u)] = 0
        return new_u



    def evolveTempWithTime_v6(self, enableDiri = True, enableNeumann = True, data = None):


        @vectorize([float64(float64, float64, float64, float64, float64, float64)])
        def update_temperature(u, ui, D, sxx, syy, szz):
            u[1:-1,1:-1,1:-1] = ui[1:-1,1:-1,1:-1] + D[1:-1,1:-1,1:-1] * \
                        ( np.nan_to_num(ui[2:, 1:-1, 1:-1] + ui[:-2, 1:-1, 1:-1] - 2 * ui[1:-1, 1:-1, 1:-1] ) * szz + \
                          np.nan_to_num(ui[1:-1, 2:, 1:-1] + ui[1:-1, :-2, 1:-1] - 2 * ui[1:-1, 1:-1, 1:-1] ) * syy + \
                          np.nan_to_num(ui[1:-1, 1:-1, 2:] + ui[1:-1, 1:-1, :-2] - 2 * ui[1:-1, 1:-1, 1:-1] ) * sxx   \
                        )
            return u

        # @vectorize(float64, float64, float64, float64, float64, float64)
        # def update_1st_derivative_temperature(u, ui, D, delta_pos, svv):
        #     pass


        # update_temperature(self.u, self.ui, self.D, np.array(self.sxx, dtype=np.float64), np.array(self.syy, dtype=np.float), np.array(self.szz, dtype=np.float))
        update_temperature(self.u, self.ui, self.D, self.sxx, self.syy, self.szz)

        # self.u[1:-1,1:-1,1:-1] = self.ui[1:-1,1:-1,1:-1] + self.D[1:-1,1:-1,1:-1] * \
        #             ( np.nan_to_num(self.ui[2:, 1:-1, 1:-1] + self.ui[:-2, 1:-1, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1] ) * self.szz + \
        #               np.nan_to_num(self.ui[1:-1, 2:, 1:-1] + self.ui[1:-1, :-2, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1] ) * self.syy + \
        #               np.nan_to_num(self.ui[1:-1, 1:-1, 2:] + self.ui[1:-1, 1:-1, :-2] - 2 * self.ui[1:-1, 1:-1, 1:-1] ) * self.sxx   \
        #             )


        # Dirichlet bundary condition with air contact surface

        if enableDiri == True:
            for d_data in self.dirichletBCCoordinates:
                pos = d_data.coords
                self.u[pos[:,0], pos[:,1], pos[:, 2]] = d_data.bcValue

        # Neumann boundary condition along z-axis

        if enableNeumann == True:
            for ndata in self.neumanBCCoordinates:
                pos = ndata.coords
                mat_id = ndata.mat_id
                bc_val = ndata.bcValue
                dfactor = self.D_dict[mat_id] * self.delta_time

                # print 'Neumann bc value ', ndata.bcValue
                # print 'mat_id %d    bcValue %d  flux   ' % (ndata.mat_id, ndata.bcValue, ndata.fluxDirection,)

                if ndata.fluxDirection == MKDirection.ZMAX:
                    # print 'ZMAX calculating'
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.szz * (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2

                    # delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    # perp_axis = [0,1,1]
                elif ndata.fluxDirection == MKDirection.ZMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.szz * (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    # delta_upos = (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    # perp_axis = [0,1,1]
                elif ndata.fluxDirection == MKDirection.YMAX:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.syy * (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # perp_axis = [1,0,1]

                elif ndata.fluxDirection == MKDirection.YMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.syy * (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # delta_upos = (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # delta_upos[np.isnan(delta_upos)] = 0
                    # perp_axis = [1,0,1]

                elif ndata.fluxDirection == MKDirection.XMAX:

                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.sxx * (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1], pos[:,2]-1] + self.delta_pos[2] * bc_val) ) *2
                    # delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] + self.delta_pos[2] * bc_val) ) *2
                    # perp_axis = [1,1,0]

                elif ndata.fluxDirection == MKDirection.XMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.sxx * (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[2] * bc_val) ) *2
                    # delta_upos = (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[2] * bc_val) ) *2
                    # perp_axis = [1,1,0]

        self.ui[:] = self.u[:] #sp.copy(self.u)


    def evolveTempWithTime_v5p1(self, enableDiri = True, enableNeumann = True, data = None):

        self.u[:] = self.ui[:]

        d1 = MKMultiMaterialSystem.sum_2nd_derivative_z(self.ui)
        d2 = MKMultiMaterialSystem.sum_2nd_derivative_y(self.ui)
        d3 = MKMultiMaterialSystem.sum_2nd_derivative_x(self.ui)

        self.u[1:-1,1:-1,1:-1] = self.ui[1:-1,1:-1,1:-1] + self.D[1:-1,1:-1,1:-1] * \
                    ( d1 * self.delta_time/self.delta_pos2[0] + \
                      d2 * self.delta_time/self.delta_pos2[1] + \
                      d3 * self.delta_time/self.delta_pos2[2]   \
                    )


        # Dirichlet bundary condition with air contact surface

        if enableDiri == True:
            for d_data in self.dirichletBCCoordinates:
                pos = d_data.coords
                self.u[pos[:,0], pos[:,1], pos[:, 2]] = d_data.bcValue

        # Neumann boundary condition along z-axis

        if enableNeumann == True:
            for ndata in self.neumanBCCoordinates:
                pos = ndata.coords
                mat_id = ndata.mat_id
                bc_val = ndata.bcValue
                dfactor = self.D_dict[mat_id] * self.delta_time

                # print 'Neumann bc value ', ndata.bcValue
                # print 'mat_id %d    bcValue %d  flux   ' % (ndata.mat_id, ndata.bcValue, ndata.fluxDirection,)

                if ndata.fluxDirection == MKDirection.ZMAX:
                    # print 'ZMAX calculating'
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.szz * (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2

                    # delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    # perp_axis = [0,1,1]
                elif ndata.fluxDirection == MKDirection.ZMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.szz * (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    # delta_upos = (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    # perp_axis = [0,1,1]
                elif ndata.fluxDirection == MKDirection.YMAX:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.syy * (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # perp_axis = [1,0,1]

                elif ndata.fluxDirection == MKDirection.YMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.syy * (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # delta_upos = (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # delta_upos[np.isnan(delta_upos)] = 0
                    # perp_axis = [1,0,1]

                elif ndata.fluxDirection == MKDirection.XMAX:

                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.sxx * (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1], pos[:,2]-1] + self.delta_pos[2] * bc_val) ) *2
                    # delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] + self.delta_pos[2] * bc_val) ) *2
                    # perp_axis = [1,1,0]

                elif ndata.fluxDirection == MKDirection.XMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.sxx * (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[2] * bc_val) ) *2
                    # delta_upos = (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[2] * bc_val) ) *2
                    # perp_axis = [1,1,0]

        self.ui[:] = self.u[:] #sp.copy(self.u)

    def evolveTempWithTime_v5(self, enableDiri = True, enableNeumann = True, data = None):

        self.u[:] = self.ui[:]

        self.u[1:-1,1:-1,1:-1] = self.ui[1:-1,1:-1,1:-1] + self.D[1:-1,1:-1,1:-1] * \
                    ( np.nan_to_num(self.ui[2:, 1:-1, 1:-1] + self.ui[:-2, 1:-1, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1] ) * self.szz + \
                      np.nan_to_num(self.ui[1:-1, 2:, 1:-1] + self.ui[1:-1, :-2, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1] ) * self.syy + \
                      np.nan_to_num(self.ui[1:-1, 1:-1, 2:] + self.ui[1:-1, 1:-1, :-2] - 2 * self.ui[1:-1, 1:-1, 1:-1] ) * self.sxx   \
                    )


        # Dirichlet bundary condition with air contact surface

        if enableDiri == True:
            for d_data in self.dirichletBCCoordinates:
                pos = d_data.coords
                self.u[pos[:,0], pos[:,1], pos[:, 2]] = d_data.bcValue

        # Neumann boundary condition along z-axis

        if enableNeumann == True:
            for ndata in self.neumanBCCoordinates:
                pos = ndata.coords
                mat_id = ndata.mat_id
                bc_val = ndata.bcValue
                dfactor = self.D_dict[mat_id] * self.delta_time

                # print 'Neumann bc value ', ndata.bcValue
                # print 'mat_id %d    bcValue %d  flux   ' % (ndata.mat_id, ndata.bcValue, ndata.fluxDirection,)

                if ndata.fluxDirection == MKDirection.ZMAX:
                    # print 'ZMAX calculating'
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.szz * (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2

                    # delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    # perp_axis = [0,1,1]
                elif ndata.fluxDirection == MKDirection.ZMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.szz * (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    # delta_upos = (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    # perp_axis = [0,1,1]
                elif ndata.fluxDirection == MKDirection.YMAX:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.syy * (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # perp_axis = [1,0,1]

                elif ndata.fluxDirection == MKDirection.YMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.syy * (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # delta_upos = (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    # delta_upos[np.isnan(delta_upos)] = 0
                    # perp_axis = [1,0,1]

                elif ndata.fluxDirection == MKDirection.XMAX:

                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.sxx * (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1], pos[:,2]-1] + self.delta_pos[2] * bc_val) ) *2
                    # delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] + self.delta_pos[2] * bc_val) ) *2
                    # perp_axis = [1,1,0]

                elif ndata.fluxDirection == MKDirection.XMIN:
                    self.u[pos[:,0], pos[:,1], pos[:,2]] += self.sxx * (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[2] * bc_val) ) *2
                    # delta_upos = (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[2] * bc_val) ) *2
                    # perp_axis = [1,1,0]

        self.ui[:] = self.u[:] #sp.copy(self.u)



    def evolveTempWithTime_v4(self, enableDiri = True, enableNeumann = True):

        # self.u[1:-1, 1:-1, 1:-1] = self.ui[1:-1, 1:-1, 1:-1]
        self.u[:] = self.ui[:]

        tmp_data = sp.copy(self.ui)

        for mid, d_val in self.D_dict.iteritems():

            p = self.model_objects[mid].position
            zs,ys,xs = self.model_objects[mid].model.model.shape

            pz = p[0]+1
            py = p[1]+1
            px = p[2]+1

            pzm = pz + (zs-1) - 1
            pym = py + (ys-1) - 1
            pxm = px + (xs-1) - 1

            self.u[pz:pzm, py:pym, px:pxm] += d_val * self.delta_time * \
                    ( \
                      (self.ui[(pz+1):(pzm+1), py:pym, px:pxm] + self.ui[(pz-1):(pzm-1), py:pym, px:pxm] - 2 * self.ui[pz:pzm, py:pym, px:pxm] ) / self.delta_pos2[0] + \
                      (self.ui[pz:pzm, (py+1):(pym+1), px:pxm] + self.ui[pz:pzm, (py-1):(pym-1), px:pxm] - 2 * self.ui[pz:pzm, py:pym, px:pxm] ) / self.delta_pos2[1] + \
                      (self.ui[pz:pzm, py:pym, (px+1):(pxm+1)] + self.ui[pz:pzm, py:pym, (px-1):(pxm-1)] - 2 * self.ui[pz:pzm, py:pym, px:pxm] ) / self.delta_pos2[2]  \
                    )


        newT = None

        for mat_id, data in self.materialContactSurface.iteritems():
            pos = data.coords
            # self.ui[pos[:,0], pos[:,1], pos[:,2]] +

            newT = self.D_dict[mat_id] * self.delta_time * \
                   (\
                       (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] +  self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]) / self.delta_pos2[0] + \
                       (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] +  self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]) / self.delta_pos2[1] + \
                       (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] +  self.ui[pos[:,0], pos[:,1], pos[:,2]-1] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]) / self.delta_pos2[2] \
                   )

            newT[np.isnan(newT)] = 0
            self.u[pos[:,0], pos[:,1], pos[:,2]] += newT
            newT = []


            # zz,yy,xx = np.meshgrid(np.linspace(1, 1+zs))
            #
            # for z in xrange(pz, pzm):
            #     for y in xrange(py, pym):
            #         for x in xrange(px, pxm):
            #             self.u[z,y,x] = self.ui[z,y,x] + self.delta_time * d_val * \
            #                               (\
            #                                   (self.ui[z+1,y,x] + self.ui[z-1,y,x] - 2 * self.ui[z,y,x])/self.delta_pos2[0] + \
            #                                   (self.ui[z,y+1,x] + self.ui[z,y-1,x] - 2 * self.ui[z,y,x])/self.delta_pos2[1] + \
            #                                   (self.ui[z,y,x+1] + self.ui[z,y,x-1] - 2 * self.ui[z,y,x])/self.delta_pos2[2]  \
            #                               )


        # for mid, d_val in self.D_dict.iteritems():
        #     tmp_data[self.main_model[:] != mid] = np.nan


        # for mat_id, data in self.D_dict:

        # Dirichlet bundary condition with air contact surface

        if enableDiri == True:
            for d_data in self.dirichletBCCoordinates:
                pos = d_data.coords
                self.u[pos[:,0], pos[:,1], pos[:, 2]] = d_data.bcValue

        # Neumann boundary condition along z-axis

        if enableNeumann == True:
            for ndata in self.neumanBCCoordinates:
                pos = ndata.coords
                self.u_neumann[:] = 0
                self.u_neumann[pos[:,0], pos[:,1], pos[:,2]] = ndata.bcValue
                mat_id = ndata.mat_id

                bc_val = ndata.bcValue

                # print 'Neumann bc value ', ndata.bcValue
                # print 'mat_id %d    bcValue %d  flux   ' % (ndata.mat_id, ndata.bcValue, ndata.fluxDirection,)

                if ndata.fluxDirection == MKDirection.ZMAX:
                    # print 'ZMAX calculating'
                    delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    perp_axis = [0,1,1]
                elif ndata.fluxDirection == MKDirection.ZMIN:
                    delta_upos = (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[0] * bc_val) ) *2
                    perp_axis = [0,1,1]
                elif ndata.fluxDirection == MKDirection.YMAX:
                    delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    perp_axis = [1,0,1]

                elif ndata.fluxDirection == MKDirection.YMIN:
                    delta_upos = (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[1] * bc_val) ) *2
                    delta_upos[np.isnan(delta_upos)] = 0
                    perp_axis = [1,0,1]

                elif ndata.fluxDirection == MKDirection.XMAX:
                    delta_upos = (-self.ui[pos[:,0], pos[:,1], pos[:,2]] + (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] + self.delta_pos[2] * bc_val) ) *2
                    perp_axis = [1,1,0]

                elif ndata.fluxDirection == MKDirection.XMIN:
                    delta_upos = (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[2] * bc_val) ) *2
                    perp_axis = [1,1,0]

                delta_upos[np.isnan(delta_upos)] = 0

                dfactor = self.D_dict[mat_id] * self.delta_time

                if perp_axis[0] == 0:
                    newT = dfactor*(delta_upos + \
                       (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] +  self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]) / self.delta_pos2[1] + \
                       (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] +  self.ui[pos[:,0], pos[:,1], pos[:,2]-1] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]) / self.delta_pos2[2] )

                elif perp_axis[1] == 0:
                    newT = dfactor * (delta_upos + \
                       (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] +  self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]) / self.delta_pos2[1] + \
                       (self.ui[pos[:,0], pos[:,1], pos[:,2]+1] +  self.ui[pos[:,0], pos[:,1], pos[:,2]-1] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]) / self.delta_pos2[2] )
                elif perp_axis[2] == 0:
                    newT = dfactor * (delta_upos + \
                       (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] +  self.ui[pos[:,0]-1, pos[:,1], pos[:,2]] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]) / self.delta_pos2[1] + \
                       (self.ui[pos[:,0], pos[:,1]+1, pos[:,2]] +  self.ui[pos[:,0], pos[:,1]-1, pos[:,2]] - 2 * self.ui[pos[:,0], pos[:,1], pos[:,2]]) / self.delta_pos2[2])

                newT[np.isnan(newT)] = 0


                self.u[pos[:,0], pos[:,1], pos[:,2]] += newT
                newT = None

        self.ui = sp.copy(self.u)




    def evolveTempWithTime_v3(self, enableDiri = True, enableNeumann = True):

        # self.u[1:-1, 1:-1, 1:-1] = self.ui[1:-1, 1:-1, 1:-1]
        self.u[:] = self.ui[:]

        tmp_data = sp.copy(self.ui)

        for mid, d_val in self.D_dict.iteritems():
            self.u[1:-1,1:-1,1:-1] += self.D[1:-1, 1:-1, 1:-1] * self.delta_time * \
                    ( (np.nan_to_num(self.ui[2:, 1:-1, 1:-1]+ self.ui[:-2, 1:-1, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1])) / self.delta_pos2[0] + \
                      (np.nan_to_num(self.ui[1:-1, 2:, 1:-1]+ self.ui[1:-1, :-2, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1])) / self.delta_pos2[1] + \
                      (np.nan_to_num(self.ui[1:-1, 1:-1, 2:]+ self.ui[1:-1, 1:-1, :-2] - 2 * self.ui[1:-1, 1:-1, 1:-1])) / self.delta_pos2[2]  \
                    )

        # for mid, d_val in self.D_dict.iteritems():
        #     tmp_data[self.main_model[:] != mid] = np.nan






        # Dirichlet bundary condition with air contact surface

        if enableDiri:
            for d_data in self.dirichletBCCoordinates:
                coords = d_data.coords
                self.u[coords[:,0], coords[:,1], coords[:, 2]] = d_data.bcValue

        # Neumann boundary condition along z-axis

        if enableNeumann:
            for ndata in self.neumanBCCoordinates:
                coords = ndata.coords
                self.u_neumann[:] = 0
                self.u_neumann[coords[:,0], coords[:,1], coords[:,2]] = ndata.bcValue

                print 'Neumann bc value ', ndata.bcValue

                # print 'mat_id %d    bcValue %d  flux   ' % (ndata.mat_id, ndata.bcValue, ndata.fluxDirection,)

                if ndata.fluxDirection == MKDirection.ZMAX:
                    ds2 = self.delta_pos2[0]
                    ds = self.delta_pos[0]
                    u_delta = self.ui[2:, 1:-1, 1:-1]
                elif ndata.fluxDirection == MKDirection.ZMIN:
                    ds2 = self.delta_pos2[0]
                    ds = self.delta_pos[0]
                    u_delta = self.ui[:-2, 1:-1, 1:-1]
                elif ndata.fluxDirection == MKDirection.YMAX:
                    ds2 = self.delta_pos2[1]
                    ds = self.delta_pos[1]
                    u_delta = self.ui[1:-1, 2:, 1:-1]
                elif ndata.fluxDirection == MKDirection.YMIN:
                    ds2 = self.delta_pos2[1]
                    ds = self.delta_pos[1]
                    u_delta = self.ui[1:-1, :-2, 1:-1]
                elif ndata.fluxDirection == MKDirection.XMAX:
                    ds2 = self.delta_pos2[2]
                    ds = self.delta_pos[2]
                    u_delta = self.ui[1:-1, 1:-1, 2:]
                elif ndata.fluxDirection == MKDirection.XMIN:
                    ds2 = self.delta_pos2[2]
                    ds = self.delta_pos[2]
                    u_delta = self.ui[1:-1, 1:-1, :-2]

                self.u[1:-1,1:-1,1:-1] +=  (self.u_neumann[1:-1,1:-1,1:-1] != 0) * self.delta_time * self.D_dict[ndata.mat_id]/ds2 *\
                                      2 * ( u_delta - (self.ui[1:-1, 1:-1, 1:-1] +  ds * self.u_neumann[1:-1, 1:-1, 1:-1]) )



        self.ui = sp.copy(self.u)


    def evolveTempWithTime_v2(self):

        # for k in self.airContactSurface.keys():
        #     coords = np.array(self.airContactSurface[k].coords, dtype=np.int)
        #     self.ui[coords[:,0], coords[:,1], coords[:, 2]] = 20
            # self.u = self.ui
        #
        # self.u[2:, 1:-1, 1:-1] = self.u[2:, 1:-1,1:-1] + self.u_neumann

        self.u[1:-1, 1:-1, 1:-1] = self.ui[1:-1, 1:-1, 1:-1] + self.D[1:-1, 1:-1, 1:-1] * self.delta_time * \
                ( (self.ui[2:, 1:-1, 1:-1] + self.ui[:-2, 1:-1, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[0] + \
                  (self.ui[1:-1, 2:, 1:-1] + self.ui[1:-1, :-2, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[1] + \
                  (self.ui[1:-1, 1:-1, 2:] + self.ui[1:-1, 1:-1, :-2] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[2]  \
                + \
                  0*(self.u_neumann[1:-1, 1:-1, 1:-1] != 0) * (2 * (self.ui[2:, 1:-1, 1:-1]) - (self.ui[1:-1, 1:-1, 1:-1] +  self.delta_pos[0] * 0*self.u_neumann[1:-1, 1:-1, 1:-1]))/ self.delta_pos2[0] \
                )


        # # Dirichlet bundary condition with air contact surface
        # for k in self.airContactSurface.keys():
        #     coords = np.array(self.airContactSurface[k].coords, dtype=np.int)
        #     self.ui[coords[:,0], coords[:,1], coords[:, 2]] = 20
            # self.u = self.ui

        # Dirichlet bundary condition with air contact surface

        for bdata in self.dirichletBCCoordinates:
            coords = bdata.coords
            self.u[coords[:,0], coords[:,1], coords[:, 2]] = bdata.bcValue


        # Neumann boundary condition along z-axis


        for ndata in self.neumanBCCoordinates:
            coords = ndata.coords
            self.u_neumann[:] = 0
            self.u_neumann[coords[:,0], coords[:,1], coords[:,2]] = ndata.bcValue

            print 'mat_id %d    bcValue %d  flux   '

            if ndata.fluxDirection == MKDirection.ZMAX:
                ds2 = self.delta_pos2[0]
                ds = self.delta_pos[0]
                u_delta = self.ui[2:, 1:-1, 1:-1]
            elif ndata.fluxDirection == MKDirection.ZMIN:
                ds2 = self.delta_pos2[0]
                ds = self.delta_pos[0]
                u_delta = self.ui[:-2, 1:-1, 1:-1]
            elif ndata.fluxDirection == MKDirection.YMAX:
                ds2 = self.delta_pos2[1]
                ds = self.delta_pos[1]
                u_delta = self.ui[1:-1, 2:, 1:-1]
            elif ndata.fluxDirection == MKDirection.YMIN:
                ds2 = self.delta_pos2[1]
                ds = self.delta_pos[1]
                u_delta = self.ui[1:-1, 2:, 1:-1]
            elif ndata.fluxDirection == MKDirection.XMAX:
                ds2 = self.delta_pos2[2]
                ds = self.delta_pos[2]
                u_delta = self.ui[1:-1, 1:-1, 2:]
            elif ndata.fluxDirection == MKDirection.XMIN:
                ds2 = self.delta_pos2[2]
                ds = self.delta_pos[2]
                u_delta = self.ui[1:-1, 1:-1, :-2]

            self.u[1:-1,1:-1,1:-1] +=  (self.u_neumann != 0) * self.delta_time * self.D_dict[ndata.mat_id]/ds2 *\
                                  2 * ( u_delta - (self.ui[1:-1, 1:-1, 1:-1] +  ds * self.u_neumann[1:-1, 1:-1, 1:-1]) )



        # self.u[1:-1, 1:-1, 1:-1] = self.u[1:-1, 1:-1, 1:-1] + (self.u_neumann[1:-1, 1:-1, 1:-1] != 0) * self.D[1:-1, 1:-1, 1:-1] * self.delta_time * \
        #                         (2 * (self.ui[2:, 1:-1, 1:-1]) - (self.ui[1:-1, 1:-1, 1:-1] -  self.delta_pos[0] * self.u_neumann[1:-1, 1:-1, 1:-1]))/ self.delta_pos2[0]


        # self.u[1:-1, 1:-1, 1:-1] = self.u[1:-1, 1:-1, 1:-1] + (2 * (self.u[2:, 1:-1, 1:-1]  - self.u[1:-1, 1:-1, 1:-1] - 2*

        # for mat_id in self.materialContactSurface.keys():
        #     pos = np.array(self.materialContactSurface[mat_id].coords, dtype=np.int)
        #     neumann_cond = self.materialContactSurface[mat_id].bcValue
        #
        #     print 'pos shape ', pos.shape
        #
        #     # for z axis only... for the moment
        #     # self.u[pos[:,0], pos[:,1], pos[:2]] =\
        #     #     self.u[pos[:,0], pos[:,1], pos[:,2]] + self.Dself.delta_time * \
        #     #     (2 * (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[0]* neumann_cond ) )/ self.delta_pos[0]\
        #     #      + (self.ui[pos[:,0], pos[:,1]+1, pos[:2]] + self.ui[pos[:,0], pos[:,1]-1, pos[:2]] - 2 * self.ui[pos[:,0], pos[:,1], pos[:2]]) / self.delta_pos2[1] \
        #     #      + (self.ui[pos[:,0], pos[:,1], pos[:2]+1] + self.ui[pos[:,0], pos[:,1], pos[:2]-1] - 2 * self.ui[pos[:,0], pos[:,1], pos[:2]]) / self.delta_pos2[2] \
        #     #      ) # plus a heat source here
        #
        #     self.u[pos[:,0], pos[:,1], pos[:2]] =\
        #         self.u[pos[:,0], pos[:,1], pos[:,2]] + self.D[pos[:,0], pos[:,1], pos[:,2]] * self.delta_time * \
        #         (2 * (self.ui[pos[:,0]+1, pos[:,1], pos[:,2]] - (self.ui[pos[:,0], pos[:,1], pos[:,2]] + self.delta_pos[0]* neumann_cond ) )/ self.delta_pos[0]\
        #          ) # plus a heat source here
        #
        # T(1,j+1)=T(1,j)+lambda*(T(2,j)-2*T(1,j)+T(2,j));
        # T(m+1,j+1)=T(m+1,j)+lambda*(T(m,j)-2*T(m+1,j)+(T(m,j)+2*dx*2));

        self.ui = sp.copy(self.u)

    def startSimulation(self, maxTimeIteration = 1000):

        for i in xrange(maxTimeIteration):
            self.evolveT()


    def plot2DSectionData(self, fig_no, data2D):

        f = plt.figure(fig_id)
        plt.clf()
        plt.cla()
        im = plt.imshow(data2D)
        plt.colorbar(im)


    def plotSurface3D(self):
        pass


    def runSimulation(self, iterationTime = 1000, initU = None, figureID = 666, ySection = 15, applyDiri = True, applyNeumann = True):
        ######################
        self.ui[:] = initU[:]

        data = []

        time_steps = iterationTime
        count = 1
        tic()

        temp_dict = {}
        for mid, d_val in self.D_dict.iteritems():
            temp_dict[mid] = np.array([])

        time_data = np.array([])

        for i in xrange(time_steps):
            # self.evolveTempWithTime_v5p1(enableDiri=applyDiri, enableNeumann=applyNeumann)
            self.evolveTempWithTime_v3(enableDiri=applyDiri, enableNeumann=applyNeumann)

            for mid, d_val in self.D_dict.iteritems():
                p = self.model_objects[mid].position
                zs,ys,xs = self.model_objects[mid].model.model.shape

                pz = p[0]
                py = p[1]
                px = p[2]

                pzm = pz + zs
                pym = py + ys
                pxm = px + xs

                temp_dict[mid] = np.append(temp_dict[mid], np.mean(self.u[pz:pzm, py:pym, px:pxm]))

            time_data = np.append(time_data, i * self.delta_time)

            progress =  np.ceil(float(i)/time_steps*100)

            if progress >= count*10:
                count += 1
                print 'Simulation completed %d' % (progress,),
                print '%'

        plot2DImage(figureID, np.nan_to_num(self.u[:,ySection,:]), xAxisLabel= 'X (mm)', yAxisLabel= 'Y (mm)')

        print 'Total simulation time = %.4e seconds with time steps of %f ms' % (time_steps * self.delta_time, self.delta_time*1000.0,)
        t2 =  toc()

        print 'Total time taken ',

        if t2 < 60:
            print toc(), 'seconds'
        else:
            print toc()/60.0, 'minutes'

        return [time_data, temp_dict]



def plot2DImage(fig_id, image0, doInvertY = True, xAxisLabel = (), yAxisLabel=()):
    ff = plt.figure(fig_id)
    plt.clf()
    plt.cla()
    ax = ff.add_subplot(111)
    im = plt.imshow(image0)
    if doInvertY:
       ax.invert_yaxis()

    plt.colorbar(im)
    ax.set_xlabel(xAxisLabel)
    ax.set_ylabel(yAxisLabel)




tm_val = 0

def tic():
    global tm_val

    tm_val = tm.time()

def toc():
    global tm_val
    return tm.time() - tm_val