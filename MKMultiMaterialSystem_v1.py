
__author__ = 'kamal'


import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from matplotlib.tri import triangulation
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
import time as tm


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


    # def getSurfaceCoordinates(self):
    #
    #     if self.model_ndim == 2:
    #         s = self.model_size;
    #
    #         x = np.linspace(0, s[1]-1, s[1], dtype = np.int)
    #         y = np.linspace(0, s[0]-1, s[0], dtype = np.int)
    #
    #         [xx,yy] = np.meshgrid(x,y)
    #
    #         xx = np.reshape(xx, (-1,1))
    #         yy = np.reshape(yy, (-1,1))
    #
    #         coords = np.zeros(xx.__len__(),2, dtype= np.int)
    #         coords[:,0] = yy.reshape((yy.__len__()))
    #         coords[:,1] = xx.reshape((xx.__len__()))
    #
    #         surface_coords = []
    #
    #         for i in xrange(xx.__len__()):
    #             if self.isSurface(coords[i,:]) == True:
    #                 surface_coords.append(coords[i,:])
    #
    #         return np.array(surface_coords,dtype = np.int)
    #
    #
    #     elif self.model_ndim == 3:
    #         s = self.model_size;
    #
    #         print 'model size = ', s
    #
    #         x = np.linspace(0, s[2]-1, s[2], dtype = np.int)
    #         y = np.linspace(0, s[1]-1, s[1], dtype = np.int)
    #         z = np.linspace(0, s[0]-1, s[0], dtype = np.int)
    #
    #         [xx,yy,zz] = np.meshgrid(x,y,z)
    #
    #         xx = np.reshape(xx, (-1,1))
    #         yy = np.reshape(yy, (-1,1))
    #         zz = np.reshape(zz, (-1,1))
    #
    #         print "max zz = ", np.max(zz[:])
    #         print "max yy = ", np.max(yy[:])
    #         print "max xx = ", np.max(xx[:])
    #
    #         # zz_idx = zz[:] < 0
    #         # zz_idx[zz != 0] = False
    #         # zz_idx[s[2]-1] = False
    #         #
    #         # yy_idx = yy[:] < 0
    #         # yy_idx[0] = True
    #         # yy_idx[s[1]-1] = True
    #         #
    #         # xx_idx = xx[:] < 0
    #         # xx_idx[0] = True
    #         # xx_idx[s[0]-1] = True
    #
    #         zz_bottom = zz[:] == 0
    #         zz_top = (zz[:] == s[0]-1)
    #
    #         yy_down = yy[:] == 0
    #         yy_up = (yy[:] == s[1]-1)
    #
    #         xx_left = xx[:] == 0
    #         xx_right = (xx[:] == s[0]-1)
    #
    #         xx_idx = np.logical_or(xx_left, xx_right)
    #
    #         xplanes = np.array([zz[xx_left], yy[xx_left], xx[xx_left]]).reshape((-1,3))
    #         xplanes = np.append(xplanes, np.array([zz[xx_right], yy[xx_right], xx[xx_right]])).reshape((-1,3))
    #
    #         yplanes = np.array([zz[yy_up], yy[yy_up], xx[yy_up]]).reshape((-1,3))
    #         yplanes = np.append(yplanes, np.array([zz[yy_down], yy[yy_down], xx[yy_down]])).reshape((-1,3))
    #
    #         zplanes = np.array([zz[zz_bottom], yy[zz_bottom], xx[zz_bottom]]).reshape((-1,3))
    #         zplanes = np.append(zplanes, np.array([zz[zz_top], yy[zz_top], xx[zz_top]])).reshape((-1,3))
    #
    #         coords = np.append(zplanes, np.append(xplanes,yplanes)).reshape((-1,3))
    #
    #         print coords
    #
    #         # zz_idx = np.logical_or(zz == 0)
    #         # yy_idx = np.logical_not(np.logical_or(yy == 0, yy == s[1]-1))
    #         # xx_idx = np.logical_not(np.logical_or(xx == 0, xx == s[0]-1))
    #         #
    #         # idx = np.logical_and(np.logical_and(zz_idx, yy_idx), xx_idx)
    #         #
    #         # print "idx shape ", idx.shape
    #         #
    #         # coords = np.transpose(np.array([zz[idx], yy[idx], xx[idx]], dtype=np.int))
    #         #
    #
    #         print 'coord shape = ', coords.shape
    #
    #         # coords = np.zeros((xx.__len__(), 3), dtype = np.int)
    #         # coords[:,0] = zz.reshape((zz.__len__()))
    #         # coords[:,1] = yy.reshape((yy.__len__()))
    #         # coords[:,2] = xx.reshape((xx.__len__()))
    #
    #         surface_coords = coords.reshape((-1,3)) # coords[idx, :]
    #
    #         # surface_coords = []
    #         #
    #         #
    #         # surface_coords.append([0,coords[0,:,:]])
    #         # surface_coords.append([s[2]-1,coords[s[2]-1,:,:])
    #         #
    #         # surface_coords.append(coords[:,0,:])
    #         # surface_coords.append(coords[:,s[1]-1,:])
    #         #
    #         # surface_coords.append(coords[:,:,0])
    #         # surface_coords.append(coords[:,:,s[0]-1])
    #         #
    #         # for i in xrange(xx.__len__()):
    #         #     if self.isSurface(coords[i,:]) == True:
    #         #         surface_coords.append(coords[i,:])
    #         # for i in xrange(xx.__len__()):
    #         #     if self.isSurface(coords[i,:]) == True:
    #         #         surface_coords.append(coords[i,:])
    #
    #         return np.array(surface_coords, dtype = np.int)
    #     else:
    #         raise Exception("ERROR: Must be a 2 or 3 dimensional model")
    #
    #     return []


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

    def __init__(self, base_model,  max_materials = 10, delta_pos = (0.01,0.01, 0.01), time_step_scaling = 1, padding = (0.2, 0.2, 0.2)):
        self.base_model = base_model
        self.system_dimension = base_model.model_ndim

        self.main_model = [] #copy.deepcopy(base_model)
        self.model_objects = []
        self.model_object_base_coord = []
        self.boundaryCondition = MKBoundaryCondition()

        self.model_objects = []

        self.neumanBC = {}
        self.dirichletBC = {}

        # self.model_objects.append(base_model)
        self.model_position = []

        self.all_bc_coords = []
        self.all_dirichlet_value = []
        self.all_neumann_value = []

        new_size = np.round(self.base_model.model_size  * (1.0 + np.array(padding)))
        mx_val = np.max(new_size)
        new_size = np.array([mx_val, mx_val, mx_val])

        self.main_model = np.zeros(new_size, dtype=np.byte);
        self.ui = np.zeros(new_size, dtype=np.float)
        self.u = np.zeros(new_size, dtype=np.float)
        self.u_neumann = np.zeros(new_size, dtype=np.float)
        self.D = np.zeros(new_size, dtype=np.float)

        self.mainModelObjectCoords = [] ## save all the object coords

        self.airContactSurface = {}
        self.materialContactSurface = {}

        self.modelDirichletBoundaryConditions = np.zeros(new_size, dtype= float)
        self.modelNeumannBoundaryConditions = np.zeros(new_size, dtype= float)

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





    def addModel(self, new_model, placement_position, place_next_to_object = -1, placement_mode = center, \
                 initial_condition = -9999, boundary_condition = 0, bc_type = 0, dirichlet_bc_value = 0, neumann_bc_value = 0):

        s = new_model.model.shape

        if self.main_model.ndim == 2:
            x = placement_position[1] + np.array([0, s[1]])
            y = placement_position[0] + np.array([0, s[0]])

            self.main_model[ y[0]:y[1], x[0]:x[1]] = new_model.model
            self.model_objects.append([new_model, np.array(placement_position)])

            self.neumanBC[new_model.model_id] = neumann_bc_value
            self.dirichletBC[new_model.model_id] = dirichlet_bc_value

            self.D[ y[0]:y[1], x[0]:x[1]] = new_model.D

            print "self.D === ", self.D

#            self.modelDirichletBoundaryConditions[ y[0]:y[1], x[0]:x[1]] = new_model.

            if initial_condition == -9999:
                self.ui[ y[0]:y[1], x[0]:x[1]] = new_model.ui
            else:
                self.ui[ y[0]:y[1], x[0]:x[1]] = initial_condition

            self.boundaryCondition.addBoundaryCondition(0, new_model.model_id, boundary_condition, bc_type)

            # self.D[new_model.model_id] = new_model.D


        elif self.main_model.ndim == 3:

            x = placement_position[2] + np.array([0, s[2]])
            y = placement_position[1] + np.array([0, s[1]])
            z = placement_position[0] + np.array([0, s[0]])

            self.main_model[z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.model
            self.model_objects.append([new_model, np.array(placement_position)])

            # self.all_bc_coords.append(new_model.)

            self.D[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.D

            self.modelDirichletBoundaryConditions[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.dirichlet_boundary_condition_data
            self.modelNeumannBoundaryConditions[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.neumann_boundary_condition_data

            if initial_condition == -9999:
                self.ui[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = new_model.ui
            else:
                self.main_model[ z[0]:z[1], y[0]:y[1], x[0]:x[1]] = initial_condition

            self.model_position.append(np.array(placement_position))

            self.boundaryCondition.addBoundaryCondition(0, new_model.model_id, boundary_condition, bc_type)

        else:
            print 'ERROR: Model must be 2 or 3 dimensional, given ', self.main_model.shape
            return



    def finalizeSystemSettings(self):

        if self.base_model.model_ndim == 3:

            self.mainModelObjectCoords = self.getMainModelObjectPositions()

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

            last_err_setting = np.seterr(divide= 'ignore')

            self.delta_time = (self.delta_pos2[0] * self.delta_pos2[1] * self.delta_pos2[2]) / (2*2*2 * self.D * ds) / self.delta_time_scale
            self.delta_time[np.isinf(self.delta_time) | np.isnan(self.delta_time)] = 0
            np.seterr(**last_err_setting)


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


    # def calculateSurfaceContactsAndCoordinates(self):
    #
    #     zs,ys,xs = self.main_model.shape
    #
    #     xx,yy,zz = np.meshgrid(np.linspace(0, xs-1, xs), \
    #                            np.linspace(0, ys-1, ys), \
    #                            np.linspace(0, zs-1, zs)  \
    #                            )
    #
    #     xx = xx.reshape((-1,1))
    #     yy = yy.reshape((-1,1))
    #     zz = zz.reshape((-1,1))
    #
    #     pos = np.append(zz, np.append(yy, xx)).reshape((3,-1)).transpose()
    #
    #     print 'pos = ', pos
    #
    #     air_contact_coord = {}
    #     mat_contact_coord = {}
    #
    #
    #     for pa in pos:
    #         p = tuple(pa)
    #
    #         # print 'p = ', p
    #         air_cont, mat_cont =  self.isSurface(p)
    #
    #         # if self.isSurface(p, withAir=True) != 0:
    #         if air_cont is not None:
    #             cur_object_id = self.main_model[p]
    #
    #             if cur_object_id in air_contact_coord.keys():
    #                 air_contact_coord[cur_object_id].id = 0
    #                 air_contact_coord[cur_object_id].coords = np.append( air_contact_coord[cur_object_id].coords, p)
    #             else:
    #                 air_contact_coord[cur_object_id] = MKSurfaceContact()
    #                 air_contact_coord[cur_object_id].id = 0
    #                 air_contact_coord[cur_object_id].coords = np.array(p)
    #
    #         if mat_cont is not None:
    #             mat_id = mat_cont #self.isSurface(p, withAir=False)
    #
    #             if cur_object_id in mat_contact_coord.keys():
    #                 air_contact_coord[cur_object_id].id = mat_id
    #                 air_contact_coord[cur_object_id].coords = np.append( mat_contact_coord[cur_object_id].coords, p)
    #             else:
    #                 if mat_id in mat_contact_coord.keys():
    #                     if mat_contact_coord[mat_id].id != cur_object_id:
    #                         mat_contact_coord[cur_object_id] = MKSurfaceContact()
    #                         mat_contact_coord[cur_object_id].id = mat_id
    #                         mat_contact_coord[cur_object_id].coords = np.array(p)
    #
    #
    #     for k in air_contact_coord.keys():
    #         air_contact_coord[k].coords =  air_contact_coord[k].coords.reshape((-1,3))
    #
    #
    #     for k in mat_contact_coord.keys():
    #         mat_contact_coord[k].coords =  mat_contact_coord[k].coords.reshape((-1,3))
    #
    #     return [air_contact_coord, mat_contact_coord]


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

        # if zs == ys and ys == xs:
        #
        #     print 'Eddges have equal dimension ', zs
        #
        #     t0 = tm.time()
        #
        #     zz,yy,xx = np.meshgrid(np.linspace(0, zs-1, zs), np.linspace(0, ys-1, ys), np.linspace(0, xs-1, xs) )
        #
        #     z0 = np.resize(zz, (1,-1))
        #     y0 = np.resize(zz, (1,-1))
        #     z0 = np.resize(zz, (1,-1))
        #
        #     pos0 = np.append(np.append(z0, y0), z0).reshape((-1, 3))
        #     idx = self.main_model[tuple(map(tuple, pos0.transpose()) )] != 0
        #
        #     print 'Time taken to calculate surface coords ', tm.time() - t0, ' seconds'
        #     return pos0[idx, :]

        # pos = np.append(np.resize(zz, (1,-1)))

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
        print 'Time taken ', t1-t0
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
                if m_id in mat_contact.keys():
                    if mat_contact[m_id] == c_id:
                        continue

                if c_id in mat_contact.keys():
                    mat_contact[c_id].id = m_id
                    mat_contact[c_id].coords = np.append(mat_contact[c_id].coords, p)
                else:
                    mat_contact[c_id] = MKSurfaceContact()
                    mat_contact[c_id].id = m_id
                    mat_contact[c_id].coords = np.array(p, dtype=np.int)


    # def isSurface(self, pos = (0,0,0)):
    #
    #     if self.main_model[pos] == 0:
    #         return None, None
    #
    #     if self.system_dimension == 2:
    #
    #         xmin = pos[1] - 1
    #         xmax = pos[1] + 1
    #
    #         if xmin < 0 or xmax >= self.main_model.shape[1]:
    #             return 0, None
    #
    #         ymin = pos[0] - 1
    #         ymax = pos[0] + 1
    #
    #         if ymin <= 0 or ymax >= self.main_model.shape[0]:
    #             return 0, None
    #
    #         if np.prod(self.main_model[ymin:ymax+1, xmin:xmax+1]) == 0:
    #             return 0, None
    #
    #     elif self.system_dimension == 3:
    #         xmin = pos[2] - 1
    #         xmax = pos[2] + 1
    #
    #         if xmin < 0 or xmax >= self.main_model.shape[2]:
    #             return 0, None
    #
    #         ymin = pos[1] - 1
    #         ymax = pos[1] + 1
    #
    #         if ymin <= 0 or ymax >= self.main_model.shape[1]:
    #             return [0, None]
    #
    #
    #         zmin = pos[0] - 1
    #         zmax = pos[0] + 1
    #
    #         if zmin <= 0 or zmax >= self.main_model.shape[0]:
    #             return 0, None
    #
    #         air_cont_mat_id = None
    #         mat_cont_mat_id = None
    #
    #         NN = np.array([[1,0,0], [-1,0,0], [0,1,0] ,[0,-1,0],[0,0,1],[0,0,-1]])
    #
    #         if np.prod(self.main_model[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1]) == 0:
    #             air_cont_mat_id =  self.main_model[pos]
    #
    #
    #         for n in NN:
    #            nn_id = self.main_model[tuple(pos+n)]
    #
    #            if self.main_model[pos] != nn_id :
    #                 mat_cont_mat_id =  nn_id # return the id of the model in contact
    #                 break
    #
    #         return air_cont_mat_id, mat_cont_mat_id
    #     #
    #     #
    #     #     if withAir == True:
    #     #         if np.prod(self.main_model[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1]) == 0:
    #     #             return 1
    #     #         else:
    #     #             return 0
    #     #
    #     #
    #     #     for p_offset in nn:
    #     #         if self.main_model[tuple(pos)] != self.main_model[tuple(pos+p_offset)] :
    #     #             return self.main_model[tuple(pos+p_offset)] # return the id of the model in contact
    #     # else:
    #     #     raise Exception('ERROR: check the dimension of the model (must be 2 or 3)')
    #     #
    #     # return  0


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


    def evolveTempWithTime(self):

        if self.system_dimension == 2: # 0 --> y, 1 --> x

            for e in self.model_objects:
                pos = e[1]
                coords = e[0].surfaceCoordinates + pos
                self.ui[coords[:,0], coords[:,1]] = e[0].dirichlet_bc_value

            self.u[1:-1, 1:-1] = self.ui[1:-1, 1:-1] + self.D[1:-1,1:-1] * self.delta_time * \
                                ( (self.ui[2:, 1:-1] + self.ui[:-2, 1:-1] - 2 * self.ui[1:-1, 1:-1]) /self.delta_pos2[0] + \
                                  (self.ui[1:-1, 2:] + self.ui[1:-1, :-2] - 2 * self.ui[1:-1, 1:-1]) /self.delta_pos2[1] )


        elif self.system_dimension == 3: # 0 --> z, 1 --> y, 2 --> x
            for e in self.model_objects:
                pos = e[1]
                coords = np.array(e[0].surfaceCoordinates + pos, dtype = np.int)
                # self.ui[coords[:,0], coords[:,1], coords[:, 2]] = e[0].dirichlet_bc_value
            #
            # v1 = self.ui[2:, 1:-1, 1:-1]
            # v2 = self.ui[:-2, 1:-1, 1:-1]
            # v3 = self.ui[1:-1, 1:-1, 1:-1]
            #
            # print "v1 shape", v1.shape
            # print "v2 shape", v2.shape
            # print "v3 shape", v3.shape
            #
            # print "D shape ", self.D.shape
            # print "delta_time ", self.delta_time
            # print 'u shape ', self.u.shape
            # print 'ui shape ', self.ui.shape



            # self.u[1:-1, 1:-1, 1:-1] = self.ui[1:-1, 1:-1, 1:-1] + self.D[1:-1, 1:-1, 1:-1] * self.delta_time[1:-1, 1:-1, 1:-1] * \
            #         ( (self.ui[2:, 1:-1, 1:-1] + self.ui[:-2, 1:-1, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[0] + \
            #           (self.ui[1:-1, 2:, 1:-1] + self.ui[1:-1, :-2, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[1] + \
            #           (self.ui[1:-1, 1:-1, 2:] + self.ui[1:-1, 1:-1, :-2] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[2]  \
            #         )

            self.u[1:-1, 1:-1, 1:-1] = self.ui[1:-1, 1:-1, 1:-1] + self.D[1:-1, 1:-1, 1:-1] * self.delta_time * \
                    ( (self.ui[2:, 1:-1, 1:-1] + self.ui[:-2, 1:-1, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[0] + \
                      (self.ui[1:-1, 2:, 1:-1] + self.ui[1:-1, :-2, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[1] + \
                      (self.ui[1:-1, 1:-1, 2:] + self.ui[1:-1, 1:-1, :-2] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[2]  \
                    )

            # self.u[1:-1, 1:-1, 1:-1] = self.ui[1:-1, 1:-1, 1:-1] + self.D[1:-1, 1:-1, 1:-1] * self.delta_time * \
            #         ( (self.ui[2:, 1:-1, 1:-1] + self.ui[:-2, 1:-1, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[0] + \
            #           (self.ui[1:-1, 2:, 1:-1] + self.ui[1:-1, :-2, 1:-1] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[1] + \
            #           (self.ui[1:-1, 1:-1, 2:] + self.ui[1:-1, 1:-1, :-2] - 2 * self.ui[1:-1, 1:-1, 1:-1]) / self.delta_pos2[2] + \
            #           (self.ui[1:-1, 2:] - 2 * self.ui[1:-1, 1:-1] + self.ui[1:-1, :-2]) /self.delta_pos2[0] )

        else:
            raise Exception('ERROR: data must of 2 or 3 dimension')


        self.ui = sp.copy(self.u)

    def calculatesTemperatures(self, tempData):

        for k in self.model_objects.keys():



    def evolveTempWithTime_v2(self):

        # Dirichlet bundary condition with air contact surface
        for k in self.airContactSurface.keys():
            coords = np.array(self.airContactSurface[k].coords, dtype=np.int)
            self.ui[coords[:,0], coords[:,1], coords[:, 2]] = 20
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


        # Dirichlet bundary condition with air contact surface
        for k in self.airContactSurface.keys():
            coords = np.array(self.airContactSurface[k].coords, dtype=np.int)
            self.ui[coords[:,0], coords[:,1], coords[:, 2]] = 20
            # self.u = self.ui

        # Neumann boundary condition along z-axis

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


class MKThermalSimulation(object):
    def __init__(self, complete_model):
        self.complete_model =  complete_model


    def evolve_ts(u, ui):
        """
        This function uses a numpy expression to
        evaluate the derivatives in the Laplacian, and
        calculates u[i,j] based on ui[i,j].
        """
        u[1:-1, 1:-1] = ui[1:-1, 1:-1] + a*dt*( (ui[2:, 1:-1] - 2*ui[1:-1, 1:-1] + ui[:-2, 1:-1])/dx2 + (ui[1:-1, 2:] - 2*ui[1:-1, 1:-1] + ui[1:-1, :-2])/dy2 )

























# class MKModel_test(object):
#     def __init__(self, min_val = (0, 0, 0), max_val = (1, 1, 1), delta_pos = (0.01, 0.01, 0.01)):
#
#         # model parameters are here
#
#         model_dimension = min_val.__len__()
#
#         self.modelPosition = []
#
#         if ~(min_val.__len__() == max_val.__len__() and max_val.__len__() == delta_pos.__len__()):
#             print "All parameters must have the same length of ", model_dimension
#             return
#
#         pv = []
#
#         for idx in xrange(model_dimension):
#             v = np.linspace(min_val[idx], delta_pos[idx], max_val[idx])
#             pv.append(v)
#
#
#         self.mesh_data = []
#
#         if model_dimension == 2:
#             x,y = np.meshgrid(pv[0], pv[1])
#             self.mesh_data.append(x)
#             self.mesh_data.append(y)
#
#             xc = np.mean( x-np.mean(x) )
#             yc = np.mean( y-np.mean(y) )
#
#             self.modelPosition.append(xc)
#             self.modelPosition.append(yc)
#
#         elif model_dimension == 3:
#             x,y,z = np.meshgrid(pv[0], pv[1], pv[2])
#             self.mesh_data.append(x)
#             self.mesh_data.append(y)
#             self.mesh_data.append(z)
#
#             xc = np.mean( x-np.mean(x) )
#             yc = np.mean( y-np.mean(y) )
#             zc = np.mean( z-np.mean(z) )
#
#             self.modelPosition.append(xc)
#             self.modelPosition.append(yc)
#             self.modelPosition.append(zc)


#
# def forceAspect(ax,aspect=1):
#     im = ax.get_images()
#     extent =  im[0].get_extent()
#     ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)