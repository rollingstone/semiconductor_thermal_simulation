__author__ = 'kamal'

import numpy as np
from scipy.optimize import curve_fit



np.divide(2.0, 4.0)

x1 = np.arange(9.0).reshape((3, 3))

x2 = np.arange(3.0)

np.divide(x1, x2)



ld_err_state = np.seterr(divide='ignore')
vv = np.array([1,3,2],dtype=np.float) / np.array([2,0,1])

print vv

vv[np.isinf(vv)] = 0

print vv







for i in xrange(5,10):
    print(i)

###############################


dd = np.linspace(0, 4, 5)
x1,y1,z1 = np.meshgrid(dd,dd,dd)

ps = [z1.reshape((-1,1))]


pp = np.vstack(zip(z1,y1,x1))
pp


p2 = np.tile(z1, (1,-1))
p2


p3 = np.resize(x1, (1,-1))
p2 = np.resize(x1, (1,-1))

pp = np.vstack((p3,p2)).reshape((2,-1)).transpose()
pp


pp = np.array([[1,2,3],[4,5,6],[7,8,9]])

for p0 in pp:
    print 'p0 = ', p0


# T(1,j+1)=T(1,j)+lambda*( T(2,j)- 2*T(1,j) + T(2,j));
# T(m+1,j+1)=T(m+1,j)+lambda*(T(m,j)-2*T(m+1,j)+(T(m,j)+2*dx*2));
import numpy as np


v1 = np.array([[1,2,3],[3,4,5],[4,2,np.nan]])
v1

v2 = np.array([[1,2,3],[3,np.nan,5],[4,2,1]])
v2


v3 = v1 + v2
v3


vm1 = np.nanmean(v3, axis=0)
vm1





idx = v1[:] == 2
idx

v1
v1[idx] = 99



dd = {}

dd[1] = 100
dd[2] = 200
dd[3] = 300

for d, v in dd.iteritems():
    print 'key ', d
    print 'value ', v


v10 = np.array([[1,2,3,4,5,6],[7,8,9,10,11,12],[13,14,15,16,17,18],[19,20,21,22,23,24]])
v10 = np.vstack(tuple(v10))
v10



import numpy as np
from numba import vectorize, float64, void


@vectorize([float64(float64,float64, float64)])
def vfunc(A,B, val):
    return np.nan_to_num(A * B + 100 * A/10.5) * val

def nfunc(A,B, val):
    return  np.nan_to_num(A * B + 100 * A/10.5) * val


a = np.random.rand(10000,10000)
b = np.random.rand(10000,10000)


v1 = np.linspace(0, 2, 3)
v2 = np.linspace(0, 1, 2)
v3 = np.linspace(0, 2, 3)

vv1,vv2, vv3 = np.meshgrid(v1, v2, v3)

vv1
vv2



import matplotlib.pyplot as plt
import numpy as np
import pickle as pk

ax = plt.subplot(111)
x = np.linspace(0, 10)
y = np.sin(x)
plt.plot(x, y)
pk.dump(ax, file('myplot.pickle', 'w'))
# Then in a separate session:


ax = pk.load(file('myplot.pickle'))
plt.show()




test_data = np.array([1,2,3,4,5,6,7])
d0 = test_data*100
np.savez('test_data_save.npz', test_data = test_data, d0 =  d0)


vv = np.load('test_data_save.npz')



def func1(x, a, c, d):
    return a*np.exp(-c*x)+d

x = 10*np.linspace(0,100,101)/100
y = np.exp(-x)

popt, pcov = curve_fit(func1, x, y, p0=(1, 1e-6, 1))


plt.figure(91)
plt.clf()
plt.cla()
plt.plot(x,y, color='red', marker = 'o')
plt.plot(x, func1(x, *popt), color='blue')
plt.show()






v1 = np.linspace(0,3,4)
v2 = np.linspace(0,5,6)
v3 = np.linspace(0,4,5)


zz,yy,xx = np.meshgrid(v3,v2,v1)


xx = xx.flatten()
yy = yy.flatten()
zz = zz.flatten()

pos = np.array([zz,yy,xx]).transpose()
pos


p2 = np.reshape(pos, (3,-1)).transpose()

pos2 = np.array([])

for iz in xrange(4):
    for iy in xrange(6):
        for ix in xrange(4):
            pos2 = np.append(pos2, [iz,iy,ix])

pos2 = np.reshape(pos2,(-1,3))



# pos = np.reshape(np.array([zz,yy,xx]),(-1,3))



from multiprocessing import Pool


def test_func(a,b,c):
    return a**10+b*c


p = Pool(processes=5)

result = p.apply_async(test_func, args=[(1,2,3)])

print result.get(timeout=2)

p.close()


def mk_func(a,b):

    print a
    print b

    a[1] = 999
    b[2] = 10000000



a = np.array([1,2,3,4,5,6])
b = np.array([10,11,12,13,14,15,16])

mk_func(a,b)

print a
print b


from multiprocessing import Pool



def f(x):
    return x*x

if __name__ == '__main__':
    pool = Pool(processes=4)              # start 4 worker processes

    print pool.map(f, [1,2,3])

    # result = pool.apply_async(f, [10])    # evaluate "f(10)" asynchronously
    # print result.get(timeout=1)           # prints "100" unless your computer i
    pool.close()



L = 10*1e-6 *100 # cm
vx_volume = L**3 ## cm3

r = 5.31 # g/cm3
c = 0.33 # J/g/C

i_val = 16e-3 # amp
v_val = 2.5 # volt

q_val = (i_val * v_val) /1.0

heat_0 = (q_val)/(r*c)
heat_vx = (q_val/vx_volume)/(r*c)

plan_vol = 20*20*vx_volume
total_vol = 20*20*10*vx_volume

heat_plane = (q_val/plan_vol)/(r*c)
heat_volume = (q_val/total_vol)/(r*c)



# heat_0 = (q_val)/(r*c)





import matplotlib.pyplot as plt
plt.ion()
class DynamicUpdate():
    #Suppose we know the x range
    min_x = 0
    max_x = 10

    def on_launch(self):
        #Set up plot
        self.figure, self.ax = plt.subplots()
        self.lines, = self.ax.plot([],[], 'o')
        #Autoscale on unknown axis and known lims on the other
        self.ax.set_autoscaley_on(True)
        self.ax.set_xlim(self.min_x, self.max_x)
        #Other stuff
        self.ax.grid()
        ...

    def on_running(self, xdata, ydata):
        #Update data (with the new _and_ the old points)
        self.lines.set_xdata(xdata)
        self.lines.set_ydata(ydata)
        #Need both of these in order to rescale
        self.ax.relim()
        self.ax.autoscale_view()
        #We need to draw *and* flush
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    #Example
    def __call__(self):
        import numpy as np
        import time
        self.on_launch()
        xdata = []
        ydata = []
        for x in np.arange(0,10,0.5):
            xdata.append(x)
            ydata.append(np.exp(-x**2)+10*np.exp(-(x-7)**2))
            self.on_running(xdata, ydata)
            time.sleep(1)
        return xdata, ydata

d = DynamicUpdate()
d()


x = np.linspace(0,3,4)
y = np.linspace(0,2,3)
z = np.linspace(0,5,6)

mat = np.zeros((6,3,4))

zz,yy,xx = np.meshgrid(z,y,x)

xf = xx.flatten().astype(np.int)
yf = yy.transpose().flatten().astype(np.int)
zf = zz.flatten().astype(np.int)

idx = np.logical_and(zf > 2 , zf < 4)


mat[zf[idx], :, :] = 100


v = np.array([zf,yf,xf]).flatten().reshape((3,-1)).transpose()
# v2 = np.array([zf,yf,xf]).reshape((3,-1)).transpose()
print v


idx2 = mat[v[:,0], v[:,1], v[:,2]] == 100

v2 = v[idx2,:]

print 'v2   '
print v2



















