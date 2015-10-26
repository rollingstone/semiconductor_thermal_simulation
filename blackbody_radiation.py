__author__ = 'kamal'

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt


def stefan_radiation(A, T, e = 0.5):
    return (e * A * (T+273)**4 * 5.67e-8)



def newton_cooling(A = 1, dT = 10, hc = 50):
    return hc * A * dT





A = 12e-8 # m2
V = 4e6*(1e-6 *100)**3  # cm3

rho = 5.3 # gm/cm3
m = rho * V # gm
C = 0.33 # J/g/C

Cm =  C*m


t0 =  stefan_radiation(A, 0)
t25 =  stefan_radiation(A, 25)
t27 =  stefan_radiation(A, 27)
t50 =  stefan_radiation(A, 50)
t100 = stefan_radiation(A, 100)
t200 = stefan_radiation(A, 200)

print 'Stefans law...'

print 't0 ', t0, 'W'
print 't25 ', t25, 'W'
print 't27 ', t27, 'W'
print 't50 ', t50, 'W'
print 't100 ', t100, 'W'
print 't200 ', t200, 'W'
print ''

print 'Difference...'
print 't25-t25', t27-t25
print 't27-t0', t27-t25
print 't50-t0', t50-t25
print 't100-t0', t100-t25
print 't200-t0', t200-t25



print 'Newton cooling...'
print ''

ct27 =  newton_cooling(A=A, dT=27)
ct50 =  newton_cooling(A=A,  dT=50)
ct100 = newton_cooling(A=A,  dT=100)
ct200 = newton_cooling(A=A,  dT=200)

print 'ct27 ', ct27, 'W'
print 'ct50', ct50, 'W'
print 'ct100', ct100, 'W'
print 'ct200', ct200, 'W'


temp = np.linspace(25,1000.0,51)

bb_radiation =  stefan_radiation(A, temp) - t25
convection =  newton_cooling(A=A, dT=temp-25)

print 'bb_radiation ', bb_radiation



plt.figure(101)
plt.clf()
plt.cla()
ax = plt.subplot()
ax.plot(temp, bb_radiation)
plt.yscale('log')
plt.show()





plt.figure(102)
plt.clf()
plt.cla()
ax = plt.subplot()
ax.plot(temp, convection)

plt.yscale('log')
plt.show()



plt.figure(103)
plt.clf()
plt.cla()
ax = plt.subplot()
h1, = ax.plot(temp, bb_radiation, color='red', linewidth = 3, label=  'Blackbody Radiation')
h2, = ax.plot(temp, convection, color = 'blue', linewidth = 3, label = 'Solid to Air heat convection')
plt.legend(handles = [h1,h2], loc=4)
plt.ylabel('Power dissipation (J/s)')
plt.xlabel('Temperature (C)')
plt.yscale('log')
plt.title('Radiation J/s vs Temperature of LED (GaAs)')
plt.show()



plt.figure(104)
plt.clf()
plt.cla()
ax = plt.subplot()
h1, = ax.plot(temp[1:], convection[1:]/bb_radiation[1:], color='red', linewidth = 3, label=  'Convection/Black Body Radiation')
# h2, = ax.plot(temp, convection, color = 'blue', linewidth = 3, label = 'Solid to Air heat convection')
plt.legend(handles = [h1])
plt.ylabel('Convection/Blackbody Radiation')
plt.xlabel('Temperature (C)')
# plt.yscale('log')
plt.show()


plt.figure(105)
plt.clf()
plt.cla()
ax = plt.subplot()
h1, = ax.plot(temp, bb_radiation/Cm, color='red', linewidth = 3, label=  'Blackbody Radiation')
h2, = ax.plot(temp, convection/Cm, color = 'blue', linewidth = 3, label = 'Solid to Air heat convection')
plt.legend(handles = [h1,h2], loc=4)
plt.ylabel('Temperature dissipation rate(C/s)')
plt.xlabel('Temperature (C)')
plt.yscale('log')
plt.title('Rate of Temperature Change vs Temperature of LED (GaAs)')
plt.show()




h = 50 # for water


