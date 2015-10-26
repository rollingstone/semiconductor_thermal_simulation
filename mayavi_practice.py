__author__ = 'kamal'


import numpy
from mayavi.mlab import *

def test_mesh():
    """A very pretty picture of spherical harmonics translated from
    the octaviz example."""
    pi = numpy.pi
    cos = numpy.cos
    sin = numpy.sin
    dphi, dtheta = pi / 250.0, pi / 250.0
    [phi, theta] = numpy.mgrid[0:pi + dphi * 1.5:dphi,
                               0:2 * pi + dtheta * 1.5:dtheta]
    m0 = 4
    m1 = 3
    m2 = 2
    m3 = 3
    m4 = 6
    m5 = 2
    m6 = 6
    m7 = 4
    r = sin(m0 * phi) ** m1 + cos(m2 * phi) ** m3 + \
        sin(m4 * theta) ** m5 + cos(m6 * theta) ** m7
    x = r * sin(phi) * cos(theta)
    y = r * cos(phi)
    z = r * sin(phi) * sin(theta)

    return mv.mesh(x, y, z, colormap="hot")



def test_points3d():
    t = numpy.linspace(0, 4 * numpy.pi, 20)
    cos = numpy.cos
    sin = numpy.sin

    x = sin(2 * t)
    y = cos(t)
    z = cos(2 * t)
    s = 2 + sin(t)

    return points3d(x, y, z, s, colormap="copper", scale_factor=.25)



m = test_points3d()



# Create the data.
from numpy import pi, sin, cos, mgrid
dphi, dtheta = pi/250.0, pi/250.0
[phi,theta] = mgrid[0:pi+dphi*1.5:dphi,0:2*pi+dtheta*1.5:dtheta]
m0 = 4; m1 = 3; m2 = 2; m3 = 3; m4 = 6; m5 = 2; m6 = 6; m7 = 4;
r = sin(m0*phi)**m1 + cos(m2*phi)**m3 + sin(m4*theta)**m5 + cos(m6*theta)**m7
x = r*sin(phi)*cos(theta)
y = r*cos(phi)
z = r*sin(phi)*sin(theta)

# View it.
from mayavi import mlab
s = mlab.mesh(x, y, z)
mlab.show()



import numpy
from mayavi.mlab import *

def test_contour_surf():
    """Test contour_surf on regularly spaced co-ordinates like MayaVi."""
    def f(x, y):
        sin, cos = numpy.sin, numpy.cos
        return sin(x + y) + sin(2 * x - y) + cos(3 * x + 4 * y)

    x, y = numpy.mgrid[-7.:7.05:0.1, -5.:5.05:0.05]
    s = contour_surf(x, y, f)
    show()
    return s

test_contour_surf()





