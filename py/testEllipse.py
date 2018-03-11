import numpy
import cv2
from ellipse import Ellipse
from opencv import getContours, getEllipse
import matplotlib.pyplot as mpl
import math

"""
Compare the ellipses obtained from Ellipse with those obtained from opencv3
"""

# define the cells where something special happens
cells = set()
for j in range(10):
    cells = cells.union({(i+j, i) for i in range(3, 14)})

# turn this into a numpy array
domain = (30, 40)
cellData = numpy.zeros(domain, numpy.int32)
# set the data to 1 in the above cells
for c in cells:
    cellData[c] = 1


# create the ellipse
ell = Ellipse(cells)

# now do the same with cv2...
test_result, test_contours = getContours(cellData)
ellipse_parameters = getEllipse(test_result,test_contours)

yc, xc = ellipse_parameters[0]
b, a = ellipse_parameters[1]
b /= 2.0
a /= 2.0
angle = ellipse_parameters[2] # angle in degrees?
print 'ellipse centre: {} {}'.format(xc, yc)
print 'ellipse axes  : {} {}'.format(a, b)
print 'ellipse angle : {}'.format(angle)
print 'ellipse area  : {}'.format(math.pi * a * b)
print 'cell area     : {}'.format(len(cells))

cos_angle = math.cos(angle*math.pi/180.)
sin_angle = math.sin(angle*math.pi/180.)

print ell
ell.show(points=[], cells=cells, show=False)

# show the ellipse from opencv2
nt = 32
ts = numpy.linspace(0., 2*math.pi, nt + 1)
dxPrime = a*numpy.cos(ts)
dyPrime = b*numpy.sin(ts)
dx =  dxPrime * cos_angle + dyPrime * sin_angle
dy = -dxPrime * sin_angle + dyPrime * cos_angle
mpl.plot(xc + dx, yc + dy, 'm--')


#mpl.contourf(test_result)
#
mpl.show()

# compare the two ellipses
# centre, axes, angle.... may be points?

