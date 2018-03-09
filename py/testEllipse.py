import numpy
import cv2
from ellipse import Ellipse

"""
Compare the ellipses obtained from Ellipse with those obtained from opencv3
"""

# define the cells where something special happens
cells = set()
for j in range(10):
    cells = cells.union({(i+j, i) for i in range(14)})

# turn this into a numpy array
domain = (30, 40)
cellData = numpy.zeros(domain, numpy.int32)
# set the data to 1 in the above cells
for c in cells:
    cellData[c] = 1


# create the ellipse
ell = Ellipse(cells)

# now do the same with cv2...
# TO DO 

# compare the two ellipses
# centre, axes, angle.... may be points?

