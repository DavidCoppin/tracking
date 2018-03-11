import numpy
import cv2
from ellipse import Ellipse
from opencv import getContours, getEllipse
import matplotlib.pyplot as mpl
import matplotlib
import math
import time

"""
Compare ellipses
"""

class EllipseTester:
    def __init__(self, label, domain_size, cells):

        self.domain_size = domain_size
        self.label = label
        self.cells = cells

        self.data = {
            'old': {'centre':[None, None], 'axes':[None, None], 'time': None, 'x': [], 'y': []},
            'new': {'centre':[None, None], 'axes':[None, None], 'time': None, 'x': [], 'y': []},
        }

        # new
        tic = time.time()
        ell = Ellipse(cells)
        toc = time.time()
        self.data['new']['time']= toc - tic

        self.data['new']['centre'] = ell.centre
        self.data['new']['axes'] = [ell.a, ell.b]
        self.data['new']['x'], self.data['new']['y'] = ell.getPolyline()

        # old
        cellData = numpy.zeros(domain_size, numpy.int32)
        for c in cells:
            cellData[c] = 1
        tic = time.time()
        test_result, test_contours = getContours(cellData)
        ellipse_parameters = getEllipse(test_result,test_contours)
        toc = time.time()
        self.data['old']['time'] = toc - tic

        # note: reversed order
        yc, xc = ellipse_parameters[0]
        b, a = ellipse_parameters[1]
        b /= 2.0
        a /= 2.0
        angle = ellipse_parameters[2] # angle in degrees

        self.data['old']['centre'] = (xc, yc)
        self.data['old']['axes'] = (a, b)

        # constructor contour
        cos_angle = math.cos(angle*math.pi/180.)
        sin_angle = math.sin(angle*math.pi/180.)
        nt = 32
        ts = numpy.linspace(0., 2*math.pi, nt + 1)
        dxPrime = a*numpy.cos(ts)
        dyPrime = b*numpy.sin(ts)
        dx =  dxPrime * cos_angle + dyPrime * sin_angle
        dy = -dxPrime * sin_angle + dyPrime * cos_angle
        self.data['old']['x'] = xc + dx
        self.data['old']['y'] = yc + dy

    def __repr__(self):
        area_new = math.pi*self.data['new']['axes'][0]*math.pi*self.data['new']['axes'][1]
        area_old = math.pi*self.data['old']['axes'][0]*math.pi*self.data['old']['axes'][1]
        res = """
EllipseTester {}
center (old new): {}    ; {}
axes   (old new): {}    ; {}
area   (old new): {}    ; {}
times  (old new): {}s   ; {}s for domain size {}
        """.format(self.label,
                   self.data['old']['centre'], self.data['new']['centre'],
                   self.data['old']['axes'], self.data['new']['axes'],
                   area_old, area_new,
                   self.data['old']['time'], self.data['new']['time'],
                   self.domain_size)
        return res

    def plot(self):

        fig, ax = mpl.subplots()
        
        # plot the cells
        patches = []
        for c in self.cells:
            i, j = c
            pts = numpy.array([[i-0.5, j-0.5], [i+0.5, j-0.5], [i+0.5, j+0.5], [i-0.5, j+0.5]])
            patch = matplotlib.patches.Polygon(pts, closed=True)
            patches.append(patch)
        p = matplotlib.collections.PatchCollection(patches, alpha=0.4)
        ax.add_collection(p)

        # plot the old ellipse
        mpl.plot(self.data['old']['x'], self.data['old']['y'], 'g-')

        # plot the new ellipse
        mpl.plot(self.data['new']['x'], self.data['new']['y'], 'r-')

        # add labels
        mpl.legend(['old', 'new'])

        # add the centres
        mpl.plot([self.data['old']['centre'][0]], [self.data['old']['centre'][1]], 'gx')
        mpl.plot([self.data['new']['centre'][0]], [self.data['new']['centre'][1]], 'r+')

        mpl.xlabel('i')
        mpl.ylabel('j')
        mpl.title(self.label)

        # show the result
        mpl.show()



#####################################################################
def testObliqueRectangle():

    domain_size = (30, 40)
    cells = set()
    for j in range(10):
        cells = cells.union({(i+j, i) for i in range(3, 14)})
    et = EllipseTester('oblique rectangle', domain_size, cells)
    print et
    et.plot()


if __name__ == '__main__':
    testObliqueRectangle()