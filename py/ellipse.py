import numpy
import math
from shapely.geometry.point import Point
from shapely import affinity
from matplotlib.patches import Polygon
"""
A Class that computes the ellipse of a cloud of points
"""

class Ellipse:

    def __init__(self, cells, min_ellipse_axis=10):
        """
        Constructor 
        @param cells set of (i,j) tuples, must have at least one cell
        @param min_ellipse_axis min axis length 
        """
        n = len(cells)
        area = float(n)

        inertia = numpy.zeros((2, 2), numpy.float64)

        # centre of cluster
        self.centre = []

        # average radii from the centre
        self.a = 0.
        self.b = 0.

        iInds = numpy.array([c[0] for c in cells], numpy.float64)
        jInds = numpy.array([c[1] for c in cells], numpy.float64)

        iCentre = iInds.sum() / area
        jCentre = jInds.sum() / area
        self.centre = numpy.array([iCentre, jCentre])

        # compute inertia tensor (symmetric)
        iInds -= iCentre
        jInds -= jCentre
        inertia[0, 0] = numpy.sum(iInds * iInds)
        inertia[0, 1] = numpy.sum(iInds * jInds)
        inertia[1, 0] = inertia[0, 1]
        inertia[1, 1] = numpy.sum(jInds * jInds)

        # the set of eigenvectors is the rotation matrix from ij space to the 
        # inertial tensor's principal axes
        eigenvals, self.axes2ijTransf = numpy.linalg.eig(inertia)
        self.ij2AxesTransf = numpy.transpose(self.axes2ijTransf)

        # angle between the principal axes and the i, j directions
        self.angle = math.atan2(self.ij2AxesTransf[0, 1], self.ij2AxesTransf[0, 0])*180./numpy.pi

        # average radii from the centre
        a, b = numpy.sqrt(eigenvals)
        self.a = max(0.5, a)
        self.b = max(0.5, b)

        # compute the total area
        area = len(cells)

	    # extend the axis to match the cluster's area
        const = math.sqrt(area /(math.pi * self.a * self.b))
        a *= const
        b *= const
        self.a *= const
        self.b *= const

        # add halo to the axes if need be
        #self.aExt = max(min_ellipse_axis, self.a)
        #self.bExt = max(min_ellipse_axis, self.b)
        # smooth transition
        amin = min_ellipse_axis
        self.aExt = self.a + amin*math.exp(-a/amin)
        self.bExt = self.b + amin*math.exp(-b/amin)


    def createEllipse(self, centre, a, b, angle):
        """
        create a shapely ellipse
        """
        circ = Point(centre).buffer(1)
        ell = affinity.scale(circ, int(a), int(b))
        ellr = affinity.rotate(ell, angle)
        return ellr


    def getPolyline(self, numSegments=32, a=None, b=None):
        """
        Return the ellipse as a segmented line
        @param numSegments number of segments
        @param a principal axis length, default to self.a
        @param b second axis length, defaults to self.b
        @return iPts, jPts arrays
        """
        if not a:
            a = self.a
        if not b:
            b = self.b
        iPts, jPts = [], []
        dt = 2 * math.pi / float(numSegments)
        for i in range(numSegments + 1):
            th = i * dt
            x = a * math.cos(th)
            y = b * math.sin(th)
            # rotate back to i,j coordinates
            ij = self.axes2ijTransf.dot([x, y])
            ij += self.centre
            iPts.append(ij[0])
            jPts.append(ij[1])
        return iPts, jPts


    def getPolylineExt(self, numSegments=32):
        """
        Return the extended ellipse as a segmented line
        @param numSegments number of segments
        @return iPts, jPts arrays
        """
        return self.getPolyline(numSegments=numSegments, a=self.aExt, b=self.bExt)


    def __repr__(self):
        """
        Print object
        """
        res = """
        Ellipse: centre = {} a = {} b = {} rotation = {} angle = {}
        """.format(self.centre, self.a, self.b, self.ij2AxesTransf, self.angle)
        return res


    def getCentre(self):
        """
        Get the barycentric centre
        """
        return self.centre


    def isPointInside(self, point):
        """
        Check if point is inside ellipse
        @param point point in j, j index space
        @return True if inside, False if outside or on the boundary
        """

        # rotate the coordinates to align them to the principal axes
        ptPrimeAbs = self.ij2AxesTransf.dot(point - self.centre)
        eps = 1.e-12

        ptPrimeAbs[0] /=  self.a + eps
        ptPrimeAbs[1] /=  self.b + eps
##        if (ptPrimeAbs[0]/(self.a + eps))**2 + (ptPrimeAbs[1]/(self.b + eps))**2 < 1.0:
        if (ptPrimeAbs[0]*ptPrimeAbs[0] + ptPrimeAbs[1]*ptPrimeAbs[1]) < 1.0:
            # inside
            return True

        return False


    def isPointInsideExt(self, point):
        """
        Check if point is inside extended ellipse
        @param point point in j, j index space
        @return True if inside, False if outside or on the boundary
        """
  
        # rotate the coordinates to align them to the principal axes
        ptPrimeAbs = self.ij2AxesTransf.dot(point - self.centre)

        ptPrimeAbs[0] /= self.aExt
        ptPrimeAbs[1] /= self.bExt
##        if (ptPrimeAbs[0]/self.aExt)**2 + (ptPrimeAbs[1]/self.bExt)**2 < 1.0:
        if (ptPrimeAbs[0]*ptPrimeAbs[0] + ptPrimeAbs[1]*ptPrimeAbs[1]) < 1.0:
            # inside
            return True

        return False


    def isEllipseInsideOf(self, otherEllipse, frac):
        """
        Check if fraction of this cluster is inside otherCluster
        @param otherEllipse
        @param frac fraction of min area of self and otherCluster
        """
        ellipse_1 = self.createEllipse(self.centre, self.a, self.b, self.angle)
        if ellipse_1.is_valid == False:
            ellipse1 = ellipse_1.buffer(0)
        else:
            ellipse1 = ellipse_1
        if ellipse1.is_empty:
            return False
        ellipse_2 = self.createEllipse(otherEllipse.centre, otherEllipse.a, otherEllipse.b, otherEllipse.angle)
        if ellipse_2.is_valid == False:
            ellipse2 = ellipse_2.buffer(0)
        else:
            ellipse2 = ellipse_2
        if ellipse2.is_empty:
            return False
#       print 'ellipse1.is_valid', ellipse1.is_valid
#        print 'ellipse1.buffer(0)', ellipse1.buffer(0), ellipse1.buffer(0).is_valid
#       print 'ellipse2.is_valid', ellipse2.is_valid
        intersect = ellipse1.intersection(ellipse2)
        min_area = min(ellipse1.area, ellipse2.area)
#        print 'intersect.area, intersect.area/min_area', intersect.area, intersect.area/min_area
        return intersect.area/min_area >= frac


    def show(self, points=[], cells={}, show=True):
        """
        Plots the ellipse
        @param points set of points to be shown as inside (stars) or outside (x)
        """
        from matplotlib import pylab
        import matplotlib

        fig, ax = pylab.subplots()

        # show the cells
        patches = []
        for c in cells:
            i, j = c
            pts = numpy.array([[i-0.5, j-0.5], [i+0.5, j-0.5], [i+0.5, j+0.5], [i-0.5, j+0.5]])
            patch = matplotlib.patches.Polygon(pts, closed=True)
            patches.append(patch)
        p = matplotlib.collections.PatchCollection(patches, alpha=0.4)
        ax.add_collection(p)

        # plot the ellipse
        iPts, jPts = self.getPolyline()
        pylab.plot(iPts, jPts, 'r-')
        iPts, jPts = self.getPolylineExt()
        pylab.plot(iPts, jPts, 'b--')

        # add points (stars if inside, crosses if outside)
        pointsInside = []
        pointsOutside = []
        for p in points:
            if self.isPointInside(p):
                pointsInside.append(p)
            else:
                pointsOutside.append(p)
        pylab.plot([p[0] for p in pointsOutside], [p[1] for p in pointsOutside], 'kx')
        pylab.plot([p[0] for p in pointsInside], [p[1] for p in pointsInside], 'cs')

        # label, title, ...
        pylab.xlabel('i')
        pylab.ylabel('j')
        if show:
            pylab.show()

#############################################################################################
def test0():
    # test zero set
    ell = Ellipse({})
    print(ell)

def test1():
    # test zero set
    ell = Ellipse({(-2, 1)})
    print(ell)

def testRectangle():
    ell = Ellipse({(i, 0) for i in range(3)}.union({(i, 1) for i in range(3)}))
    print(ell)

    pts = []
    # these points should be inside
    pt = numpy.array([1., 0.5])
    assert(ell.isPointInside(pt))
    pts.append(pt)

    pt = numpy.array([1.8, 0.5])
    assert(ell.isPointInside(pt))
    pts.append(pt)

    pt = numpy.array([1., 0.99])
    assert(ell.isPointInside(pt))
    pts.append(pt)

    # these points should be outside
    pt = numpy.array([1.82, 0.5])
    #assert(not ell.isPointInside(pt))
    pts.append(pt)

    pt = numpy.array([1., 1.01])
    #assert(not ell.isPointInside(pt))
    pts.append(pt)

    #ell.show(pts)


def testRectangleSlanted():
    import random
    random.seed(1234)

    cells = {(i, 0) for i in range(4)}.union({(i - 1, 1) for i in range(4)})
    ell = Ellipse(cells)
    print(ell)

    # create lots of random points
    pts = [numpy.array([-2. + 6*random.random(), -1. + 3*random.random()]) for i in range(1000)]

    ell.show(pts, cells)


def testRandom():
    import random
    random.seed(1234)
    ell = Ellipse({(random.randint(0, 200), random.randint(0, 100)) for i in range(500)})
    print(ell)


def testMinEllipseAxis(axis=1):
    cells = {(i, 0) for i in range(4)}.union({(i - 1, 1) for i in range(4)})
    ell = Ellipse(cells, min_ellipse_axis=axis)
    print 'testMinEllipseAxis angle (deg) ', ell.angle
    print 'transf ij -> axes ', ell.ij2AxesTransf
    print 'ellipse axes: ', ell.a, ell.b
    print 'ellipse ext axes: ', ell.aExt, ell.bExt
    ell.show(cells=cells)
    

def testMinEllipseAreaBig():
    cells = {(i, 0) for i in range(4)}.union({(i - 1, 1) for i in range(4)})
    ell = Ellipse(cells, min_ellipse_axis=120)
    ell.show()

def testExt45Deg():
    """
    Check that the extended ellipse has the same pricipal axes directions as the
    original ellipse
    """
    # very flat oblique cells
    cells = {(i, i) for i in range(20)}
    ell = Ellipse(cells, min_ellipse_axis=5)
    print 'angle should about 45 +- 90 deg: ', ell.angle
    ell.show()

def testExt63Deg():
    """
    Check that the extended ellipse has the same pricipal axes directions as the
    original ellipse
    """
    # very flat oblique cells
    cells = {(i, 2*i) for i in range(20)}
    ell = Ellipse(cells, min_ellipse_axis=5)
    print 'angle should about 63 +- 90 deg: ', ell.angle
    ell.show()

def testExtHorizontal():
    cells = {(i, 10) for i in range(20)}
    ell = Ellipse(cells, min_ellipse_axis=5)
    print 'angle should about 0 +- 90 deg: ', ell.angle
    ell.show()

def testExtVertical():
    cells = {(10, i) for i in range(20)}
    ell = Ellipse(cells, min_ellipse_axis=5)
    print 'angle should about 90 +- 90 deg: ', ell.angle
    ell.show()

def testExtLeftOblique():
    cells = {(50 - i, 2*i) for i in range(20)}
    ell = Ellipse(cells, min_ellipse_axis=5)
    print 'angle should about 117 +- 90 deg: ', ell.angle
    ell.show()




if __name__ == '__main__':
    testExt45Deg()
    testExt63Deg()
    testExtHorizontal()
    testExtVertical()
    testExtLeftOblique()
    testMinEllipseAxis(axis=1)
    testMinEllipseAxis(axis=2)
    testMinEllipseAxis(axis=5)
    testMinEllipseAxis(axis=10)
    testMinEllipseAxis(axis=100)
    test0()
    test1()
    testRectangle()
    testRectangleSlanted()
    testRandom()
