import numpy
import math

class Ellipsis:

    def __init__(self, cells={}):
        """
        Constructor 
        @param cells set of (i,j) tuples
        """
        n = len(cells)
        nm = float(n)

        # centre of the cluster
        self.centre = []
        if n > 0:
            iCentre = numpy.sum([c[0] for c in cells]) / nm
            jCentre = numpy.sum([c[1] for c in cells]) / nm
            self.centre = numpy.array([iCentre, jCentre])

        # compute inertia tensor (symmetric)
        inertia = numpy.zeros((2, 2), numpy.float64)
        for i in range(2):
            for j in range(2):
                inertia[i, j] = numpy.sum([(c[i] - self.centre[i])*(c[j] - self.centre[j]) for c in cells])

        # the set of eigenvectors is the rotation matrix from ij space to the 
        # inertial tensor's principal axes
        eigenvals, self.ij2AxesTransf = numpy.linalg.eig(inertia)

        # from the axes to ij (inverse of the above)
        self.axes2ijTransf = numpy.transpose(self.ij2AxesTransf)

        # average radii from the centre
        self.a = math.sqrt(eigenvals[0] / nm)
        self.b = math.sqrt(eigenvals[1] / nm)


    def getEllipseAsPolyline(self, numSegments=32):
        """
        Return the ellipsis as a segmented line
        @return iPts, jPts arrays
        """
        iPts, jPts = [], []
        dt = 2 * math.pi / float(numSegments)
        for i in range(numSegments + 1):
            th = i * dt
            x = self.a * math.cos(th)
            y = self.b * math.sin(th)
            # rotate back to i,j coordinates
            ij = self.axes2ijTransf.dot([x, y])
            ij += self.centre
            iPts.append(ij[0])
            jPts.append(ij[1])
        return iPts, jPts


    def __repr__(self):
        """
        Print object
        """
        res = """
        Ellpsis: centre = {} a = {} b = {} rotation = {}
        """.format(self.centre, self.a, self.b, self.ij2AxesTransf)
        return res


    def getCentre(self):
        """
        Get the barycentric centre
        """
        return self.centre


    def isPointInside(self, point):
        """
        Check if a point is inside an ellipsis
        @param point point in j, j index space
        @return True if inside, False if outside or on the boundary
        """
        
        # subtract the centre
        point -= self.centre

        # distance of point to the axes
        ptPrimeAbs = abs(self.ij2AxesTransf.dot(point))

        if (ptPrimeAbs[0] < self.a) and (ptPrimeAbs[1] < self.b):
            return True

        return False


#############################################################################################
def test0():
    # test zero set
    ell = Ellipsis({})
    print(ell)

def test1():
    # test zero set
    ell = Ellipsis({(-2, 1)})
    print(ell)

def testRectangle():
    ell = Ellipsis({(i, 0) for i in range(3)}.union({(i, 1) for i in range(3)}))
    # these points should be inside
    pt = numpy.array([1., 0.5])
    assert(ell.isPointInside(pt))
    pt[0] = 1.8; pt[1] = 0.5
    assert(ell.isPointInside(pt))
    pt[0] = 1.; pt[1] = 0.99
    assert(ell.isPointInside(pt))
    # these points should be outside
    pt[0] = 1.82; pt[1] = 0.5
    assert(not ell.isPointInside(pt))
    pt[0] = 1.; pt[1] = 1.01
    assert(not ell.isPointInside(pt))
    
    print(ell)

def testRectangleSlanted():
    ell = Ellipsis({(i, 0) for i in range(4)}.union({(i - 1, 1) for i in range(4)}))
    print(ell)

def testRandom():
    import random
    random.seed(1234)
    ell = Ellipsis({(random.randint(0, 200), random.randint(0, 100)) for i in range(500)})
    print(ell)

if __name__ == '__main__':
    test0()
    test1()
    testRectangle()
    testRectangleSlanted()
    testRandom()


