import numpy
import math

class Ellipsis:

    def __init__(self, cells={}):
        """
        Constructor 
        @param cells set of (i,j) tuples
        """
        self.centre = self.getCentre()
        inertia = numpy.zeros((2, 2), numpy.float64)
        inertia[0, 0] = numpy.sum([(c[0] - self.centre[0])**2 for c in cells])
        inertia[0, 1] = numpy.sum([(c[0] - self.centre[0])*(c[1] - self.centre[1]) for c in cells])
        inertia[1, 0] = inertia[0, 1]
        inertia[1, 1] = numpy.sum([(c[1] - self.centre[1])**2 for c in cells])

        # the set of eigenvectors is the rotation matrix from ij space to the 
        # principal axes
        eigenvals, self.ij2AxesTransf = numpy.linalg.eig(inertia)

        # from the axes to ij
        self.axes2ijTransf = numpy.transpose(self.ij2AxesTransf)

        # average radii from the centre
        nm = float(len(cells))
        self.a = numpy.sqrt(eigenvals[0] / nm)
        self.b = numpy.sqrt(eigenvals[1] / nm)

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
        iCentre, jCentre = None, None
        n = len(self.cells)
        if n > 0:
            iCentre = numpy.sum([c[0] for c in self.cells]) / float(n)
            jCentre = numpy.sum([c[1] for c in self.cells]) / float(n)
        return numpy.array([iCentre, jCentre])


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
	pass


if __name__ == '__main__':
    test0()

