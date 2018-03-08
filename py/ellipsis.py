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

        # diagonalize
        self.eigenvals, eigenvecs = numpy.linalg.eig(inertia)

        # number of masses
        nm = float(len(cells))

        # average radii from the centre
        self.a = numpy.sqrt(eigenvals[0] / nm)
        self.b = numpy.sqrt(eigenvals[1] / nm)

        # orientation given by the first eigenvector (the other must 
        # be perpendicular)
        vecA = eigenvecs[:, 0]
        self.angle = math.atan2(vecA[1], vecA[0])

        cosa = math.cos(angle)
        sina = math.sin(angle)

        # transformation from principle axes directions to i,j
        self.axes2ijTransf = numpy.array([[cosa, -sina], [sina, cosa]])

        # transformation from i, j to principle axes directions
        self.ij2AxesTransf = numpy.array([[cosa, sina], [-sina, cosa]])

    def getCentre(self):
    	return self.centre

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
        Ellpsis: centre = {} a = {} b = {} angle = {}
        """.format(self.centre, self.a, self.b, self.angle)
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

    def isPointInsideEllipse(self, point, 
               ellipsisCentre, ellipsisA, ellipsisB, ellipsisAngle):
        """
        Check if a point is inside an ellipse
        @return True if inside, False if outside or on the boundary
        """
        
        # subtract the centre
        point -= self.centre

        # distance of point to the axes
        ptPrimeAbs = abs(self.ij2AxesTransf.dot(point))

        return False


#############################################################################################
def test0():
	pass


if __name__ == '__main__':
    test0()

