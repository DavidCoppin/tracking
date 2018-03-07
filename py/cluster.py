import numpy

class Cluster:

    def __init__(self, cells={}):
        """
        Constructor 
        @param cells set of (i,j) tuples
        """
        self.cells = cells
        self.centre = []
        self.ellipsisA = None
        self.ellipsisB = None
        self.ellipsisAngle = None
        self.iMin, self.jMin, self.iMax, self.jMax = None, None, None, None
    
    def update(self):
        self.centre, self.ellipsisA, self.ellipsisB, self.angle = self.getEllipsis()
        self.iMin = numpy.minimum([c[0] for c in self.cells])
        self.jMin = numpy.minimum([c[1] for c in self.cells])
        self.iMax = numpy.maximum([c[0] for c in self.cells])
        self.jMax = numpy.minimum([c[1] for c in self.cells])

    def merge(self, otherCluster):
        self.cells.union(otherCluster.cells)
        self.update()

    def overlaps(self, otherCluster):

        # DAVID TO VERIFY!!!!

        # quick check...
        res = False
        if otherCluster.iMax < self.iMin or otherCluster.iMin > self.iMax:
            return res
        if otherCluster.jMax < self.jMin or otherCluster.jMin > self.jMax:
            return res

        # the clusters overlap if the ellipses overlap
        


        return True

         
    def write_file(self, filename):
        pass

    def getNumCells(self):
        return len(self.cells)

    def __getCentre(self):
        """
        Get the barycentric centre
        """
        iCentre, jCentre = None, None
        n = len(self.cells)
        if n > 0:
            iCentre = numpy.sum([c[0] for c in self.cells]) / float(n)
            jCentre = numpy.sum([c[1] for c in self.cells]) / float(n)
        return numpy.array([iCentre, jCentre])

    def getEllipsis(self):

        centre = self.__getCentre()
        inertia = numpy.zeros((2, 2), numpy.float64)
        inertia[0, 0] = numpy.sum([(c[0] - centre[0])^2 for c in self.cells])
        inertia[0, 1] = numpy.sum([(c[0] - centre[0])*(c[1] - centre[1]) for c in self.cells])
        inertia[1, 0] = inertia[0, 1]
        inertia[1, 1] = numpy.sum([(c[1] - centre[1])^2 for c in self.cells])

        # diagonalize
        eigenvals, eigenvecs = numpy.linalg.eig(intertia)
        a = numpy.sqrt(eigenvals[0])
        b = numpy.sqrt(eigenvals[1])
        vecA = eigenvecs[:, 0]
        angle = math.atan2(vecA[1], vecA[0])

        return centre, a, b, angle



#############################################################################################
def test1():
    b1 = Blob(1, {(1,)})
    b2 = Blob(1, {(2,)})
    assert(b1.overlaps(b1))
    assert(not b1.overlaps(b2))
    assert(b1.shift(0, 1).overlaps(b2))
    assert(not b1.shift(0, -1).overlaps(b2))
    assert(b1.merge(b2)) # yes can be merged
    assert(b1.getNumCells() == 2)

if __name__ == '__main__':
    test1()

