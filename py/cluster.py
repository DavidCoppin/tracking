import numpy

class Cluster:

    def __init__(self, cells={}):
        """
        Constructor 
        @param cells set of (i,j) tuples
        """
        # set of i,j cells 
        self.cells = cells

        # location of the centre of cells
        self.centre = []

        # ellipsis parameters
        self.ellipsisA = None
        self.ellipsisB = None
        self.ellipsisAngle = None

        # min/max indices of the box containing the set of points
        self.iMin, self.jMin, self.iMax, self.jMax = None, None, None, None
    
    def update(self):

        self.centre, self.ellipsisA, self.ellipsisB, self.angle = \
                 self.getEllipsis()

        self.iMin = numpy.minimum([c[0] for c in self.cells])
        self.jMin = numpy.minimum([c[1] for c in self.cells])
        self.iMax = numpy.maximum([c[0] for c in self.cells])
        self.jMax = numpy.minimum([c[1] for c in self.cells])

    def merge(self, otherCluster):
        self.cells.union(otherCluster.cells)
        self.update()

    def overlaps(self, otherCluster):

        # quick check...
        res = False
        if otherCluster.iMax < self.iMin or otherCluster.iMin > self.iMax:
            return res
        if otherCluster.jMax < self.jMin or otherCluster.jMin > self.jMax:
            return res

        # the clusters overlap if the ellipses overlap.
        # TO IMPLEMENT

        return True

    def toArray(self, bounds=[]):
        """
        Convert to numpy array
        @bounds (iMin, iMax, jMin, jMax)
        @return array of coordinates, array of zeros and ones
        """
        iMin, iMax, jMin, jMax = self.iMin, self.iMax, self.jMin, self.jMax
        if bounds:
            iMin, iMax, jMin, jMax = bounds
        iCoords = numpy.arange(iMin, iMax + 1)
        jCoords = numpy.arange(jMin, jMax + 1)
        ijValues = numpy.zeros((len(iCoords), len(jCoords)), numpy.int32)
        ijValues[ (c for c in self.cells) ] = 1

        return iCoords, jCoords, ijValues

    def show(self):
        """
        Plots the cluster
        """
        from matplotlib import pylab
        iCoords, jCoords, ijValues = self.toArray()
        pylab.matshow(ijValues)
        pylab.show()
         
    def write_file(self, filename):
        """
        Write to netcdf file
        @param filename file name
        """
        import netCDF4
        iCoords, jCoords, ijValues = self.toArray()
        nc = netCDF4.Dataset(filename, 'w', format="NETCDF4")
        iDim = nc.createDimension('i', size=ijCoords.shape[0])
        jDim = nc.createDimension('j', size=ijCoords.shape[1])
        var = nc.createVariable('nb', 'i4', dimensions=(iDim, jDim))
        iDim[:] = iCoords
        jDim[:] = jCoords
        var[...] = ijValues
        nc.close()


    def getNumCells(self):
        return len(self.cells)

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

# private methods 

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



#############################################################################################
def test0():
    # should be able to create a cluster with nothing in it
    cluster = Cluster()
    cluster.update()
    cluster.write_file('test0.nc')

if __name__ == '__main__':
    test0()

