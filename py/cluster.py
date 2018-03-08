import numpy
import math

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
        self.box = [[None, None], [None, None]]
    
    def update(self):

        if len(self.cells) > 0:

            self.centre, self.ellipsisA, self.ellipsisB, self.ellipsisAngle = \
                 self.getEllipsis()

            for dim in range(0, 1):
                self.box[0][dim] = numpy.min([c[dim] for c in self.cells])
                self.box[1][dim] = numpy.max([c[dim] for c in self.cells])


    def merge(self, otherCluster):
        """
        Merge this cluster with another
        @param otherCluster
        @return new cluster that is the union of this and otherCluster
        """
        self.cells.union(otherCluster.cells)
        self.update()

    def overlaps(self, otherCluster):
        """
        Find our if this clsuter overlaps with otherCluster
        @param otherCluster
        @return True if there is overlap, False otherwise
        """
        # quick check if the boxes don't overlap...
        noOverlap = True
        for dim in range(0, 1):
            # max is smaller than this min
            noOverlap &= otherCluster.box[1][dim] < self.box[0][dim]
            noOverlap &= self.box[1][dim] < otherCluster.box[0][dim]
            # min is bigger than other max
            noOverlap &= otherCluster.box[0][dim] > self.box[1][dim]
            noOverlap &= self.box[0][dim] > otherCluster.box[1][dim]

        if noOverlap:
            return False

        # the clusters one centre is inside the other ellipse
        if self.__isPointInsideEllipse(self.centre, 
                otherCluster.centre, 
                otherCluster.ellipsisA, otherCluster.ellipsisB,
                otherCluster.ellipsisAngle):
            return True

        if self.__isPointInsideEllipse(otherCluster.centre, 
                self.centre, 
                self.ellipsisA, self.ellipsisB,
                self.ellipsisAngle):
            return True

        return False

    def toArray(self, bounds=[]):
        """
        Convert this cluster to numpy (dense) array
        @bounds [[iMin, jMin], [iMax, jMax]]
        @return array of coordinates, array of zeros and ones
        """
        if len(self.cells) <= 0:
            # no op
            return numpy.array([]), numpy.array([]), numpy.array([])

        if not bounds:
            bounds = self.box
        iCoords = numpy.arange(bounds[0][0], bounds[1][0] + 1)
        jCoords = numpy.arange(bounds[0][1], bounds[1][1] + 1)
        ijValues = numpy.zeros((len(iCoords), len(jCoords)), numpy.int32)
        iMin, jMin = bounds[0]
        print [c for c in self.cells]
        print [(c[0] - iMin, c[1] - jMin) for c in self.cells]
        ijValues[ [(c[0] - iMin, c[1] - jMin) for c in self.cells] ] = 1

        return iCoords, jCoords, ijValues

    def getEllipseAsPolyline(self, numSegments=32):
        """
        Return the ellipsis as a segmented line
        @return iPts, jPts arrays
        """
        # rotation of the i j coords to the principal axes of the inertia matrix
        cosa = math.cos(self.ellipsisAngle)
        sina = math.sin(self.ellipsisAngle)
        transf = numpy.array([[cosa, sina], [-sina, cosa]])
        invTrans = numpy.linalg.inv(transf)

        iPts, jPts = [], []
        dt = 2 * math.pi / float(numSegments)
        for i in range(numSegments + 1):
            th = i * dt
            x = self.ellipsisA * math.cos(th)
            y = self.ellipsisB * math.cos(th)
            # rotate back to i,j coordinates
            ij = invTrans.dot([x, y])
            ij + self.centre
            iPts.append(ij[0])
            jPts.append(ij[1])

        return iPts, jPts


    def show(self):
        """
        Plots the cluster
        """
        from matplotlib import pylab
        bounds = [[self.box[0][1] - 1, self.box[0][1] - 1], 
                  [self.box[1][0] + 1, self.box[1][1] + 1]]
        iCoords, jCoords, ijValues = self.toArray(bounds=bounds)
        # show the cluster
        pylab.matshow(ijValues)
        # show the ellipsis
        iPts, jPts = self.getEllipseAsPolyline()
        pylab.plot(iPts, jPts, 'm-')
        pylab.show()

         
    def writeFile(self, filename):
        """
        Write to netcdf file
        @param filename file name
        """
        import netCDF4
        iCoords, jCoords, ijValues = self.toArray()
        nc = netCDF4.Dataset(filename, 'w', format="NETCDF4")
        iDim = nc.createDimension('iDim', size=iCoords.shape[0])
        jDim = nc.createDimension('jDim', size=jCoords.shape[0])
        iVar = nc.createVariable('i', 'i4', dimensions=('iDim',))
        jVar = nc.createVariable('j', 'i4', dimensions=('jDim',))
        nbVar = nc.createVariable('nb', 'i4', dimensions=('iDim', 'jDim'))
        iVar[:] = iCoords
        jVar[:] = jCoords
        nbVar[:, :] = ijValues
        nc.close()


    def getNumCells(self):
        """
        Get the number of cells
        @return number
        """
        return len(self.cells)

    def getEllipsis(self):
        """
        Get the ellipsis parameters that best represent this cluster
        @return centre, radius1, radius2, angle
        """
        centre = self.__getCentre()
        inertia = numpy.zeros((2, 2), numpy.float64)
        inertia[0, 0] = numpy.sum([(c[0] - centre[0])**2 for c in self.cells])
        inertia[0, 1] = numpy.sum([(c[0] - centre[0])*(c[1] - centre[1]) for c in self.cells])
        inertia[1, 0] = inertia[0, 1]
        inertia[1, 1] = numpy.sum([(c[1] - centre[1])**2 for c in self.cells])

        # diagonalize
        eigenvals, eigenvecs = numpy.linalg.eig(inertia)
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

    def __isPointInsideEllipse(self, point, 
        ellipsisCentre, ellipsisA, ellipsisB, ellipsisAngle):
        # TO DO 
        pass



#############################################################################################
def test0():
    # should be able to create a cluster with nothing in it
    cluster = Cluster()
    cluster.update()
    # not sure why writing a zero dimensioned variable is not working with netcdf4
    #cluster.writeFile('test0.nc')

def test1():
    # should be able to create a cluster with nothing in it
    cluster = Cluster({(-1, -2)})
    cluster.update()
    #cluster.show()
    cluster.writeFile('test1.nc')

def testHorizLine():
    # should be able to create a cluster with nothing in it
    cluster = Cluster({(-1, -2), (0, -2), (1, -2), (2, -2)})
    cluster.update()
    cluster.show()


if __name__ == '__main__':
    test0()
    test1()
    testHorizLine()

