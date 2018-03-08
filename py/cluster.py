import numpy
import math
from ellipsis import Ellipsis

class Cluster:

    def __init__(self, cells={}):
        """
        Constructor 
        @param cells set of (i,j) tuples
        """
        # set of i,j cells 
        self.cells = cells

        # ellipsis representing the "average" distribution
        # of cells
        self.ellpsis = None

        # min/max indices of the box containing the set of points
        self.box = [[None, None], [None, None]]
    

    def getNumCells(self):
        """
        Get the number of cells
        @return number
        """
        return len(self.cells)


    def update(self):

        if len(self.cells) > 0:

            self.ellipsis = Ellipsis(self.cells)

            for dim in range(0, 2):
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
        for dim in range(0, 2):
            # max is smaller than this min
            noOverlap &= otherCluster.box[1][dim] < self.box[0][dim]
            noOverlap &= self.box[1][dim] < otherCluster.box[0][dim]
            # min is bigger than other max
            noOverlap &= otherCluster.box[0][dim] > self.box[1][dim]
            noOverlap &= self.box[0][dim] > otherCluster.box[1][dim]

        if noOverlap:
            # the boxes are disjoint, no chance that there is overlap
            return False

        # the clusters one centre is inside the other ellipse
        if self.ellipsis.isPointInside(otherCluster.ellipsis.getCentre()):
            return True

        if otherCluster.ellipsis.isPointInside(self.ellipsis.getCentre()):
            return True

        return False

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
        for c in self.cells:
            ijValues[c[0] - iMin, c[1] - jMin] = 1

        return iCoords, jCoords, ijValues

    def __repr__(self):
        """
        Print object
        """
        res = """
        Cluster: num cells = {} box = {} ellipsis centre = {} a = {} b = {} transf = {}
        """.format(len(self.cells), self.box, \
            self.ellipsis.centre, self.ellipsis.a, self.ellipsis.b, \
            self.ellipsis.ij2AxesTransf)
        return res


#############################################################################################
def test0():
    # should be able to create a cluster with nothing in it
    cluster = Cluster()
    cluster.update()
    # not sure why writing a zero dimensioned variable is not working with netcdf4
    #cluster.writeFile('test0.nc')

def test1():
    cluster = Cluster({(-1, -2)})
    cluster.update()
    print('test1 {}'.format(cluster))
    #cluster.writeFile('test1.nc')

def testHorizLine():
    cluster = Cluster({(-1, -2), (0, -2), (1, -2), (2, -2)})
    cluster.update()
    print('testHorizLine {}'.format(cluster))

def testDipole():
    cluster = Cluster({(-2, 0), (2, 0)})
    cluster.update()
    print('testDipole {}'.format(cluster))

def testRectangle():
    cluster = Cluster({(i, 0) for i in range(3)}.union({(i, 1) for i in range(3)}))
    cluster.update()
    print('testRectangle {}'.format(cluster))

def testRectangleSlanted():
    cluster = Cluster({(i, 0) for i in range(4)}.union({(i - 1, 1) for i in range(4)}))
    cluster.update()
    print('testRectangleSlanted {}'.format(cluster))

def testRandom():
    import random
    random.seed(1234)
    cluster = Cluster({(random.randint(0, 200), random.randint(0, 100)) for i in range(500)})
    cluster.update()
    print('testRandom {}'.format(cluster))


if __name__ == '__main__':
    test0()
    test1()
    testHorizLine()
    testDipole()
    testRectangle()
    testRectangleSlanted()
    testRandom()

