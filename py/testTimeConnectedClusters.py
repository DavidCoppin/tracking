import numpy
import matplotlib
import time
from time_connected_clusters import TimeConnectedClusters, extractClusters

"""
Test TimeConnectedClusters
"""

def testDavid():
    """
    Checking that we can create a time clonnected cluster with no cluster
    """
    lats = numpy.arange(0, 20)
    lons = numpy.arange(0, 30)

    data = numpy.zeros((len(lats), len(lons)), numpy.float32)

    # create a cluster
    data[4:4+5, 2] = 1
    data[5:5+3, 4] = 1
    data[4, 3] = 1
    data[8, 3] = 1

    for i in range(5):
        data[4+i, 6+i] = 2
        data[8-i, 9+i] = 2
    data[5, 7:14] = 2

    clusters = extractClusters(numpy.flipud(data), thresh_min=0.0, thresh_max=0.8)

    tcc = TimeConnectedClusters()
    # add a couple of times
    tcc.addTime(clusters)
    tcc.addTime(clusters)

    print tcc

    # write to file
    #tcc.writeFile('david.nc')




if __name__ == '__main__':
    testDavid()
