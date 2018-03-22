import numpy
import matplotlib
import time
from time_connected_clusters import TimeConnectedClusters, extractClusters
from cluster import Cluster

"""
Test TimeConnectedClusters
"""

def testRectangle():
    rectangle = {(2, 3), (3, 3), (4, 3), (2, 4), (3, 4), (4, 4)}
    clusters = [Cluster(rectangle)]

    tcc = TimeConnectedClusters()
    tcc.addTime(clusters)
    print tcc
    tcc.addTime(clusters)
    print tcc


def testTwoRectangles():
    rectangle = {(2, 3), (3, 3), (4, 3), (2, 4), (3, 4), (4, 4)}
    clusters = [Cluster(rectangle)]

    tcc = TimeConnectedClusters()
    tcc.addTime(clusters)
    print tcc


def testDigits():
    """
    Checking that we can create a time clonnected cluster with no cluster
    """
    lats = numpy.arange(0, 20)
    lons = numpy.arange(0, 50)

    data = numpy.zeros((len(lats), len(lons)), numpy.float64)

    # create a cluster
    data[10:10+5, 10] = 1

    """
    data[2, 5:5+4] = 2
    data[4, 5:5+4] = 2
    data[6, 5:5+4] = 2
    data[3, 5] = 2
    data[5, 8] = 2
    """

    print data
    clusters = extractClusters(numpy.flipud(data), thresh_min=0., thresh_max=10.)
    print clusters

    tcc = TimeConnectedClusters()

    # add a couple of times
    tcc.addTime(clusters)
    #tcc.addTime(clusters)

    print tcc

    # write to file
    #tcc.writeFile('digits.nc')




if __name__ == '__main__':
    testRectangle()
    #testDigits()
