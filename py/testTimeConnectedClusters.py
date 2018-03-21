import numpy
import matplotlib
import time
from time_connected_clusters import TimeConnectedClusters, extractClusters

"""
Test TimeConnectedClusters
"""

def testNoCluster(lat_indx_slice, lon_indx_slice):
    """
    Checking that we can create a time clonnected cluster with no cluster
    """
    lats = numpy.arange(lat_indx_slice.start, lat_indx_slice.stop)
    lons = numpy.arange(lon_indx_slice.start, lon_indx_slice.stop)
    data = numpy.zeros((len(lats), len(lons)), numpy.int32)
    clusters = extractClusters(numpy.flipud(data), thresh_min=0.0, thresh_max=0.8)

    tcc = TimeConnectedClusters()
    # add a couple of times
    tcc.addTime(clusters)
    tcc.addTime(clusters)

    # write to file
    tcc.writeFile('no_cluster.nc')




if __name__ == '__main__':
    testNoCluster(lat_indx_slice=slice(200, 300), lon_indx_slice=slice(50, 350))
