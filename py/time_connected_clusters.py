from cluster import Clust
import numpy as np
import netCDF4

"""
Manages clusters across time in such a way that we one can easily extract all the clusters of a give Id and time index
"""

class TimeClusterSystem:

    def __init__(self):
        """
        No argument constructor
        """
        # flat list of clusters
        self.clusters = []

        # list of dictionaries. Each dictionary represents clusters across time which share the same 
        # ID. 
        self.cluster_connectivity = []

        # min/max values in index space
        self.i_minmax = []
        self.j_minmax = []
        self.t_max = -1


    def addTime(self, data, thresh_low, thresh_high):
        """
        Add a time from image data
        @param data image data
        @param thresh_low low precip threshold
        @param threshold_high high precip threshold
        """

        # read the data and create a list of clusters


        # iterate over each of the new clusters and determine the connection between each new cluster
        # and existing clusters stored at previous time steps.

        pass

    def reduce(self, cluster_list):
        """
        Reduce the list of clusters by merging overlapping clusters
        @param cluster_list input cluster list
        @return output cluster list
        """
        return []

    def writeFile(self, filename, i_minmax=[], j_minmax=[]):
        """
        Write data to netcdf file
        @param filename file name
        @param i_minmax min/max lat indices
        @param j_minmax min/max lon indices
        """
        f = netCDF4.Dataset(filename, 'w')

        # create dimensions

        iMin, iMax = self.i_minmax
        if i_minmax:
            iMin = min(self.i_minmax[0], i_minmax[0])
            iMax = max(self.i_minmax[1], i_minmax[1])
        iDim = f.createDimension('iDim', size=iMax-iMin)

        jMin, jMax = self.j_minmax
        if j_minmax:
            jMin = min(self.j_minmax[0], j_minmax[0])
            jMax = max(self.j_minmax[1], j_minmax[1])
        jDim = f.createDimension('jDim', size=jMax-jMin)

        # inifinte dimension
        tDim = f.createDimension('tDim', size=None)

        # create variables

        i_index = f.createVariable('i_index', 'i4', ('iDim',))
        j_index = f.createVariable('j_index', 'i4', ('jDim',))
        t_index = f.createVariable('t_index', 'i4', ('tDim',))

        # check ordering!!
        nb = f.createVariable('nb', 'i4', ('tDim', 'jDim', 'iDim'))

        # write the data
        i_index[:] = np.arange(iMin, iMax + 1)
        j_index[:] = np.arange(jMin, jMax + 1)
        t_index[:] = np.arange(0, self.t_max + 1)

        # check ordering!!
        data = np.zeros((self.t_max, jMax - jMin, iMax - iMin), np.int32)

        for cluster_id in range(len(self.cluster_connectivity)):
            for time_index in range(self.t_max):
                clusters = self.getClusters(cluster_id, time_index)
                for cl in clusters:
                    iCoords, jCoords, ijVals = cl.toArray()
                    data[time_index, iCoords, jCoords] = ijVals
        nb[:] = data



    def getNumberOfIds(self):
        """
        Get the number of independent clusters
        @return number
        """
        return len(cluster_connectivity)



    def getClusters(self, cluster_id, time_index):
        """
        Get the clusters of given ID at time index
        @param cluster_id
        @param time_index
        @return list of cluster
        """
        return [self.clusters[i] for i in self.cluster_connectivity[cluster_id][time_index]]



###############################################################################
def test(filename):
    import netCDF4
    
    f = netCDF4.Dataset(filename, 'r')
    data = np.flipud(f.variables['CMORPH'][:, lat_slice, lon_slice])
    # TO DO 
    

if __name__ == '__main__':
    filename = 'Cmorph-2010_01_10.nc'
    test(filename)
