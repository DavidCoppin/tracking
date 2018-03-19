from cluster import Clust
import numpy as np
import netCDF4

"""
Manages clusters across time in such a way that we one can easily extract all the clusters of a give Id and time index
"""

class TimeConnectedClusters:

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
        new_clusters = self.exactClustersFromData(data, thresh_low, thresh_high)
        new_clusters = self.reduce(new_clusters)

        # special case if initial time
        if self.t_max < 0:
            self.clusters += new_clusters
            self.t_max += 1
            for cluster_id in range(len(new_clusters)):
                self.cluster_connectivity[cluster_id] = {self.t_max: [cluster_id]}
            # done
            return

        # add the new clusters
        start_id = len(self.clusters)
        self.clusters += new_clusters
        self.t_max += 1

        num_new_clusters = len(new_clusters)
        num_clusters = len(self.cluster_connectivity)

        #
        # establish the connectivity with past clusters
        #

        # track forward
        previous_time_cluster_ids = [self.cluster_connectivity[i].get(self.t_max - 1, []) \
                                        for i in range(num_clusters)]
        for i in range(num_new_clusters):
            cli = new_clusters[i]
            for j in previous_time_cluster_ids:
                clj = self.clusters[j]
                if cli.isCentreInsideOf(clj):
                    # add this cluster
                    self.cluster_connectivity[j][self.t_max] = self.cluster_connectivity[j].get(self.t_max, []) + [start_id + i]

        # track backward
        for i in range(num_new_clusters):
            cli = new_clusters[i]
            for t_index in range(self.t_max - 1, -1 , -1):
                for j in range(len(self.cluster_connectivity):
                    # TO DO 





        # backward tracking

        # iterate over each of the new clusters and determine the connection between ecah cluster
        # and existing clusters stored at previous time steps.

        num_ids = len(self.cluster_connectivity)

        current_time_index = self.t_max
        for cluster_id in range(num_ids):

            # forward tracking 
            current_clusters = [self.clusters[i] for i in range(self.cluster_connectivity[cluster_id].get(current_time_index, []))]




    def extractClusters(self, data, thresh_low, thresh_high):
        """
        """
        ma_data = np.ma.masked_where(data <= thresh_min, data).astype(np.uint8)
        
        # building threshold
        tmp_data = ma_data.filled(0)
        tmp_data[np.where(tmp_data!=0)] = 255
        bw_data = tmp_data.astype(np.uint8)
        
        # building markers and borders
        ma_conv = np.ma.masked_where(data <= thresh_max, data)
        tmp_conv = ma_conv.filled(0)
        tmp_conv[np.where(tmp_conv!=0)] = 255
        bw_conv = tmp_conv.astype(np.uint8)
        markers= ndimage.label(bw_conv, structure=np.ones((3, 3)))[0]
        border = cv2.dilate(bw_data, None, iterations=5)
        border -= cv2.erode(border, None)
        markers[border == 255] = 255

        # watershed
        ws = watershed(-data, markers, mask=bw_data)

        # load into clusters
        res = []
        for idVal in range(1, ws.max()):
            iVals, jVals = np.where(ws == idVal)
            numVals = len(iVals)
            if numVals > 0:
                cells = [(iVals[i], jVals[i]) for i in range(len(iVals))]
                # store this cluster as a list with one element (so far). Each 
                # element will have its own ID
                res.append([Cluster(cells)])

        return res



    def reduce(self, cluster_list):
        """
        Reduce the list of clusters by merging overlapping clusters
        @param cluster_list input cluster list
        @return output cluster list
        """
        res = copy.deepcopy(cluster_list)
        n = len(cluster_list)
        remove_indices = []

        for i in range(n):
            cli = res[i]
            for j in range(i + 1, n):
                clj = cluster_list[j]
                if cli.isCentreInsideOf(clj) and clj.isCentreInsideOf(cli):
                    # merge and tag clj for removal
                    cli += clj
                    remove_indices.append(j)

        remove_indices.revert()
        for i in remove_indices:
            del res[i]

        return res

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
