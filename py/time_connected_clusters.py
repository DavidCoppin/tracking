from cluster import Cluster
import numpy as np
import netCDF4
import copy

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
        self.cluster_connect = []

        # min/max values in index space
        self.i_minmax = []
        self.j_minmax = []
        self.t_index = 0


    def fuse(self, id1, id2):
        """
        Fuse time clusters, the result will be stored in id1 and id2 will be removed entirely
        @param id1 Id of track 1, track 1 will contain track 2
        @param id2 Id of track 2, track 2 will be removed
        """
        for t_index in range(self.t_index):
            clusters = self.cluster_connect[id1].get(t_indx, []) + \
                       self.cluster_connect[id2].get(t_indx, [])
            if clusters:
                self.cluster_connect[id1][t_index] = clusters

        del self.cluster_connect[id2]


    def copy(self, id1, t_index, index_list, id2):
        """
        Copy the clusters from one Id to another for a given time index
        @param id1 the Id source
        @param t_index time index
        @param index_list list of indices of cluster to be copied
        @param id2 destination of the clusters to be moved
        """
        self.cluster_connect[id2][t_index] = self.cluster_connect[id2].get(t_index, []) \
           + [self.cluster_connect[id1][t_index][i] for i in index_list]


    def move(self, id1, t_index, index_list, id2):
        """
        Move the clusters form one Id to another for a given time index, the clusters
          will be deleted from the source track
        @param id1 the Id source track
        @param t_index time index
        @param index_list list of indices of cluster to be copied
        @param id2 destination track Id where the clusters will be moved        
        """
        self.copy(id1, t_index, cluster_index_list, id2)
        # remove the clusters form Id1, starting form the last 
        cluster_index_list.sort()
        cluster_index_list.reverse()
        for i in cluster_index_list:
            del self.cluster_connect[id1][t_index][i]


    def addTime(self, new_clusters):
        """
        Add time entry
        @param new_clusters list of new clusters
        """

        # merge overlapping clusters
        new_clusters = self.reduce(new_clusters)

        index = len(self.clusters)

        # special case if first time step
        if self.t_index == 0:
            # there is no past
            # just create a new track for each cluster
            for new_cl in new_clusters:

                self.clusters.append(new_cl)
                self.cluster_connect.append({self.t_index: [index]})
                index += 1

            # done
            return 

        # assign the new clusters to existing tracks
        for new_cl in new_clusters:
            self.clusters.append(new_cl)
            # find out if this cluster belongs to an existing track
            connected_clusters = []
            track_ids = []
            for track_id in range(len(self.cluster_connect)):
                old_cluster_inds = self.cluster_connect[track_id].get(self.t_index - 1, [])
                for old_cl in [self.clusters[i] for i in old_cluster_inds]:
                    # are new_cl and old_cl overlapping?
                    if old_cl.isCentreInsideOf(new_cl): # and new_cl.isCentreInsideOf(old_cl):
                        connected_clusters.append(old_cl)
                        track_ids.append(track_id)
            # this cluster could be assigned to any of the tracks in track_ids
            if len(track_ids) == 0:
  	        self.cluster_connect.append({self.t_index: [index]})
            else:
                # choose the track for which the distance between new_cl is smallest to 
                # to any of the clusters of that track at t - dt
                dists = np.array([new_cl.getDistance(cl) for cl in connected_clusters])
                i = np.argmin(dists)
                self.cluster_connect[track_ids[i]][self.t_index] = \
                  self.cluster_connect[track_ids[i]].get(self.t_index, []) + [index] 

        # test if new_cl.isCentreInsideOf(old_cl): if so then fuse 


    def extractClusters(self, data, thresh_low, thresh_high):
        """
        Extract clusters from an image data
        @param data
        @param thresh_low
        @param thresh_high
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

        remove_indices.reverse()
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
        t_index[:] = np.arange(0, self.t_index + 1)

        # check ordering!!
        data = np.zeros((self.t_index, jMax - jMin, iMax - iMin), np.int32)

        for cluster_id in range(len(self.cluster_connect)):
            for time_index in range(self.t_index):
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
        return len(cluster_connect)



    def getClusters(self, track_id, time_index):
        """
        Get the clusters of given ID at time index
        @param track_id
        @param time_index
        @return list of clusters
        """
        return [self.clusters[i] for i in self.cluster_connect[track_id].get(time_index, [])]


    def __repr__(self):
        res = """
TimeConnectedCluster: num of clusters   {}
                      num of time steps {}
                      num of Ids        {}
                      clusters          {}
                      connectivity      {}
        """.format(len(self.clusters), self.t_index + 1, \
              len(self.cluster_connect), [cl.cells for cl in self.clusters], \
              self.cluster_connect)
        return res



###############################################################################
def testNoCluster():
    tcc = TimeConnectedClusters()
    print 'No cluster'
    print tcc

def testOneCluster():
    tcc = TimeConnectedClusters()
    c0 = Cluster({(1, 1), (2, 1), (2, 2)})
    tcc.addTime([c0])
    print 'One cluster'
    print tcc

def testReduceNonOverlapping():
    tcc = TimeConnectedClusters()
    c0 = Cluster({(1, 1), (2, 1), (2, 2)})
    c1 = Cluster({(1, 3), (1, 4), (2, 4)})
    reduced_list = tcc.reduce([c0, c1])
    print 'input list non-overlapping', [c0, c1]
    print 'output list non-overlapping', reduced_list

def testReduceOverlapping():
    tcc = TimeConnectedClusters()
    c0 = Cluster({(1, 1), (2, 1), (2, 2)})
    c1 = Cluster({(1, 1), (1, 2), (2, 2)})
    reduced_list = tcc.reduce([c0, c1])
    print 'input list overlapping ', [c0, c1]
    print 'output list overlapping ', reduced_list

def testTwoClustersAtTime0():
    tcc = TimeConnectedClusters()
    c0 = Cluster({(1, 1), (2, 1), (2, 2)})
    c1 = Cluster({(1, 1), (1, 2), (2, 2)})
    tcc.addTime([c0, c1])
    print 'Two clusters'
    print tcc    

if __name__ == '__main__':
    testNoCluster()
    testOneCluster()
    testReduceNonOverlapping()
    testReduceOverlapping()
    #testTwoClustersAtTime0()
