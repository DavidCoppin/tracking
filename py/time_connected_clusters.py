from cluster import Cluster
import numpy as np
import netCDF4
import copy
from scipy import ndimage
import cv2
from skimage.morphology import watershed
from skimage.feature import peak_local_max
    

def reduce(cluster_list):
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
                # add clj to cli and tag clj for removal
                cli += clj
                remove_indices.append(j)

    # remove the tagged clusters
    remove_indices.reverse()
    for i in remove_indices:
        del res[i]

    return res


def extractClusters(data, thresh_min, thresh_max):
    """
    Extract clusters from an image data
    @param data
    @param thresh_min
    @param thresh_max
    @return list of clusters
    """

    """
    image = np.logical_and(data >= thresh_min, data <= thresh_max)
    print image
    # Now we want to separate the two objects in image
    # Generate the markers as local maxima of the distance to the background
    distance = ndimage.distance_transform_edt(image)
    local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3, 3)),
                            labels=image)
    markers = ndimage.label(local_maxi)[0]
    labels = watershed(-distance, markers, mask=image)
    """


    ma_data = np.ma.masked_where(data <= thresh_min, data).astype(np.uint8)
        
    # building threshold
    tmp_data = ma_data.filled(fill_value=0)
    tmp_data[np.where(tmp_data !=0)] = 255
    bw_data = tmp_data.astype(np.uint8)
        
    # building markers and borders
    ma_conv = np.ma.masked_where(data <= thresh_max, data)
    tmp_conv = ma_conv.filled(fill_value=0)
    tmp_conv[np.where(tmp_conv !=0)] = 255
    bw_conv = tmp_conv.astype(np.uint8)
    markers = ndimage.label(bw_conv, structure=np.ones((3, 3)))[0]
    border = cv2.dilate(bw_data, None, iterations=5)
    border -= cv2.erode(border, None)
    markers[border == 255] = 255

    # labels each feature
    labels = watershed(-data, markers, mask=bw_data)

    # load into clusters
    res = []
    for idVal in range(1, labels.max()):
        iVals, jVals = np.where(labels == idVal)
        numVals = len(iVals)
        if numVals > 0:
            cells = [(iVals[i], jVals[i]) for i in range(len(iVals))]
            # store this cluster as a list with one element (so far). Each 
            # element will have its own ID
            res.append(Cluster(cells))

    return res


"""
Manages clusters across time in such a way that we one can easily extract all the clusters of a give Id and time index
"""

class TimeConnectedClusters:

    # used to initialize min/max index values
    LARGE_INT = 999999999

    def __init__(self):
        """
        No argument constructor
        """
        # flat list of clusters
        self.clusters = []

        # list of dictionaries. Each dictionary represents clusters across time which share the same 
        # ID. 
        self.cluster_connect = []

        # current time index
        self.t_index = 0


    def fuse(self, id1, id2):
        """
        Fuse time clusters, the result will be stored in id1 and id2 will be removed entirely
        @param id1 Id of track 1, track 1 will contain track 2
        @param id2 Id of track 2, track 2 will be removed
        """
        print '*** fusing ', id1, id2
        for t_index in range(self.t_index):
            clusters = self.cluster_connect[id1].get(t_index, []) + \
                       self.cluster_connect[id2].get(t_index, [])
            if clusters:
                self.cluster_connect[id1][t_index] = clusters

        del self.cluster_connect[id2]


    def addTime(self, new_clusters):
        """
        Add time entry
        @param new_clusters list of new clusters
        """

        # merge overlapping clusters
        new_clusters = reduce(new_clusters)

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
            self.t_index += 1
            return 

        # assign the new clusters to existing tracks

        for new_cl in new_clusters:

            new_cl_track_id = -1

            # add the cluster
            self.clusters.append(new_cl)

            # find out if this cluster belongs to an existing track. This is equivalent to 
            # forward tracking
            num_tracks = self.getNumberOfTracks()
            connected_clusters = []
            connected_track_ids = []
            for track_id in range(num_tracks):
                old_cluster_inds = self.cluster_connect[track_id].get(self.t_index - 1, [])
                for old_cl in [self.clusters[i] for i in old_cluster_inds]:
                    # is the centre of old_cl inside the ellipse of new_cl?
                    if old_cl.isCentreInsideOf(new_cl):
                        connected_clusters.append(old_cl)
                        connected_track_ids.append(track_id)

            # this cluster could be assigned to any of the tracks in track_ids
            if len(connected_track_ids) == 0:
                # this cluster is on its own
                # create a new entry 
                self.cluster_connect.append({self.t_index: [index]})
                new_cl_track_id = self.getNumberOfTracks()
            else:
                # choose the track for which the distance between new_cl is smallest to 
                # to any of the clusters of that track at t - dt
                dists = np.array([new_cl.getDistance(cl) for cl in connected_clusters])
                i = np.argmin(dists)
                new_cl_track_id = connected_track_ids[i]
                self.cluster_connect[new_cl_track_id][self.t_index] = self.cluster_connect[new_cl_track_id].get(self.t_index, []) + [index]

            # backward tracking
            track_ids_to_fuse = []
            for track_id in range(num_tracks):

                if track_id == new_cl_track_id:
                    # skip self
                    continue

                old_cluster_inds = self.cluster_connect[track_id].get(self.t_index - 1, [])
                for old_cl in [self.clusters[i] for i in old_cluster_inds]:

                    # is the centre of new_cl inside the ellipse of old_cl?
                    if new_cl.isCentreInsideOf(old_cl):
                        # the two tracks take the same Id 
                        self.fuse(track_id, new_cl_track_id)

            # update the cluster index
            index += 1

        # update the time index
        self.t_index += 1


    def getMinMaxIndices(self):
        """
        Get the low/high end box indices
        @return i_min, j_min, i_max, j_max
        """
        i_min, j_min, i_max, j_max = self.LARGE_INT, self.LARGE_INT, -self.LARGE_INT, -self.LARGE_INT
        for track in self.cluster_connect:
            for cl_indx_list in track.values():
                i_min = min(i_min, min([self.clusters[k].box[0][0] for k in cl_indx_list]))
                j_min = min(j_min, min([self.clusters[k].box[0][1] for k in cl_indx_list]))
                i_max = max(i_max, max([self.clusters[k].box[1][0] for k in cl_indx_list]))
                j_max = max(j_max, max([self.clusters[k].box[1][1] for k in cl_indx_list]))
        return i_min, j_min, i_max, j_max



    def writeFile(self, filename, i_minmax=[], j_minmax=[]):
        """
        Write data to netcdf file
        @param filename file name
        @param i_minmax min/max lat indices
        @param j_minmax min/max lon indices
        """
        f = netCDF4.Dataset(filename, 'w')

        # create dimensions

        iMin, jMin, iMax, jMax = self.getMinMaxIndices()

        if i_minmax:
            iMin = min(iMin, i_minmax[0])
            iMax = max(iMax, i_minmax[1])
        num_i = iMax - iMin + 1
        iDim = f.createDimension('iDim', size=num_i)

        if j_minmax:
            jMin = min(jMin, j_minmax[0])
            jMax = max(jMax, j_minmax[1])
        num_j = jMax - jMin + 1
        jDim = f.createDimension('jDim', size=num_j)

        # inifinte dimension
        tDim = f.createDimension('tDim', size=None)

        # create variables

        i_index = f.createVariable('i_index', 'i4', ('iDim',))
        j_index = f.createVariable('j_index', 'i4', ('jDim',))
        t_index = f.createVariable('t_index', 'i4', ('tDim',))

        # check ordering!!
        nb_var = f.createVariable('nb', 'i4', ('tDim', 'jDim', 'iDim'))

        # write the data
        i_index[:] = np.arange(iMin, iMax + 1)
        j_index[:] = np.arange(jMin, jMax + 1)
        t_index[:] = np.arange(0, self.t_index)

        # data buffer, check ordering!!
        data = np.zeros((self.t_index, num_j, num_i), np.int32)

        for track_id in range(self.getNumberOfTracks()):
            for t_index in range(self.t_index):
                clusters = self.getClusters(track_id, t_index)
                for cl in clusters:
                    n_cells = cl.getNumCells()
                    tis = [t_index] * n_cells
                    jis = [c[1] - jMin for c in cl.cells]
                    iis = [c[0] - iMin for c in cl.cells]
                    # check ordering!!
                    data[tis, jis, iis] = track_id + 1
        # now write all the data in one go
        nb_var[:] = data

        f.close()
        

    def getNumberOfTracks(self):
        """
        Get the number of tracks
        @return number
        """
        return len(self.cluster_connect)


    def getLastTimeStep(self):
        """
        Get the last time step (also number of time steps)
        @return number
        """
        return self.t_index



    def getClusters(self, track_id, time_index):
        """
        Get the clusters of given ID at time index
        @param track_id track Id
        @param time_index
        @return list of clusters
        """
        return [self.clusters[i] for i in self.cluster_connect[track_id].get(time_index, [])]


    def __repr__(self):
        res = """
TimeConnectedCluster: num of clusters   {}
                      num of time steps {}
                      num of tracks     {}
                      clusters          {}
                      connectivity      {}
        """.format(len(self.clusters), self.t_index, \
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
    c0 = Cluster({(1, 1), (2, 1), (2, 2)})
    c1 = Cluster({(1, 3), (1, 4), (2, 4)})
    reduced_list = reduce([c0, c1])
    print 'input list non-overlapping', [c0, c1]
    print 'output list non-overlapping', reduced_list

def testReduceOverlapping():
    c0 = Cluster({(1, 1), (2, 1), (2, 2)})
    c1 = Cluster({(1, 1), (1, 2), (2, 2)})
    reduced_list = reduce([c0, c1])
    print 'input list overlapping ', [c0, c1]
    print 'output list overlapping ', reduced_list

def testTwoClustersAtTime0():
    tcc = TimeConnectedClusters()
    c0 = Cluster({(1, 1), (2, 1), (2, 2)})
    c1 = Cluster({(1, 1), (1, 2), (2, 2)})
    tcc.addTime([c0, c1])
    print tcc
    c3 = Cluster({(4,3), (7,3)})
    tcc.addTime([c3])
    print tcc
    tcc.writeFile('twoClusters.nc')
    print 'Two clusters'
    print tcc    

if __name__ == '__main__':
    testNoCluster()
    testOneCluster()
    testReduceNonOverlapping()
    testReduceOverlapping()
    testTwoClustersAtTime0()
