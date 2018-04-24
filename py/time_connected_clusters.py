from cluster import Cluster
from ellipse import Ellipse
import numpy as np
import netCDF4
import copy
import math
import cPickle
import functools
    

def __reduceOne(cluster_list, frac):
    """
    Reduce the list of clusters by merging overlapping clusters
    @param cluster_list in/out cluster list
    @param frac: threshold for overlapping ellipses to be reduced
    @return True if the list was reduced
    """
    n = len(cluster_list)
    for i in range(n):
        cli = cluster_list[i]
        for j in range(i + 1, n):
            clj = cluster_list[j]

            if cli.isCentreInsideOf(clj) and clj.isCentreInsideOf(cli):
                cli += clj
                del cluster_list[j]
                return True

    return False


def reduce(cluster_list, frac=1.0):
    """
    Fully reduce until no more reduction can be applied 
    @param cluster_list in/out cluster list
    @param frac: threshold for overlapping ellipses to be reduced
    """
    go = True
    while go:
        go = __reduceOne(cluster_list, frac)



"""
Manages clusters across time in such a way that we one can easily extract all the clusters 
of a given Id and time index
"""

class TimeConnectedClusters:

    # used to initialize min/max index values
    LARGE_INT = 999999999

    def __init__(self):
        """
        Constructor
        """
        # flat list of clusters
        self.clusters = []

        # precipitation data
        self.precip = []

        # list of dictionaries. Each dictionary represents clusters across time which share the same 
        # ID. 
        self.cluster_connect = []

        # current time index
        self.t_index = 0

        # time
        self.time = []


    def getPrecip(self, var, time):
        """
        Store precipitation data and append it
        """
        if np.shape(self.precip)[0]==0:
            self.precip = var
            self.time = time
        else:
            self.precip = np.vstack([self.precip, var])
            self.time = np.append(self.time, time)

    def fuse(self, track_ids):
        """
        Fuse track Ids. Only the first track ID will survive, all other track Ids
        will be folded into the first one
        @param track_ids set of track Ids
        @return track Ids marked for deletion
        """

        tr_ids = list(track_ids)
        tr_ids.sort()
        track_id0 = tr_ids[0]
        cl_conn0 = self.cluster_connect[track_id0]

        for track_id in tr_ids[1:]:
            for t_index, cl_list in self.cluster_connect[track_id].items():
                cl_conn0[t_index] = cl_conn0.get(t_index, []) + cl_list

        # return the track Ids marked for deletion
        n = len(tr_ids)
        return [tr_ids[i] for i in range(n - 1, 0, -1)]


    def addTime(self, new_clusters, frac):
        """
        Add time entry
        @param new_clusters list of new clusters
        @param frac TO DESCRIBE
        """
        # merge overlapping clusters, this will reduce the number of 
        # clusters to track but should have no influence on the 
        # final result
        reduce(new_clusters, frac)

        # current number of clusters
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


        new_track_ids = self._forwardTracking(new_clusters)

        new_track_ids_to_fuse = self._backwardTracking(new_track_ids)

        # reduce list of tracks to fuse by looking at common track Ids between
        # the elements of new_track_ids_to_fuse. These elements will be tagged for 
        # removal in delete_elem
        self._fuseAll(new_track_ids_to_fuse)
 
        # done with assigning, update the time index
        self.t_index += 1


    def _forwardTracking(self, new_clusters):
        """
        Forward tracking: assign new clusters to existing tracks
        @param new_clusters new clusters
        @return set of track Ids to which the new clusters belong
        """

        # current number of clusters
        index = len(self.clusters)

        # set of track Ids to which the new clusters will be assigned to
        new_track_ids = set() # new cluster index to track id

        for new_cl_index in range(len(new_clusters)):

            new_cl = new_clusters[new_cl_index]

            # the tack Id that we need ot assign this cluster to
            new_track_id = -1

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
                    # is the centre of new_cl inside the ellipse of old_cl?
                    if new_cl.isCentreInsideOfExt(old_cl) or old_cl.isCentreInsideOfExt(new_cl):
                        connected_clusters.append(old_cl)
                        connected_track_ids.append(track_id)

            # this cluster could be assigned to any of the tracks in track_ids
            if len(connected_track_ids) == 0:
                # this cluster is on its own
                # create a new entry 
                new_track_id = self.getNumberOfTracks()
                self.cluster_connect.append({self.t_index: [index]})
            else:
                # choose the track for which the distance between new_cl is smallest to 
                # to any of the clusters of that track at t - dt
                dists = np.array([new_cl.getDistance(cl) for cl in connected_clusters])
                i = np.argmin(dists)
                new_track_id = connected_track_ids[i]
                self.cluster_connect[new_track_id][self.t_index] = \
                               self.cluster_connect[new_track_id].get(self.t_index, []) + [index]

            new_track_ids.add(new_track_id)

            # update self.cluster's index
            index += 1

        return new_track_ids


    def getBigClusterAt(self, track_id, t_index):
        """
        Construct a big cluster from all the track_id clusters at time t_index
        @param track_id track Id
        @param t_index time index
        @return one big cluster representing the merge of all smaller clusters
        """
        clusters = [self.clusters[i] for i in self.cluster_connect[track_id].get(t_index, [])]
        if not clusters:
            return None
        all_cells = functools.reduce(lambda x, y: x.union(y), [cl.cells for cl in clusters])
        return Cluster(all_cells)


    def _backwardTracking(self, new_track_ids):
        """
        Backward tracking: 
        @param new_track_ids list of track Ids to which 
                             _forwardTracking associated
                             the clusters
        @return list of new track Ids that will need to be fused
        """
        # go through each new track and see if the track should 
        # be merged with another track. Two tracks are tagged for a fuse if 
        # cluster at time t - dt is inside the group of clusters at time t

        # create big clusters for each of the track_ids that are present at
        # the previous time step
        old_big_clusters = {}
        for track_id in range(self.getNumberOfTracks()):
        	cluster_ids = self.cluster_connect[track_id].get(self.t_index - 1, None)
        	if cluster_ids:
        		old_big_clusters[track_id] = self.getBigClusterAt(track_id, self.t_index - 1)

        # find the tracks to fuse
        new_track_ids_to_fuse = []

        # iterate over each the tracks the new clusters belong to
        for new_track_id in new_track_ids:

            # big cluster in new_track_id at the present time
            big_cluster = self.getBigClusterAt(new_track_id, self.t_index)

            track_ids_to_fuse = set()
            for track_id in old_big_clusters:

                if track_id == new_track_id:
                    # skip self
                    continue

                # get the big cluster in track_id at the previous time
                old_big_cluster = old_big_clusters[track_id]
#                if old_big_clusters[track_id].isClusterInsideOf(big_cluster, frac=0.8):
                if old_big_clusters[track_id].isCentreInsideOfExt(big_cluster):
                	# tag this cluster for later fuse with new_track_id
                    track_ids_to_fuse.add(track_id)

            if track_ids_to_fuse:
            	# add new_track_id to the set
                track_ids_to_fuse.add(new_track_id)
                # store
                new_track_ids_to_fuse.append(track_ids_to_fuse)

        return new_track_ids_to_fuse


    def _fuseAll(self, new_track_ids_to_fuse):
        """
        Apply fuse method to a list of track branches to merge
        @param new_track_ids_to_fuse list of sets
        """

        # create union of tracks to merge and tage the merged tracks
        # removal
        delete_elem = set()
        n = len(new_track_ids_to_fuse)
        for i in range(n):
            li = new_track_ids_to_fuse[i]
            for j in range(i + 1, n):
                lj = new_track_ids_to_fuse[j]
                if li.intersection(lj):
                    # there is an intersection, join the two
                    new_track_ids_to_fuse[i] = li.union(lj)
                    # and tag lj for removal
                    delete_elem.add(j)

        # remove the tagged elements, working our way backwards
        delete_elem = list(delete_elem)
        delete_elem.sort(reverse=True)
        for i in delete_elem:
            del new_track_ids_to_fuse[i]

        # now fuse
        delete_track_ids = []
        for ids_to_fuse in new_track_ids_to_fuse:
            delete_track_ids += self.fuse(ids_to_fuse)

        # now delete the track Ids that were fused
        delete_track_ids.sort(reverse=True)
        for i in delete_track_ids:
            del self.cluster_connect[i]



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


    def getCells(self, track_id, t_index):
        """
        Get the cluster cellss
        @param track_id track Id
        @param t_index time index
        @return the i and j cells as separate arrays
        """
        ijs = []
        for cl in self.getClusters(track_id, t_index):
            ijs += list(cl.cells)
        iis = np.array([ij[0] for ij in ijs])
        jjs = np.array([ij[1] for ij in ijs])
        return iis, jjs


    def removeTracksByValidMask(self, valid_mask, frac):
        """
        Remove the tracks that never overlap with the valid mask
        @param valid_mask is 1 over the regions where we keep the tracks
        @param frac threshold for the min fraction of overlap (0 <= frac <= 1)
        """

        # store the track ids that need to be removed
        remove_track_ids = []

        # iterate over each track id
        for track_id in range(self.getNumberOfTracks()):

            found_overlap = False

            # iterate over each time step
            for t_index in self.cluster_connect[track_id]:

                iis, jjs = self.getCells(track_id, t_index)

                numOverlap = valid_mask[iis, jjs].sum()
                if numOverlap >= frac*len(iis):
                    found_overlap = True
                    break
                

            if not found_overlap:
                remove_track_ids.append(track_id)
        
        # remove the tracks that were tagged for removal
        # walking our way back from the end. Should we remove 
        # the clusters first?
        remove_track_ids.sort(reverse=True)
        for track_id in remove_track_ids:
            del self.cluster_connect[track_id]


    def writeFile(self, filename, unit, lat, lon, i_minmax=[], j_minmax=[]):
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
        num_i = iMax - iMin 
        f.createDimension('lat', size=num_i)

        if j_minmax:
            jMin = min(jMin, j_minmax[0])
            jMax = max(jMax, j_minmax[1])
        num_j = jMax - jMin 
        f.createDimension('lon', size=num_j)

        # infinite dimension
        time = f.createDimension('time', size=None)

        # create variables

        i_index = f.createVariable('lat', 'f', ('lat',) ,zlib=True,complevel=9,\
                     least_significant_digit=4)
        f.variables['lat'].units='degrees_north'
        f.variables['lat'].standard_name='latitude'
        f.variables['lat'].axis='Y'
        f.variables['lat'].long_name='latitude'

        j_index = f.createVariable('lon', 'f', ('lon',), zlib=True,complevel=9,\
                     least_significant_digit=4)
        f.variables['lon'].units='degrees_east'
        f.variables['lon'].standard_name='longitude'
        f.variables['lon'].axis='X'
        f.variables['lon'].long_name='longitude'

        t_index = f.createVariable('time', 'f', ('time',), zlib=True,complevel=9,\
                     least_significant_digit=4)
        f.variables['time'].axis='T'
        f.variables['time'].calendar='standard'
        f.variables['time'].units=unit

        # check ordering!!
        precip = f.createVariable('cprec','f',('time','lat','lon'),zlib=True,complevel=9,\
                    least_significant_digit=4)
        f.variables['cprec'].girdtype='lonlat'
        f.variables['cprec'].code=999
        f.variables['cprec'].long_name='detected coastal precipitation'
        f.variables['cprec'].standard_name= 'detected_precipitation'
        f.variables['cprec'].short_name='cprec'
        f.variables['cprec'].units='mm/h'

        nb_var = f.createVariable('nb', 'i4', ('time', 'lat', 'lon'))
        f.variables['nb'].girdtype='lonlat'
        f.variables['nb'].code=0
        f.variables['nb'].long_name='identification number of the clusters'
        f.variables['nb'].standard_name= 'number of the clusters'
        f.variables['nb'].short_name='nb'
        f.variables['nb'].units='no unit'

        # write the data
        i_index[:] = lat[:] #np.arange(iMin, iMax + 1)  #lat[:]
        j_index[:] = lon[:] #np.arange(jMin, jMax + 1)  #lon[:]
        t_index[:] = self.time[:] #np.arange(0, self.t_index)

        # data buffer, check ordering!!
        data = np.zeros((self.t_index, num_i, num_j), np.int32)

        for track_id in range(self.getNumberOfTracks()):
            for t_index in range(self.t_index):
                clusters = self.getClusters(track_id, t_index)
                for cl in clusters:
                    n_cells = cl.getNumberOfCells()
                    tis = [t_index] * n_cells
                    jis = [c[1] - jMin for c in cl.cells]
                    iis = [c[0] - iMin for c in cl.cells]
                    # check ordering!!
                    data[tis, iis, jis] = track_id + 1
        # now write all the data in one go
        mask=np.zeros((self.t_index, num_i, num_j))
        mask[np.where(data != 0)] = 1
        precip[:] = self.precip*mask 
        nb_var[:] = data

        f.close()

    def findCluster(self, cluster_id):
        """
        Find the track number and time index of a cluster
        @param cluster_id cluster Id
        @return track_id, time_index
        @note returns -1, -1 if the cluster id was not found
        """
        found_track_id = -1
        found_time_index = -1
        for track_id in range(self.getNumberOfTracks()):
            for time_index, cluster_ids in self.cluster_connect[track_id].items():
                if cluster_id in cluster_ids:
                    found_track_id = track_id
                    found_time_index = time_index
        return found_track_id, found_time_index


    def showEllipses(self, track_id, time_inds=[]):
        """
        Show all the cluster ellipses of track track_id at time time_index
        @param track_id track ID
        @param time_inds time index list, use None to specify all
        """
        from matplotlib import pylab

        # colour map
        def rgb(x):
            r = max(0., math.sin(math.pi*(x - 0.5)))
            g = max(0., math.cos(math.pi*x))
            b = math.sin(math.pi*x)
            return r, g, b

        track = self.cluster_connect[track_id]
        t_inds = track.keys()
        if time_inds:
            t_inds = time_inds

        for i in range(len(t_inds)):
            ti = t_inds[i]
            x = float(i)/float(len(t_inds))
            color = rgb(x)

            # get the ellipses
            ellipses = [self.clusters[j].ellipse \
                            for j in track.get(ti, [])]

            if len(ellipses) == 0:
                print 'WARNING: no cluster at time index '.format(ti)

            for j in range(len(ellipses)):

                el = ellipses[j]
                iPts, jPts = el.getPolyline()
                iPtsExt, jPtsExt = el.getPolylineExt()
                pylab.plot(jPts, iPts, '-', c=color)
                pylab.plot(jPtsExt, iPtsExt, '.', c=color)
                xc, yc = el.getCentre()
                pylab.plot(yc, xc, '+', c=color)

        pylab.title('Track {} (netcdf {}) time frames {} -> {}'.format(track_id, 
                                                            track_id + 1,
                                                            t_inds[0]+1, 
                                                            t_inds[-1]+1))
        pylab.show()


    def save(self, filename="timeConnectedClusters.pckl"):
        """
        Save object to file 
        @param filename file name
        """
        f = open(filename, 'w')
        cPickle.dump(self, f)


    def getNumberOfTracks(self):
        """
        Get the number of tracks
        @return number
        """
        return len(self.cluster_connect)


    def getNumberOfTimeSteps(self):
        """
        Get the number of time steps
        @return number
        """
        return self.t_index



    def getClusters(self, track_id, time_index):
        """
        Get the clusters of given track Id at time index
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
    tcc.addTime([c0],0.8)
    print 'One cluster'
    print tcc

def testReduceNonOverlapping():
    c0 = Cluster({(1, 1), (2, 1), (2, 2)})
    c1 = Cluster({(1, 3), (1, 4), (2, 4)})
    cs = [c0, c1]
    print 'input list non-overlapping', cs
    reduce(cs)
    print 'output list non-overlapping', cs

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
    tcc.addTime([c0, c1], 0.8)
    print tcc
    c3 = Cluster({(4,3), (7,3)})
    tcc.addTime([c3], 0.8)
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
