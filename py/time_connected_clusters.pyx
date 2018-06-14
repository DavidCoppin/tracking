#cython: profile=False

from cluster import Cluster
from ellipse import Ellipse
import numpy as np
import netCDF4
import copy
import math
import cPickle
import functools
import tempfile
import os
import gzip
import sys


## it is faster to pass in the values instead of the arrays
## the rotation matrix has the form [[tr00, tr01], [-tr01, tr00]]
## it's important that this function be inlined with the calling function
## for performance
## only need the top row of the rotation matrix if we uyse symetry
cdef bint _isPointInsideEllipse(double a, double b,
	                       double tr00, double tr01,
                           double centreX, double centreY,
                           double pointX, double pointY):
    """
    Check if a point is inside an ellipse
    @param a x radius of ellipse in rotated coordinates
    @param b y radius of ellipse in rotated coordinates
    @param tr00 element of rotation matrix
    @param tr01 element of rotation matrix
    @param centreX x coordinate of centre
    @param centreY y coordinate of centre
    @param pointX x cpoordinat of point
    @param pointY y coordinate of point
    @return True if point is inside, False otherwise
    """
    pointX -= centreX
    pointY -= centreY
    # rotate the coordinates to align them to the principal axes
    cdef double ptXPrime =(+tr00 * pointX + tr01 * pointY) / a
    cdef double ptYPrime =(-tr01 * pointX + tr00 * pointY) / b

    return (ptXPrime*ptXPrime + ptYPrime*ptYPrime < 1.0)


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
        eli = cli.ellipse
        eli_transf = eli.ij2AxesTransf
        eli_centre = eli.centre
        eli_a = eli.a
        eli_b = eli.b
        for j in range(i + 1, n):
            clj = cluster_list[j]
            elj = clj.ellipse
            elj_transf = elj.ij2AxesTransf
            elj_centre = elj.centre

            isCliInsideClj = _isPointInsideEllipse(elj.a, elj.b, elj_transf[0,0], elj_transf[0,1],
                                                   elj_centre[0], elj_centre[1], 
                                                   eli_centre[0], eli_centre[1])
            isCljInsideCli = _isPointInsideEllipse(eli_a, eli_b, eli_transf[0,0], eli_transf[0,1],
                                                   eli_centre[0], eli_centre[1], 
                                                   elj_centre[0], elj_centre[1])


#            if cli.isCentreInsideOf(clj) and clj.isCentreInsideOf(cli):
            if isCliInsideClj and isCljInsideCli:
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
        # number of clusters
        self.num_clusters = 0

        # list of dictionaries. Each dictionary represents clusters across time which share the same
        # ID. Each element is a track {t_index0: {'area':area0, 'clusters':[cluster0, cluster1,..]},
        # t_index1: { ...}} where area and clusters are keys in a dictionary for each t_index.
        # t_index# is the time index and cluster# is a Cluster instance, and area# is the total area
        # covered by a track at time #
        self.cluster_connect = []

        # current time index
        self.t_index = 0


    def addTime(self, new_clusters, frac, frac_decrease):
        """
        Add time entry
        @param new_clusters list of new clusters
        @param frac TO DESCRIBE
        @param frac_decrease: decrease in size above which a track is cut (in percent) 
        """
        # merge overlapping clusters, this will reduce the number of
        # clusters to track but should have no influence on the
        # final result
        reduce(new_clusters, frac)

        # special case if first time step
        if self.t_index == 0:
            # there is no past
            # just create a new track for each cluster
            for new_cl in new_clusters:
                area = new_cl.getNumberOfCells()
                self.cluster_connect.append({self.t_index: {'area': area, 'clusters': [new_cl]}})
                self.num_clusters += 1

            # done
            self.t_index += 1
            return


        new_track_ids = self._forwardTracking(new_clusters)

        new_track_ids_to_fuse = self._backwardTracking(new_track_ids)

        # reduce list of tracks to fuse by looking at common track Ids between
        # the elements of new_track_ids_to_fuse. These elements will be tagged for
        # removal in delete_elem
        self._fuseAll(new_track_ids_to_fuse)

        # cut tracks when large decrease in area
        self.cutTracks(frac_decrease)

        # done with assigning, update the time index
        self.t_index += 1


    def _forwardTracking(self, new_clusters):
        """
        Forward tracking: assign new clusters to existing tracks
        @param new_clusters new clusters
        @return set of track Ids to which the new clusters belong to
        """

        # set of track Ids to which the new clusters will be assigned to
        new_track_ids = set() # new cluster index to track id

        for new_cl_index in range(len(new_clusters)):
            new_cl = new_clusters[new_cl_index]
            new_el = new_cl.ellipse
            new_transf = new_el.ij2AxesTransf[0, :]
            new_centre = new_el.centre
            new_aExt = new_el.aExt
            new_bExt = new_el.bExt

            # the track Id that we need to assign this cluster to
            new_track_id = -1

            # find out if this cluster belongs to an existing track. This is equivalent to
            # forward tracking
            num_tracks = self.getNumberOfTracks()
            connected_clusters = []
            connected_track_ids = []
            for track_id in range(num_tracks):
                if self.t_index-1 in self.cluster_connect[track_id]:
                    # track exists at t_index - 1
                    old_clusters = self.cluster_connect[track_id][self.t_index - 1].get('clusters')
                else :
                    old_clusters = []
                for old_cl in old_clusters:
                    old_el = old_cl.ellipse
                    old_transf = old_el.ij2AxesTransf[0, :]
                    old_centre = old_el.centre
                    # is the centre of new_cl inside the ellipse of old_cl?
                    isNewClInsideOldCl = _isPointInsideEllipse(old_el.aExt, old_el.bExt, 
                    	                                       old_transf[0], old_transf[1],
                                                               old_centre[0], old_centre[1],
                                                               new_centre[0], new_centre[1])
                    # is the centre of old_cl inside the ellipse of new_cl?
                    isOldClInsideNewCl = _isPointInsideEllipse(new_aExt, new_bExt, 
                    	                                       new_transf[0], new_transf[1],
                                                               new_centre[0], new_centre[1],
                                                               old_centre[0], old_centre[1])


#                    if new_cl.isCentreInsideOfExt(old_cl) or old_cl.isCentreInsideOfExt(new_cl):
                    if isNewClInsideOldCl or isOldClInsideNewCl:
                        connected_clusters.append(old_cl)
                        connected_track_ids.append(track_id)

            area = new_cl.getNumberOfCells()
            if len(connected_track_ids) == 0:
                # this cluster is on its own
                # create a new entry
                new_track_id = self.getNumberOfTracks()
                self.cluster_connect.append({self.t_index: {'area': area, 'clusters': [new_cl]}})
            else:
                # choose the track for which the distance between new_cl any of the clusters
                # of that track at t - dt is smallest
                dists = np.array([new_cl.getDistance(cl) for cl in connected_clusters])
                i = np.argmin(dists)
                new_track_id = connected_track_ids[i]
                if self.t_index in self.cluster_connect[new_track_id]:
                    # time index is already in the track, append the new cluster and update area
                    self.cluster_connect[new_track_id][self.t_index]['clusters'].append(new_cl)
#                    self.cluster_connect[new_track_id][self.t_index]['clusters'] = \
#                            self.cluster_connect[new_track_id][self.t_index]['clusters'] + [new_cl]
                    self.cluster_connect[new_track_id][self.t_index]['area'] = \
                            self.cluster_connect[new_track_id][self.t_index].get('area') + area
                else:
                    # create a new entry
                    self.cluster_connect[new_track_id][self.t_index] = {'area': area, \
                            'clusters': [new_cl]}
#                print 'self.cluster_connect[new_track_id]', self.cluster_connect[new_track_id]
            new_track_ids.add(new_track_id)

            # update number of clusters
            self.num_clusters += 1

        return new_track_ids


    def getBigClusterAt(self, track_id, t_index):
        """
        Construct a big cluster from all the track_id clusters at time t_index
        @param track_id track Id
        @param t_index time index
        @return one big cluster representing the merge of all smaller clusters
        """
        clusters = self.cluster_connect[track_id][t_index].get('clusters', [])
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
            if self.t_index-1 in self.cluster_connect[track_id]:
                  old_big_clusters[track_id] = self.getBigClusterAt(track_id, self.t_index - 1)

        # find the tracks to fuse
        new_track_ids_to_fuse = []

        # iterate over each the tracks the new clusters belong to
        for new_track_id in new_track_ids:

            # big cluster in new_track_id at the present time
            big_cluster = self.getBigClusterAt(new_track_id, self.t_index)
            big_ellipse = big_cluster.ellipse
            big_transf = big_ellipse.ij2AxesTransf
            big_centre = big_ellipse.centre
            big_aExt = big_ellipse.aExt
            big_bExt = big_ellipse.bExt

            track_ids_to_fuse = set()
            for track_id in old_big_clusters:

                if track_id == new_track_id:
                    # skip self
                    continue

                # get the big cluster in track_id at the previous time
                old_big_cluster = old_big_clusters[track_id]
                obc_centre = old_big_cluster.ellipse.centre

                # is the old big cluster inside the big clusters ellipse?
                isOBCInsideBC = _isPointInsideEllipse(big_aExt, big_bExt, big_transf[0,0], big_transf[0,1],
                                                           big_centre[0], big_centre[1],
                                                           obc_centre[0], obc_centre[1])

#                if old_big_clusters[track_id].isClusterInsideOf(big_cluster, frac=0.8):
#                if old_big_clusters[track_id].isCentreInsideOfExt(big_cluster):
                if isOBCInsideBC:
                    # tag this cluster for later fuse with new_track_id
                    track_ids_to_fuse.add(track_id)

            if track_ids_to_fuse:
            	# add new_track_id to the set
                track_ids_to_fuse.add(new_track_id)
                # store
                new_track_ids_to_fuse.append(track_ids_to_fuse)

        return new_track_ids_to_fuse


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
                if t_index in self.cluster_connect[track_id] and t_index not in cl_conn0 :
                    cl_conn0[t_index] = {'area': cl_list['area'], 'clusters': cl_list['clusters']}
                elif t_index in self.cluster_connect[track_id] and t_index in cl_conn0 :
                    cl_conn0[t_index]['area'] = cl_conn0[t_index]['area'] + cl_list['area']
                    cl_conn0[t_index]['clusters'] = cl_conn0[t_index]['clusters'] + cl_list['clusters']
                else :
                    pass

        # return the track Ids marked for deletion
        n = len(tr_ids)
        return [tr_ids[i] for i in range(n - 1, 0, -1)]


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


    def cutTracks(self, frac):
        """
        Cut tracks when area decreases by more than frac %. If cut, end previous track 
        and create a new one
        @param frac: fraction
        """
        num_tracks = self.getNumberOfTracks()
        for track_id in range(num_tracks):
            # check if track exist for last 2 time steps
            if self.t_index in self.cluster_connect[track_id] and \
                    self.t_index -1 in self.cluster_connect[track_id] :
                # compare area
                if self.cluster_connect[track_id][self.t_index-1]['area'] > 490 and \
                        (self.cluster_connect[track_id][self.t_index]['area'] < \
#               if (self.cluster_connect[track_id][self.t_index]['area'] < \
                        (1-frac) * self.cluster_connect[track_id][self.t_index - 1]['area']):
                    # create new track with track at t_index and delete it at t_index -1
                    print 'self.cluster_connect[track_id][self.t_index-1][area]', \
                            self.cluster_connect[track_id][self.t_index-1]['area']
                    print 'self.cluster_connect[track_id][self.t_index][area]', \
                            self.cluster_connect[track_id][self.t_index]['area']
                    self.cluster_connect.append({self.t_index: {'area': \
                            self.cluster_connect[track_id][self.t_index]['area'], 'clusters': \
                            self.cluster_connect[track_id][self.t_index]['clusters']}})
                    del self.cluster_connect[track_id][self.t_index]


    def getMinMaxIndices(self):
        """
        Get the low/high end box indices
        @return i_min, j_min, i_max, j_max
        """
        i_min, j_min, i_max, j_max = self.LARGE_INT, self.LARGE_INT, -self.LARGE_INT, -self.LARGE_INT
        for track in self.cluster_connect:
            for cl_list in track.values():
                i_min = min(i_min, min([cl.box[0][0] for cl in cl_list['clusters']]))
                j_min = min(j_min, min([cl.box[0][1] for cl in cl_list['clusters']]))
                i_max = max(i_max, max([cl.box[1][0] for cl in cl_list['clusters']]))
                j_max = max(j_max, max([cl.box[1][1] for cl in cl_list['clusters']]))
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


    def getTrackSize(self, track_id, time_index):
        num_cells = 0
        for cluster in self.cluster_connect[track_id][time_index].get('clusters', []):
            num_cells += cluster.getNumberOfCells()
        return num_cells


    def getMaxTrackArea(self, track_id):
        """
        Get the maximum area of a track
        @param track_id: track id
        @return the maximum area
        """
        all_track = self.cluster_connect[track_id].values()
        hi = -sys.maxint-1
        for x in (item['area'] for item in all_track):
            hi = max(x,hi)
        return hi


    def removeTrack(self, track_id):
        """
        Remove a track
        @param track_id track Id
        """
        self.cluster_connect.pop(track_id)


    def getStartEndTimes(self, track_id):
        """
        Get the start/end time indices
        @param track_id track Id
        @return t_beg, t_end
        """
        t_inds = self.cluster_connect[track_id].keys()
        t_inds.sort()
        return t_inds[0], t_inds[-1]


    def getTrackSize(self, track_id, time_index):
        num_cells = 0
        for cluster in self.cluster_connect[track_id][time_index].get('clusters', []):
            num_cells += cluster.getNumberOfCells()
        return num_cells


    def checkNoSynoptic(self, max_cells, length_time, track_id):
        """
        Get the maximum Area along a track
        @param track_id: track Id
        @return max_area
        """
        no_synoptic = True
        # Remove long AND big clusters
        length = len(self.cluster_connect[track_id])
        if length < length_time:
            pass
        else :
            max_area = self.getMaxTrackArea(track_id)
            if max_area >= max_cells :
                no_synoptic = False
                print 'length_time, max_cells', length_time, max_cells
                print 'no_synoptic, length, max_area', no_synoptic, length, max_area
        return no_synoptic


    def harvestTracks(self, prefix, i_minmax, j_minmax, mask, frac, max_cells, \
                       length_time, pickle_index, dead_only=False):
        """
        Harvest tracks and remove from list
        @param prefix to be prepended to the file name
        @param i_minmax min/max lat indices
        @param j_minmax min/max lon indices
        @param dead_only only harvest tracks that are no longer alive, otherwise
                         harvest all the tracks
        """
        t_index_min = self.LARGE_INT
        t_index_max = -self.LARGE_INT
        tracks_to_harvest = []
        good_tracks_to_harvest = []
        for track_id in range(len(self.cluster_connect)):
            t_beg, t_end = self.getStartEndTimes(track_id)
            if not dead_only or t_end < self.t_index - 1:
                t_index_min = min(t_index_min, t_beg)
                t_index_max = max(t_index_max, t_end)
                tracks_to_harvest.append(track_id)
                # keep only tracks that are above islands at some time
                # and that are not synoptic
                if dead_only or (self.checkTrackOverMask(mask, frac, track_id) \
                         and self.checkNoSynoptic(max_cells, length_time, track_id)):
                    good_tracks_to_harvest.append(track_id)

        # write the tracks to file
#        print "harvesting tracks: pickle index: %d" % pickle_index
        prfx = prefix + '_{}_{}_'.format(t_index_min, t_index_max)
        sufx = '_%d' % pickle_index
        num_times = t_index_max - t_index_min + 1
#        print 'good_tracks_to_harvest', good_tracks_to_harvest
        self.saveTracks(good_tracks_to_harvest, num_times, i_minmax, j_minmax,
                        prfx, suffix=sufx)

        # remove the harvested tracks
        tracks_to_harvest.sort(reverse=True)
#        print '... harvesting and removing tracks {}'.format(good_tracks_to_harvest)
#        print '... removing tracks {}'.format(tracks_to_harvest)
        for track_id in tracks_to_harvest:
            self.removeTrack(track_id)


    def checkTrackOverMask(self, mask, frac, track_id):
        """
        Keep the tracks the overlap with the mask
        @param valid_mask is 1 over the regions where we keep the tracks
        @param frac threshold for the min fraction of overlap (0 <= frac <= 1)
        """
        found_overlap = False
        for t_index in self.cluster_connect[track_id]:
            iis, jjs = self.getCells(track_id, t_index)
            numOverlap = mask[iis, jjs].sum()
            if numOverlap >= frac*len(iis):
                found_overlap = True
                break
        return found_overlap


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
        # walking our way back from the end.
        remove_track_ids.sort(reverse=True)
        for track_id in remove_track_ids:
            self.removeTrack(track_id)



    def toArray(self, num_times, i_minmax=[], j_minmax=[], track_id_list=None):
        """
        Convert clusters in tcc into 3D array
        @param num_times number of time indices
        @param i_minmax min/max lat indices
        @param j_minmax min/max lon indices
        @param track_id_list list of track Ids, use None to select all tracks
        """

        # get the sizes
        iMin, jMin, iMax, jMax = self.getMinMaxIndices()

        if i_minmax:
            iMin = min(iMin, i_minmax[0])
            iMax = max(iMax, i_minmax[1])
        num_i = iMax - iMin #+ 1

        if j_minmax:
            jMin = min(jMin, j_minmax[0])
            jMax = max(jMax, j_minmax[1])
        num_j = jMax - jMin #+ 1

        track_ids = [i for i in range(self.getNumberOfTracks())]
        if track_id_list:
            track_ids = track_id_list

        # data buffer, check ordering!!
        data = np.zeros((num_times, num_i, num_j), np.int32)

        for track_id in track_ids:
            for time_index in range(num_times):
                clusters = self.getClusters(track_id, time_index)
                for cl in clusters:
                    n_cells = cl.getNumberOfCells()
                    tis = [time_index] * n_cells
                    jis = [c[1] - jMin for c in cl.cells]
                    iis = [c[0] - iMin for c in cl.cells]
                    # check ordering!!
                    data[tis, iis, jis] = track_id + 1
        return data


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
                if cluster_id in cluster_ids['clusters']:
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
            ellipses = [cl.ellipse for cl in track[ti].get('clusters', [])]

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


    def saveTracks(self, track_id_list, num_times, i_minmax, j_minmax, prefix, suffix=''):
        """
        Save the tracks to file
        @param track_id_list list of track Ids
        @param num_times number of time indices
        @param i_minmax min/max lat indices
        @param j_minmax min/max lon indices
        @param prefix prefix of the file
        """
        if not track_id_list:
            # nothing to do
            return
        # toArray seems to produce that just pickling the clusters
        #data = self.toArray(num_times, i_minmax, j_minmax, track_id_list)
        data = [self.cluster_connect[tid] for tid in track_id_list]
        with tempfile.NamedTemporaryFile(prefix=prefix, dir=os.getcwd(), delete=False, suffix=suffix) as f:
#            print "saving tracks to file: %s" % f.name
            with gzip.GzipFile(fileobj=f) as gzf:
                cPickle.dump(data, gzf)


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
        return self.cluster_connect[track_id][time_index].get('clusters', [])


    def __repr__(self):
        res = """
TimeConnectedCluster: num of clusters   {}
                      num of time steps {}
                      num of tracks     {}
                      connectivity      {}
        """.format(self.num_clusters, self.t_index, \
              len(self.cluster_connect), \
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
    print 'Two clusters'
    print tcc

if __name__ == '__main__':
    testNoCluster()
    testOneCluster()
    testReduceNonOverlapping()
    testReduceOverlapping()
    testTwoClustersAtTime0()
