import numpy
import matplotlib
import time
from time_connected_clusters import TimeConnectedClusters
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
    tcc.writeFile('rectangle.nc', i_minmax=(0, 10), j_minmax=(0, 8))


def testIndRectangles():
    rect1 = {(2, 3), (3, 3), (4, 3), (2, 4), (3, 4), (4, 4)}
    rect2 = {(6, 3), (7, 3), (8, 3), (6, 4), (7, 4), (8, 4)}
    rect3 = {(3, 4), (4, 4), (5, 4), (3, 5), (4, 5), (5, 5)}

    tcc = TimeConnectedClusters()
    tcc.addTime([Cluster(rect1), Cluster(rect2)])
    print tcc
    tcc.addTime([Cluster(rect3)])
    print tcc
    tcc.writeFile('independant_rectangles.nc', i_minmax=(0, 10), j_minmax=(0, 8))


def testTwoMergingRectangles():
    print '='*70
    print 'testTwoMergingRectangles'
    print '-'*70
    rect0 = {(2, 3), (3, 3), (4, 3), (2, 4), (3, 4), (4, 4),(2, 5), (3, 5), (4, 5)}
    rect1 = {(6, 3), (7, 3), (8, 3), (6, 4), (7, 4), (8, 4),(6, 5), (6, 5), (6, 5)}
    rect2 = {(3, 3), (4, 3), (5, 3), (6, 3), (7, 3), (3, 4), (4, 4), (5, 4), (6, 4), 
             (7, 4), (3, 5), (4, 5), (5, 5), (6, 5), (7, 5)}
    c0, c1, c2 = Cluster(rect0), Cluster(rect1), Cluster(rect2)
    if c0.isCentreInsideOf(c2):
        print 'c0 is inside c2'
    if c1.isCentreInsideOf(c2):
        print 'c1 is inside c2'
    if c1.isCentreInsideOf(c0):
        print 'c1 is inside c0'
    if c2.isCentreInsideOf(c0):
        print 'c2 is inside c0'
    # Should have only one track because cluster 0 and 1 at t=0 merge at t=1. Should 
    # pass into fuse    # Prob: cluster 0 at=0 does not become cluster 1
    tcc = TimeConnectedClusters()
    tcc.addTime([c0, c1])
    print tcc
    tcc.addTime([c2])
    print tcc
    tcc.writeFile('two_merging_rectangles.nc', i_minmax=(0, 10), j_minmax=(0, 8))
    assert(tcc.getNumberOfTracks() == 1)


def testOnlyFuse():
    print '='*70
    print 'testOnlyFuse'
    print '-'*70
    rect0 = {(2, 3), (3, 3), (2, 4), (3, 4), (2, 5), (3, 5)}
    rect1 = {(7, 3), (8, 3), (7, 4), (8, 4), (7, 5), (8, 5)}
    rect2 = {(2, 3), (3, 3), (4, 3), (5, 3), (6, 3), (7, 3), (8, 3), 
             (2, 4), (3, 4), (4, 4), (5, 4), (6, 4), (7, 4), (8, 4), 
             (2, 5), (3, 5), (4, 5), (5, 5), (6, 5), (7, 5), (8, 5)}
    c0 = Cluster(rect0)
    c1 = Cluster(rect1)
    c2 = Cluster(rect2)
    # Simplest fuse test: no forward tracking possible (center of rect3 can't 
    # be inside ellipse of cluster 1 or 2, but centers of clusters 1 and 2 should be 
    # inside ellipse of cluster 3
    # Result expected: all clusters should get same id: 1 ==> connectivity [{0: [0, 1], 1: [2]}]
    # Problem : only c2 appears in netcdf file at t=0 while it should be c0 and c1. 
    # Probably overwritten (lat,lon) of c0 and c1 by (lat,lon) of c2
    tcc = TimeConnectedClusters()
    print 'time step {}: adding cluster with centre {} and {}'.format(tcc.getNumberOfTimeSteps(),
                                                                       c0.getCentre(), 
                                                                       c1.getCentre())
    tcc.addTime([c1, c0])
    print tcc

    if c2.isCentreInsideOf(c0):
        print 'c2 is inside c0'
    if c2.isCentreInsideOf(c1):
        print 'c2 is inside c1'
    if c0.isCentreInsideOf(c2):
        print 'c0 is inside c2'
    if c1.isCentreInsideOf(c2):
        print 'c1 is inside c2'
    print 'time step {}: adding cluster with centre {}'.format(tcc.getNumberOfTimeSteps(),
                                                                       c2.getCentre())
    tcc.addTime([c2])
    print tcc
    tcc.writeFile('only_fuse.nc', i_minmax=(0, 10), j_minmax=(0, 8))
    assert(tcc.getNumberOfTracks() == 1)
    assert(tcc.getNumberOfTimeSteps() == 2)


def testOnlySplit():
    rect1 = {(2, 3), (3, 3), (2, 4), (3, 4), (2, 5), (3, 5)}
    rect2 = {(7, 3), (8, 3), (7, 4), (8, 4), (7, 5), (8, 5)}
    rect3 = {(2, 3), (3, 3), (4, 3), (5, 3), (6, 3), (7, 3), (8, 3), (2, 4), (3, 4), (4, 4), (5, 4), (6, 4), (7, 4), (8, 4), (2, 5), (3, 5), (4, 5), (5, 5), (6, 5), (7, 5), (8, 5)}
    # Simplest split test : no backward tracking possible (center of rect3 can't be inside ellipse of cluster 1 or 2, but centers of clusters 1 and 2 should be inside ellipse of cluster 3
    # Result expected: all clusters should get same id: 1 ==> connectivity [{0: [0], 1: [1, 2]}]
    # Problem : cluster 2 keeps its id
    tcc = TimeConnectedClusters()
    tcc.addTime([Cluster(rect3)])
    print tcc
    tcc.addTime([Cluster(rect1), Cluster(rect2)])
    print tcc
    tcc.writeFile('only_split.nc', i_minmax=(0, 10), j_minmax=(0, 8))


def testSplittingInTwo():
    rect1 = {(3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9)}
    circ2 = {(3, 8), (3, 9), (4, 7), (4, 8), (4, 10), (5, 7), (5, 8), (5, 9), (5, 10), (6, 8), (6, 9)}
    circ3 = {(3, 3), (4, 2), (4, 3), (4, 4), (5, 3)}
    tcc = TimeConnectedClusters()
    tcc.addTime([Cluster(rect1)])
    print tcc
    tcc.addTime([Cluster(circ2), Cluster(circ3)])
    print tcc
    tcc.writeFile('splitting_in_two.nc', i_minmax=(0, 10), j_minmax=(0, 8))
    # expected result: 
    # Because centre of circ2 and circ3 are inside ellipse of rect1, the result
    # should be [{0: [rect1], 1: [circ2, circ3]} 
    # Currently fusing, test if rect1's centre is inside ellipse of circ2 or circ3 (should not be). 
    # But even if it is, fuse should be active because it will try to fuse with itself


def extractClusters(data, thresh_min, thresh_max):
    """ 
    Extract clusters from an image data
    @param data
    @param thresh_min
    @param thresh_max
    @return list of clusters
    """
    # remove data below minimum threshold
    ma_data = numpy.ma.masked_where(data <= thresh_min, data)
    # building black and white image with lower threshold to create borders for watershed
    tmp_data = ma_data.filled(fill_value=0)
    tmp_data[numpy.where(tmp_data !=0)] = 255 
    bw_data = tmp_data.astype(numpy.uint8)
    border = cv2.dilate(bw_data, None, iterations=5)
    border -= cv2.erode(border, None)

    # remove data below minimum threshold
    ma_conv = numpy.ma.masked_where(data <= thresh_max, data)
    # building black and white image with high threshold to serve as markers for watershed..
    tmp_conv = ma_conv.filled(fill_value=0)
    tmp_conv[numpy.where(tmp_conv !=0)] = 255 
    bw_conv = tmp_conv.astype(numpy.uint8)
    markers = ndimage.label(bw_conv, structure=numpy.ones((3, 3)))[0]
    # add border on image with high threshold to tell the watershed where it should fill in
    markers[border == 255] = 255 
    # labels each feature
    labels = watershed(-data, markers, mask=bw_data)

    # load into clusters
    res = []
    for idVal in range(1, labels.max()+1):
        iVals, jVals = numpy.where(labels == idVal)
        numVals = len(iVals)
        if numVals > 0:
            cells = [(iVals[i], jVals[i]) for i in range(len(iVals))]
            # store this cluster as a list with one element (so far). Each.
            # element will have its own ID
            res.append(Cluster(cells))
    return res


def testOverlap():
    print '='*70
    print 'testOverlap'
    print '-'*70
    # Expected result: since c0 is inside c1, would expect/like that it is directly included
    # in track with c1
    rect0 = {(4, 3), (5, 3), (6, 3)}
    rect1 = {(3, 3), (7, 3), (3, 4), (4, 4), (5, 4), (6, 4),
             (7, 4), (3, 5), (4, 5), (5, 5), (6, 5), (7, 5)}
    c0, c1 = Cluster(rect0), Cluster(rect1)
    if c0.isCentreInsideOf(c1):
        print 'c0 is inside c1'
    if c1.isCentreInsideOf(c0):
        print 'c1 is inside c0'
    tcc = TimeConnectedClusters()
    tcc.addTime([c0, c1])
    print tcc
    tcc.addTime([c0])
    print tcc
    tcc.writeFile('overlap.nc', i_minmax=(0, 10), j_minmax=(0, 8))
    assert(tcc.getNumberOfTracks() == 1)


def testDigits():
    """
    Checking that we can create a time connected cluster from image
    """
    lats = numpy.arange(0, 20)
    lons = numpy.arange(0, min_area)
    data = numpy.zeros((len(lats), len(lons)), numpy.float64)
    # create a cluster
    data[10:10+5, 10] = 5
    data[2:4, 5:5+4] = 3
    data[4:6, 5:5+4] = 1
    data[6:8, 5:5+4] = 3
    data[3, 12] = 5
    data[5, 15] = 2
#    print data
    # Version with Class
#    clusters = extractClusters.extractPoints(numpy.flipud(data), thresh_min=0., thresh_max=2.5)
    # Version inside testTimeConnecteClusters
    clusters = extractClusters(numpy.flipud(data), thresh_min=0., thresh_max=2.5)
    print clusters

    tcc = TimeConnectedClusters()

    # add a couple of times
    tcc.addTime(clusters)
    #tcc.addTime(clusters)

    print tcc

    # write to file
    tcc.writeFile('digits.nc', i_minmax=(0, len(lats)), j_minmax=(0, len(lons)))


def testSeveralFuse():
    print '='*70
    print 'testSeveralFuse'
    print '-'*70
    rect0 = {(2, 3), (3, 3), (2, 4), (3, 4), (2, 5), (3, 5)}
    rect1 = {(7, 3), (8, 3), (7, 4)}
    rect2 = {(8, 4), (7, 5), (8, 5)}
    rect3 = {(7, 3), (8, 3), (7, 4), (8, 4), (7, 5), (8, 5)}
    rect4 = {(2, 3), (3, 3), (4, 3), (5, 3), (6, 3), (7, 3), (8, 3),
             (2, 4), (3, 4), (4, 4), (5, 4), (6, 4), (7, 4), (8, 4),
             (2, 5), (3, 5), (4, 5), (5, 5), (6, 5), (7, 5), (8, 5)}
    rect5 = {(10, 5), (11, 5), (12, 5), (10, 6), (11, 6), (12, 6)}
    rect6 = {(10, 7), (11, 7), (12, 7), (10, 6), (11, 6), (12, 6)}
    c0 = Cluster(rect0)
    c1 = Cluster(rect1)
    c2 = Cluster(rect2)
    c3 = Cluster(rect3)
    c4 = Cluster(rect4)
    c5 = Cluster(rect5)
    c6 = Cluster(rect6)

    # Description: test different fuses at different time steps
    # Expected result: 4 tracks after t0, 3 tracks after t1, 2 after t2 with
    #       connectivity [{0: [0, 2], 1: [5]}, {0: [1], 1: [4]}, {0: [3], 1: [6]}] at t1
    #       connectivity [{0: [0, 2, 1], 1: [5, 4], 2: [7]}, {0: [3], 1: [6], 2: [8]}] at t2

    #
    # time_index = 0
    tcc = TimeConnectedClusters()
    print 'time step {}: adding clusters with centres {} {} {} {}'.format(tcc.getNumberOfTimeSteps(),
                                                                         c1.getCentre(),
                                                                         c0.getCentre(),
                                                                         c2.getCentre(),
                                                                         c5.getCentre())
    tcc.addTime([c1, c0, c2, c5])
    assert(tcc.getNumberOfTracks() == 4)
    print tcc

    #
    # time_index = 1
    print 'time step {}: adding clusters with centres {} {} {}'.format(tcc.getNumberOfTimeSteps(),
                                                                         c0.getCentre(),
                                                                         c3.getCentre(),
                                                                         c6.getCentre())
    tcc.addTime([c0, c3, c6])
    assert(tcc.getNumberOfTracks() == 3)
    print tcc
    tr_id5, t_indx5 = tcc.findCluster(5)
    tr_id4, t_indx4 = tcc.findCluster(4)
    tr_id6, t_indx6 = tcc.findCluster(6)
    print 'cluster 5 is in track {} at time index {}'.format(tr_id5, t_indx5)
    print 'cluster 4 is in track {} at time index {}'.format(tr_id4, t_indx4)
    print 'cluster 6 is in track {} at time index {}'.format(tr_id6, t_indx6)
    assert(t_indx5 == t_indx4 == t_indx6 == 1)
    assert(tr_id5 != tr_id4)
    assert(tr_id4 != tr_id6)
    assert(tr_id5 != tr_id6)

    #
    # time_index = 2
    print 'time step {}: adding cluster with centres {} {}'.format(tcc.getNumberOfTimeSteps(),
                                                                     c4.getCentre(),
                                                                     c5.getCentre())
    c7 = c4
    c8 = c5
    tcc.addTime([c7, c8])
    #       connectivity [{0: [0, 2, 1], 1: [5, 4], 2: [7]}, {0: [3], 1: [6], 2: [8]}] at t2
    assert(tcc.getNumberOfTracks() == 2)
    tr_id0, t_indx0 = tcc.findCluster(0)
    tr_id1, t_indx1 = tcc.findCluster(1)
    tr_id2, t_indx2 = tcc.findCluster(2)
    tr_id3, t_indx3 = tcc.findCluster(3)
    tr_id4, t_indx4 = tcc.findCluster(4)
    tr_id5, t_indx5 = tcc.findCluster(5)
    tr_id6, t_indx6 = tcc.findCluster(6)
    tr_id7, t_indx7 = tcc.findCluster(7)
    tr_id8, t_indx8 = tcc.findCluster(8)
    assert(t_indx8 == 2 and t_indx7 == 2)
    assert(tr_id7 != tr_id8)
    assert(tr_id5 == tr_id4 and t_indx5 == t_indx4 == 1)
    assert(tr_id6 != tr_id5 and t_indx6 == 1)
    # this fails
    assert((tr_id0 == tr_id1 == tr_id2) and (t_indx0 == t_indx1 == t_indx2))

    print tcc
    tcc.writeFile('several_fuse.nc', i_minmax=(0, 15), j_minmax=(0, 10))


def testProbMinFuse():
    print '='*70
    print 'testProbMinFuse'
    print '-'*70
    # Description: test fusing bug between colliding clusters
    # works fine with rect0 = {(7, 3), (8, 3)} and rect1 = {(7, 4), (8, 4), (7, 5), (8, 5)}
    # but does not work with 
    # rect0 = {(7, 4), (8, 4), (7, 5), (8, 5)} and 
    #rect1 = {(7, 3), (8, 3)},
    # which is just replacing rect0 by rect1 and vice versa
    # --> case ok: c2 inside c1 but not inside c0
    #     bug: c2 inside c0 but not inside c1
    # ==> problem probably when appending track in addTime: after tcc.addTime([c2]), track and time
    # for c1 are switched to -1  (cf print cluster 1 is in track -1 at time index -1)
    # Expected result: 2 tracks after t0, 1 tracks after t1 with
    #       connectivity [{0: [0, 1], 1: [2]}] at t1
    rect0 = {(7, 4), (8, 4), (7, 5), (8, 5)}
    rect1 = {(7, 3), (8, 3)}
    rect2 = {(7, 3), (8, 3), (7, 4), (8, 4), (7, 5), (8, 5)}
    c0 = Cluster(rect0)
    c1 = Cluster(rect1)
    c2 = Cluster(rect2)
    if c0.isCentreInsideOf(c2):
        print 'c0 is inside c2'
    if c1.isCentreInsideOf(c2):
        print 'c1 is inside c2'
    if c2.isCentreInsideOf(c0):
        print 'c2 is inside c0'
    if c2.isCentreInsideOf(c1):
        print 'c2 is inside c1'

    tcc1 = TimeConnectedClusters()
    tcc2 = TimeConnectedClusters()
    print 'time step {}: adding clusters with centres {} {}'.format(tcc1.getNumberOfTimeSteps(),
                                                                         c0.getCentre(),
                                                                         c1.getCentre())
    tcc1.addTime([c0, c1])
    tcc2.addTime([c1, c0]) # change the order
    assert(tcc1.getNumberOfTracks() == 2)
    assert(tcc2.getNumberOfTracks() == 2)
    tr_id0, t_indx0 = tcc1.findCluster(0)
    tr_id1, t_indx1 = tcc1.findCluster(1)
    print 'cluster 0 is in track {} at time index {}'.format(tr_id0, t_indx0)
    print 'cluster 1 is in track {} at time index {}'.format(tr_id1, t_indx1)
    print '>>> tcc1'
    print tcc1
    print '<<< tcc2'
    print tcc2

    #
    # time_index = 1
    print 'time step {}: adding clusters with centres {}'.format(tcc1.getNumberOfTimeSteps,
                                                                         c2.getCentre())
    tcc1.addTime([c2])
    tcc2.addTime([c2])
    tr_id0, t_indx0 = tcc1.findCluster(0)
    tr_id1, t_indx1 = tcc1.findCluster(1)
    tr_id2, t_indx2 = tcc1.findCluster(2)
    print 'cluster 0 is in track {} at time index {}'.format(tr_id0, t_indx0)
    print 'cluster 1 is in track {} at time index {}'.format(tr_id1, t_indx1)
    print 'cluster 2 is in track {} at time index {}'.format(tr_id2, t_indx2)
    print '>>> tcc1'
    print tcc1
    print '<<< tcc2'
    print tcc2
    # run some checks
    assert(tcc1.getNumberOfTracks() == tcc2.getNumberOfTracks() == 1)
    assert(t_indx0 == 0)
    assert(t_indx1 == 0)
    assert(t_indx2 == 1)
    assert(tr_id0 == 0 and tr_id1 == 0 and tr_id2 == 0)
    tcc1.writeFile('prob_min_fuse1.nc', i_minmax=(0, 15), j_minmax=(0, 10))
    tcc2.writeFile('prob_min_fuse2.nc', i_minmax=(0, 15), j_minmax=(0, 10))


def testSplitMulti():
    print '='*70
    print 'testSplitMulti'
    print '-'*70
    rect0 = {(2, 3), (3, 3), (4, 3), (5, 3), (6, 3), (7, 3), (8, 3),
             (2, 4), (3, 4), (4, 4), (5, 4), (6, 4), (7, 4), (8, 4),
             (2, 5), (3, 5), (4, 5), (5, 5), (6, 5), (7, 5), (8, 5),
             (2, 6), (3, 6), (4, 6), (5, 6), (6, 6), (7, 6), (8, 6)}
    rect1 = {(2, 3), (3, 3), (2, 4), (3, 4), (2, 5), (3, 5), (2, 6), (3, 6)}
    rect2 = {(7, 3), (8, 3), (7, 4)}
    rect3 = {(8, 4), (7, 5), (8, 5), (7, 6), (8, 6)}
    rect4 = {(2, 3), (3, 3), (2, 4), (2, 2), (3, 2)}
    rect5 = {(3, 5), (2, 6), (3, 6), (2, 7), (3, 7)}
    rect6 = {(7, 3), (8, 3), (7, 4), (7, 2), (8, 2), (7, 1)}
    rect7 = {(8, 5), (7, 6), (8, 6), (7, 7), (8, 7)}
    c0 = Cluster(rect0)
    c1 = Cluster(rect1)
    c2 = Cluster(rect2)
    c3 = Cluster(rect3)
    c4 = Cluster(rect4)
    c5 = Cluster(rect5)
    c6 = Cluster(rect6)
    c7 = Cluster(rect7)
    # Description: test consecutive splits at different time steps
    # Expected result: 1 tracks all along
    tcc = TimeConnectedClusters()
#    print 'time step {}: adding clusters with centres {}'.format(tcc.getNumberOfTimeSt>
#                                                                         c0.getCentre())
    tcc.addTime([c0])
    print tcc
#    print 'time step {}: adding clusters with centres {} {} {}'.format(tcc.getNumberOfTimeSteps>
#                                                                         c1.getCentre(),
#                                                                         c2.getCentre(),
#                                                                         c3.getCentre())
    tcc.addTime([c1, c2, c3])
    print tcc
#    print 'time step {}: adding cluster with centres {} {} {} {}'.format(tcc.getNumberOfTimeSteps(),
#                                                                     c4.getCentre(),
#                                                                     c5.getCentre(),
#                                                                     c6.getCentre(),
#                                                                     c7.getCentre())
    tcc.addTime([c4, c5, c6, c7])
    print tcc
    tcc.writeFile('split_multi.nc', i_minmax=(0, 15), j_minmax=(0, 10))
    assert(tcc.getNumberOfTracks() == 1)


def testMovingClusters():
    print '='*70
    print 'testMovingClusters'
    print '-'*70
    # Expected result: after new ellipse fitted to min area, c1, c2 should be inside c0
    # and c4, c5, c6 inside c3 ==> 2 tracks expected
    rect0 = {(3, 2), (4, 2), (3, 3), (4, 3)}
    rect1 = {(6, 2), (7, 2), (6, 3), (7, 3)}
    rect2 = {(5, 4), (6, 4), (5, 5), (6, 5)}
    rect3 = {(10, 10), (11, 10)}
    rect4 = {(10, 13), (11, 13)}
    rect5 = {(12, 12), (13, 12)}
    rect6 = {(13, 10), (14, 10)}
    min_axis = 6 # 10000 # 50
    c0, c1, c2 = Cluster(rect0, min_ellipse_axis=min_axis), Cluster(rect1, min_ellipse_axis=min_axis), Cluster(rect2, min_ellipse_axis=min_axis)
    c3, c4, c5, c6 = Cluster(rect3, min_ellipse_axis=min_axis), Cluster(rect4, min_ellipse_axis=min_axis), Cluster(rect5, min_ellipse_axis=min_axis), Cluster(rect6, min_ellipse_axis=min_axis)
    if c1.isCentreInsideOf(c0):
        print 'c1 is inside c0'
    if c2.isCentreInsideOf(c0):
        print 'c2 is inside c0'
    if c4.isCentreInsideOf(c3):
        print 'c4 is inside c3'
    if c5.isCentreInsideOf(c3):
        print 'c5 is inside c3'
    if c6.isCentreInsideOf(c3):
        print 'c6 is inside c3'
    # Should have only one track because cluster 0 and 1 at t=0 merge at t=1. Should.
    # pass into fuse    # Prob: cluster 0 at=0 does not become cluster 1
    tcc = TimeConnectedClusters()
    tcc.addTime([c0, c3])
    print tcc
    tcc.addTime([c1, c2, c4, c5, c6])
    print tcc
    tcc.writeFile('moving_clusters.nc', i_minmax=(0, 16), j_minmax=(0, 16))
    assert(tcc.getNumberOfTracks() == 2)


if __name__ == '__main__':
    #testRectangle()
    #testIndRectangles()
    #testTwoMergingRectangles()
    #testOnlyFuse()
    #testOnlySplit()
    #testSplittingInTwo()
    #testSplitMulti()
    #testSeveralFuse()
    #testProbMinFuse()
    #testOverlap()
    testMovingClusters()

