import numpy
import matplotlib
import time_connected_clusters
from cluster import Cluster
import matplotlib.pyplot as plt
import random
import copy


"""
Test reduce 
"""

def plot(clusters_before, clusters_after, title):
    for cl in clusters_before:
        iPts, jPts= cl.ellipse.getPolyline()
        x, y = cl.ellipse.getCentre()
        plt.plot(iPts, jPts, 'g-', [x], [y], 'g+')
    for cl in clusters_after:
        iPts, jPts= cl.ellipse.getPolyline()
        x, y = cl.ellipse.getCentre()
        plt.plot(iPts, jPts, 'r-', [x], [y], 'r+')
    plt.title(title)
    plt.show()


def testReduce1():
    """
    All the clusters merge into a single cluster
    """
    rect1 = {(2+i, 3) for i in range(10)}.union({(2+i, 4) for i in range(10)})
    # same rectangle but slightly shifted to the right
    rect2 = {(3+i, 3) for i in range(10)}.union({(3+i, 4) for i in range(10)})
    # now shift up
    rect3 = {(3+i, 4) for i in range(10)}.union({(3+i, 5) for i in range(10)})
    clusters = [Cluster(rect1), Cluster(rect2), Cluster(rect3)]

    clusters_before = copy.deepcopy(clusters)
    time_connected_clusters.reduce(clusters)
    plot(clusters_before, clusters, 'testReduce1')

    assert(len(clusters) == 1)

def testReduce2():
    """
    Three clusters merge into a 2 clusters
    """
    rect1 = {(2+i, 3) for i in range(10)}.union({(2+i, 4) for i in range(10)})
    # same rectangle but slightly shifted to the right
    rect2 = {(3+i, 3) for i in range(10)}.union({(3+i, 4) for i in range(10)})
    # now shift up
    rect3 = {(3+i, 5) for i in range(10)}.union({(3+i, 5) for i in range(10)})
    clusters = [Cluster(rect1), Cluster(rect2), Cluster(rect3)]

    clusters_before = copy.deepcopy(clusters)
    time_connected_clusters.reduce(clusters)
    plot(clusters_before, clusters, 'testReduce2')

    assert(len(clusters) == 2)

def testReduce3():
    """
    One cluster is just a tiny bit outside
    """
    rect1 = {(2+i, 3) for i in range(10)}.union({(2+i, 4) for i in range(10)})
    # same rectangle but slightly shifted to the right
    rect2 = {(3+i, 3) for i in range(10)}.union({(3+i, 4) for i in range(10)})
    # now shift up
    rect3 = {(3+i, 4) for i in range(10)}.union({(3+i, 5.1) for i in range(10)})
    clusters = [Cluster(rect1), Cluster(rect2), Cluster(rect3)]

    clusters_before = copy.deepcopy(clusters)
    time_connected_clusters.reduce(clusters)
    plot(clusters_before, clusters, 'testReduce3')

    assert(len(clusters) == 2)

def testReduce4():
    """
    One cluster is slightly inside
    """
    rect1 = {(2+i, 3) for i in range(10)}.union({(2+i, 4) for i in range(10)})
    # same rectangle but slightly shifted to the right
    rect2 = {(3+i, 3) for i in range(10)}.union({(3+i, 4) for i in range(10)})
    # now shift up
    rect3 = {(3+i, 4) for i in range(10)}.union({(3+i, 5.05) for i in range(10)})
    clusters = [Cluster(rect1), Cluster(rect2), Cluster(rect3)]

    clusters_before = copy.deepcopy(clusters)
    time_connected_clusters.reduce(clusters)
    plot(clusters_before, clusters, 'testReduce4')

    assert(len(clusters) == 1)


def testLarge(numClusters):
    """
    Many clusters at random locations
    """
    random.seed(1234)
    template = {(2+i, 3) for i in range(10)}.union({(2+i, 4) for i in range(10)})
    clusters = []
    for i in range(numClusters):
        offset = 10*random.random(), 10*random.random()
        cells = {(offset[0] + c[0], offset[1] + c[1]) for c in template}
        clusters.append(Cluster(cells))

    clusters_before = copy.deepcopy(clusters)
    time_connected_clusters.reduce(clusters)

    print 'testLarge: num clusters before {} after {}'.format(len(clusters_before), len(clusters))

    plot(clusters_before, clusters, 'testLarge')




if __name__ == '__main__':
    testReduce1()
    testReduce2()
    testReduce3()
    testReduce4()
    testLarge(20)

