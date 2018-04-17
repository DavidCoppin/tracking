import numpy
import matplotlib
import time_connected_clusters
from cluster import Cluster
import matplotlib.pyplot as plt


"""
Test reduce 
"""

def plot(clusters, title):
    for cl in clusters:
        iPts, jPts= cl.ellipse.getPolyline()
        x, y = cl.ellipse.getCentre()
        plt.plot(iPts, jPts, '-', [x], [y], '+')
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

    plot(clusters, 'testReduce1: before')
    time_connected_clusters.reduce(clusters)
    plot(clusters, 'testReduce1: after')
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

    plot(clusters, 'testReduce2: before')
    time_connected_clusters.reduce(clusters)
    plot(clusters, 'testReduce2: after')
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

    plot(clusters, 'testReduce3: before')
    time_connected_clusters.reduce(clusters)
    plot(clusters, 'testReduce3: after')
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

    plot(clusters, 'testReduce4: before')
    time_connected_clusters.reduce(clusters)
    plot(clusters, 'testReduce4: after')
    assert(len(clusters) == 1)


if __name__ == '__main__':
    testReduce1()
    testReduce2()
    testReduce3()
    testReduce4()

