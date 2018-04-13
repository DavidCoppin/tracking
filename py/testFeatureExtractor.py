import numpy
from feature_extractor import FeatureExtractor
from cluster import Cluster

def testABC():
    print '=' * 70
    print 'testABC'
    print '-' * 70

    a = {(i, i) for i in range(10)}
    a = a.union({(i + 10, 10 - i) for i in range(10)})
    a = a.union({(5 + i, 4) for i in range(10)})

    b = {(40, i) for i in range(10)}
    b = b.union({(40 + i, 10) for i in range(4)})
    b = b.union({(40 + i, 5) for i in range(4)})
    b = b.union({(40 + i, 0) for i in range(4)})
    b = b.union({(44, 6 + i) for i in range(4)})
    b = b.union({(44, 1 + i) for i in range(4)})

    c = {(21, 20 + i) for i in range(9)}
    c = c.union({(21, 30 + i) for i in range(9)})
    c = c.union({(20, 21 + i) for i in range(8)})

    ca = Cluster(a)
    cb = Cluster(b)
    cc = Cluster(c)

    bounds = [(0, 0), (60, 100)]

    ia, ij, aa = ca.toArray(bounds)
    ib, jb, ab = cb.toArray(bounds)
    ic, jc, ac = cc.toArray(bounds)

    # combine the letters, setting zeros and ones
    data = numpy.logical_or(aa, ab, ac)
    nz = numpy.sum(data)
    ntot = data.shape[0] * data.shape[1]
    print 'data has {} non-zeros values out of {}'.format(nz, ntot)

    fe = FeatureExtractor(data, thresh_low=0, thresh_high=4)
    cs = fe.getClusters()
    print 'fe', fe
    print 'number of features: {}'.format(len(cs))
    for c in cs:
        print c
    assert(len(cs) == 3)



if __name__ == '__main__':
    testABC()
