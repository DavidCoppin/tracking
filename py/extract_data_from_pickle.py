#from time_connected_clusters import TimeConnectedClusters
import cPickle
import matplotlib.pyplot as mpl

#tcc = TimeConnectedClusters()
f1 = open("cm_t0-56_5Pguf1")
tcc1 = cPickle.load(f1)
test = tcc1.toArray(56, i_minmax=(0, 600), j_minmax=(0, 500))



