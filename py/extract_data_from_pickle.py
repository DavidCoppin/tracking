#from time_connected_clusters import TimeConnectedClusters
import cPickle
import matplotlib.pyplot as mpl
import argparse

parser = argparse.ArgumentParser(description='extract indices for clusters')
parser.add_argument('-f', dest='filename', default='', help='Pickle files')
args = parser.parse_args()

for nb in len(filename):
    num = [int(s) for s in args.filename.split('_') if s.isdigit()]
    if num[0]<48:
        print 'num', num, filename

f1 = open(args.filename)
tracks = cPickle.load(f1)

#print 'len(tracks)', len(tracks)
#print 'tracks[0]', tracks[0]
#for cl in tracks[0][49]: print cl.toArray()


