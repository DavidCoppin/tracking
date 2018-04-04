import argparse
import cPickle

parser = argparse.ArgumentParser(description='plot ellipses of time connected clusters object')
parser.add_argument('-n', dest='track_id', type=int, default=0, help='Track ID')
parser.add_argument('-t', dest='time_inds', type=str, default="", help='Time indices')
parser.add_argument('-f', dest='filename', default='timeConnectedClusters.pckl', help='Pickle file name')
args = parser.parse_args()

f = open(args.filename)
tcc = cPickle.load(f)
tinds = []
if args.time_inds:
	t_min, t_max = args.time_inds.split(':')
	tinds = range(int(t_min), int(t_max))
tcc.showEllipses(args.track_id, tinds)



