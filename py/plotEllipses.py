import argparse
import pickle

parser = argparse.ArgumentParser(description='plot ellipses of time connected clusters object')
parser.add_argument('-n', dest='track_id', type=int, default=0, help='Track ID')
parser.add_argument('-t', dest='time_index', type=int, default=-1, help='Time index or -1 for all times')
parser.add_argument('-f', dest='filename', default='timeConnectedClusters.pckl', help='Pickle file name')
args = parser.parse_args()

f = open(args.filename)
tcc = pickle.load(f)
tindx = args.time_index
if tindx < 0:
	tindx = None
tcc.showEllipses(args.track_id, tindx)



