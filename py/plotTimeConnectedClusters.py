import argparse
import pickle

parser = argparse.ArgumentParser(description='plot content of time connected clusters object')
parser.add_argument('-i', dest='track_id', type=int, default=0, help='Track ID')
parser.add_argument('-t', dest='time_index', type=int, default=0, help='Time index or None for all times')
parser.add_argument('-p', dest='filename', default='timeConnectedClusters.pckl', help='Pickle file name)
args = parser.parse_args()

f = open(args.filename)
tcc = pickle.load(f)
tcc.showEllipses(args.track_id, args.time_index)



