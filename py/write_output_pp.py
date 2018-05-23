import numpy as np
import argparse
from output_from_pickle import OutputFromPickle, readTxt, createTxt
import sys

def writeOutputPP(inputdir, outputdir):
    id = 0
    track_id ={}
    filenames = readTxt(str(inputdir)+'filenames.txt')
    lat_tot = readTxt(str(inputdir)+'lat_tot.txt')
    lon_tot = readTxt(str(inputdir)+'lon_tot.txt')
    lat = list(map(lambda x: float(x.replace(",", "")), lat_tot))
    lon = list(map(lambda x: float(x.replace(",", "")), lon_tot))
#   list_prefix = ['io-cm', 'pacific', 'america', 'africa']
#    list_prefix = ['io-cm', 'pacific', 'america']
#    list_prefix = ['africa']
    list_prefix = ['png']
    for nb_day in xrange(len(filenames)):
        print 'write_output', filenames[nb_day]
        ofp = OutputFromPickle(nb_day, lat, lon, inputdir, outputdir, list_prefix, track_id, id)
        files = ofp.selectPickles()
        if len(files)==0:
            print 'no files in writeOutputPP'
            sys.exit()
        files2 = files.sort()
        print 'files', files, len(files)
        ofp.extractTracks(files)
        ofp.writeFile('final', filenames[nb_day])

        # Delete pickle that will not be used anymore
        ofp.deletePickles()
        id = ofp.id
        track_id = ofp.track_id



#################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Post-processing')
    parser.add_argument('-i', dest='inputdir', default='pickles/', help='directory where pickles \
                         and txt files are stored')
    parser.add_argument('-o', dest='outputdir', default='output/', help='directory where netcdf \
                             are stored')
    args = parser.parse_args()
    writeOutputPP(args.inputdir, args.outputdir)
