'''
@description: Method leading the post-processing of pickles into global output
'''

import numpy as np
import argparse
from output_from_pickle import OutputFromPickle, readTxt, createTxt
import sys

def writeOutputPP(inputdir, outputdir, list_prefix):
    """
    @param inputdir: directory where pickles are stored
    @param outputdir: directory where outputs will be stored
    """
    id = 0
    track_id ={}
    filenames = readTxt(str(inputdir)+'filenames.txt')
    lat_tot = readTxt(str(inputdir)+'lat_tot.txt')
    lon_tot = readTxt(str(inputdir)+'lon_tot.txt')
    lat = list(map(lambda x: float(x.replace(",", "")), lat_tot))
    lon = list(map(lambda x: float(x.replace(",", "")), lon_tot))
    if len(list_prefix)==1:
        str_lat_lon = readTxt(str(inputdir)+'lat-lon_'+str(list_prefix[0])+'.txt')
        lat_lon = [int(str_lat_lon[0]), int(str_lat_lon[1]), int(str_lat_lon[2]), \
                    int(str_lat_lon[3])]
    else :
        lat_lon = [0, -1, 0, -1]
    for nb_day in xrange(len(filenames)):
        print 'write_output', filenames[nb_day]
        ofp = OutputFromPickle(nb_day, lat, lon, inputdir, outputdir, list_prefix, track_id, id)
        files = ofp.selectPickles()
        if len(files)==0:
            print 'no files in writeOutputPP'
            sys.exit()
        files2 = files.sort()
        print 'files', files, len(files)
        sys.exit()
        ofp.extractTracks(files)
        ofp.writeFile('always_unzip', filenames[nb_day], lat_lon)

        # Delete pickle that will not be used anymore
#        ofp.deletePickles()
        id = ofp.id
        track_id = ofp.track_id



#################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Post-processing')
    parser.add_argument('-i', dest='inputdir', default='pickles/', help='directory where pickles \
                         and txt files are stored')
    parser.add_argument('-o', dest='outputdir', default='output/', help='directory where netcdf \
                             are stored')
#    parser.add_argument('-p','--prefix', nargs='+', help='<Required> Set flag', required=True)
    parser.add_argument('-p', dest='prefix', nargs='+', default=['io-cm', 'pacific', 'america', \
                           'africa'])
    args = parser.parse_args()
    writeOutputPP(args.inputdir, args.outputdir, args.prefix)
