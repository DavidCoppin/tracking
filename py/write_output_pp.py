'''
@description: Method leading the post-processing of pickles into global output
'''

import numpy as np
import argparse
from output_from_pickle import OutputFromPickle, readTxt, createTxt
from datetime import datetime
import sys, os
import cPickle

def writeOutputPP(date, inputdir, outputdir, list_prefix, suffix, restart_dir):
    """
    @param inputdir: directory where pickles are stored
    @param outputdir: directory where outputs will be stored
    """
    restart = False
    restart_file = None
    id = 0
    track_id ={}
    dict_pickles = {}
    filenames = readTxt(str(inputdir)+'filenames.txt')
    print 'filenames', filenames
    if len(list_prefix)==1:
        lat_tot = readTxt(str(inputdir)+'lat_tot_'+str(list_prefix[0])+'.txt')
        lon_tot = readTxt(str(inputdir)+'lon_tot_'+str(list_prefix[0])+'.txt')
        str_lat_lon = readTxt(str(inputdir)+'lat-lon_'+str(list_prefix[0])+'.txt')
        lat_lon = [int(str_lat_lon[0]), int(str_lat_lon[1]), int(str_lat_lon[2]), \
                              int(str_lat_lon[3])]
    else:
        lat_tot = readTxt(str(inputdir)+'lat_tot.txt')
        lon_tot = readTxt(str(inputdir)+'lon_tot.txt')
        lat_lon = np.zeros(4)
    lat = list(map(lambda x: float(x.replace(",", "")), lat_tot))
    lon = list(map(lambda x: float(x.replace(",", "")), lon_tot))
    if restart_dir is not None:
        restart_file = os.path.join(restart_dir, "info_pp.pkl")
        if os.path.exists(restart_file):
            # restart the post_processing where it left off
            f = open(restart_file)
            restart_data = cPickle.load(f)

            # unpack
            lastname = restart_data['filename']
            offset = restart_data['nb_day']
            track_id = restart_data['track_id']
            id = restart_data['id']
            print 'restart lastname, offset, track_id, id', lastname, offset, track_id, id
            print 'date', date
            for i in range(len(filenames)):
                if str(date) in str(filenames[i]):
                    print 'found in str(filenames[i])', str(filenames[i]), i
            sys.exit()

    for nb_day in xrange(len(filenames)):
        print 'write_output', filenames[nb_day]
        ofp = OutputFromPickle(nb_day, lat, lon, inputdir, outputdir, list_prefix, dict_pickles,\
                                track_id, id)
        files = ofp.selectPickles()
        if len(files)==0:
            print 'no files in writeOutputPP'
        files2 = files.sort()
        print 'files', files, len(files)
        ofp.extractTracks(files)
        ofp.writeFile(str(suffix), filenames[nb_day], lat_lon)

        # Delete pickle that will not be used anymore
#        ofp.deletePickles()
        dict_pickles = ofp.dict_pickles
        track_id = ofp.track_id
        id = ofp.id

    # save informations about post-processing
    ofp.saveInfoPP(filenames[nb_day], nb_day)


#################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Post-processing')
    parser.add_argument('-d', dest='date', default='2010_02_21', help='End date YYYY_MM_DD')
    parser.add_argument('-i', dest='inputdir', default='pickles/', help='directory where pickles \
                         and txt files are stored')
    parser.add_argument('-o', dest='outputdir', default='output/', help='directory where netcdf \
                             are stored')
#    parser.add_argument('-p','--prefix', nargs='+', help='<Required> Set flag', required=True)
    parser.add_argument('-p', dest='prefix', nargs='+', default=['io-cm', 'pacific', 'america', \
                           'africa'])
    parser.add_argument('-s', dest='suffix', default='final', help='suffix for netcdf file')
    parser.add_argument('-restart_dir', default=None, help="Directory for storing and loading \
                               restart files")
    args = parser.parse_args()
#    try:
#        lyear = datetime.strptime(args.date,'%Y-%m-%d')
#    except IndexError,ValueError:
#        sys.stdout.write(helpstring+'\n')
#        sys.exit()
    writeOutputPP(args.date, args.inputdir, args.outputdir, args.prefix, args.suffix, \
                   args.restart_dir)
