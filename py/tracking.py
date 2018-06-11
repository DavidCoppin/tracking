from netCDF4 import Dataset as nc
import argparse
import numpy as np
import matplotlib.pyplot as mpl
import time
import glob
from time_connected_clusters import TimeConnectedClusters
from feature_extractor import FeatureExtractor
from cluster import Cluster
from coastal_mapping import CoastalMapping
from output_from_pickle import OutputFromPickle
from write_output_pp import createTxt, readTxt
import configparser
import sys,os,string
import gzip
import bz2
import cPickle
from datetime import datetime,timedelta as td


def tracking(fyear, lyear, minmax_lons, minmax_lats, suffix, restart_dir,
             restart_interval, harvestPeriod):
    # are we restarting
    restart = False
    restart_file = None
    pickle_index = 0
    if restart_dir is not None:
        restart_file = os.path.join(restart_dir, "clusters_restart_%s.pkl" % suffix)
        if os.path.exists(restart_file):
            # restart the simulation where it left off
            restart = True
            print "Loading restart file: %s" % restart_file

            # load restart file
#            with open(restart_file) as fh:
            with gzip.GzipFile(restart_file) as fh:
                restart_data = cPickle.load(fh)

            # unpack
            tcc = restart_data['tcc']
            minmax_lats = restart_data['minmax_lats']
            minmax_lons = restart_data['minmax_lons']
            fyear = restart_data['fyear']
            list_filename = restart_data['list_filename']
            pickle_index = restart_data['pickle_index']

            # check finish date is after restart date
            if fyear > lyear:
                print "Error: finish date (%s) is earlier than restart date (%s)" % (lyear, fyear)
                sys.exit(1)

            print "Restarting from: %s" % fyear

            # we need to delete any pickle files that were written after this restart file was written,
            # otherwise they will be duplicated as we compute them again
            config_full = configparser.ConfigParser()
            config_full.read('config.cfg')
            C = config_full['clusters']
            targetdir = os.path.expandvars(C.get('targetdir'))
            pickles = glob.glob(os.path.join(targetdir, "%s_*" % suffix))
            for pickle in pickles:
                pickle_index_file = int(pickle.split("_")[-1])
                if pickle_index_file > pickle_index:
                    print "Deleting pickle file that was created after the restart file: %s" % pickle
                    os.remove(pickle)

            # increment the pickle index
            pickle_index += 1

    # if not restarting, create empty TimeConnectedClusters instance
    if not restart:
        tcc = TimeConnectedClusters()
        list_filename = []

    # prepare for writing restart files
    if restart_file is None or restart_interval is None:
        # don't write restart files
        restart_interval = None
        restart_file = None

    else:
        # make the restart directory, if required
        if not os.path.exists(restart_dir):
            os.makedirs(restart_dir)
        print "Saving restart files to: %s" % restart_file
        print "Saving restart files every %d days" % restart_interval

    # run tracking
    _tracking_main(tcc, list_filename, fyear, lyear, minmax_lons, minmax_lats,
                   suffix, harvestPeriod, restart_file, restart_interval,
                   pickle_index)


def _tracking_main(tcc, list_filename, fyear, lyear, minmax_lons, minmax_lats,
                   suffix, harvestPeriod, restart_file, restart_interval,
                   pickle_index):

    ##########################################################################
    # Import arguments from config.cfg or fix default
    ##########################################################################
    config_full = configparser.ConfigParser()
    config_full.read('config.cfg')
    print config_full.sections()
    C = config_full['clusters']

    # Read values from config.cfg
    data_path = os.path.expandvars(C.get('data_path'))
    lsm = os.path.expandvars(C.get('lsm_path'))
    print 'lsm', lsm
    targetdir = os.path.expandvars(C.get('targetdir'))
    reso = C.getint('reso')
    print 'reso', reso
    min_prec = C.getfloat('min_prec', 0)
    max_prec = C.getfloat('max_prec', 3.0)
    print 'min_prec, max_prec', min_prec, max_prec
    szone = C.getint('szone', 8)
    lzone = C.getint('lzone', 50)
    print 'szone, lzone', szone, lzone
    frac_mask = C.getfloat('frac_mask', 1.0)
    frac_ellipse = C.getfloat('frac_ellipse', 1.0)
    print 'frac_mask, frac_ellipse', frac_mask, frac_ellipse
    min_axis = C.getint('min_axis', 6)
    print 'min_axis', min_axis
    min_size = C.getint('min_size', 0)
    max_size = C.getint('max_size', 800000)
    print 'min_size, max_size', min_size, max_size
    max_cells = C.getint('max_cells', 4500)
    t_life = C.getint('t_life', 5)
    print 'max_cells, t_life', max_cells, t_life
    save = C.getboolean('save')
    #########################################################################

    # Get lat-lon limits of region and save it post_processing
    lon_slice = slice(minmax_lons[0], minmax_lons[1])
    lat_slice = slice(minmax_lats[0], minmax_lats[1])
    if os.path.isfile(str(targetdir)+'lat-lon_'+str(suffix)+'.txt'):
        pass
    else:
        createTxt(str(targetdir)+'lat-lon_'+str(suffix)+'.txt', [minmax_lats[0], \
                   minmax_lats[1], minmax_lons[0], minmax_lons[1]])

    # Create the two coastal masks
    cm = CoastalMapping(lsm, np.int(reso), lat_slice, lon_slice, np.int(szone), \
                         np.int(lzone), np.int(min_size), np.int(max_size))

    # difference between start and end dates
    delta = lyear - fyear

    #########################################################################
    # Loop over days
    #########################################################################
    for nb_day in xrange(delta.days + 1):
        date = fyear + td(days=nb_day)
        filename=os.path.join(str(data_path)+'Cmorph-' \
               + str(date.year) + '_' + str(date.month).zfill(2) + '_'\
               + str(date.day).zfill(2) + '.nc.bz2')

        # Open the file
        filename = filename.replace('--','-').replace('__','_')
        print 'filename', filename
        list_filename = np.append(list_filename, filename)
        zipfile = bz2.BZ2File(filename)
        data_unzip = zipfile.read()
        newfilename = filename[:-4]
        open(newfilename, 'wb').write(data_unzip)
        try:
            f = nc(newfilename)
        except RuntimeError:
            f = nc(newfilename.replace('-','_'))
        if lon_slice.start < lon_slice.stop:
            all_data = f.variables["CMORPH"][:, lat_slice, lon_slice]
        else:
            all_data1 = f.variables["CMORPH"][:, lat_slice,lon_slice.start:]
            all_data2 = f.variables["CMORPH"][:, lat_slice,:lon_slice.stop]
            all_data = np.concatenate((all_data1, all_data2), axis=2)

        # Store once and for all info: lat, lon for post-processing
        if nb_day == 0:
            lat = f.variables['lat'][minmax_lats[0]:minmax_lats[1]]
            if os.path.isfile(str(targetdir)+'lon_tot.txt'):
                pass
            else:
                # Store for zone
                lat_tot_zone = f.variables['lat'][minmax_lats[0]:minmax_lats[1]]
                lon_tot_zone = f.variables['lon'][minmax_lons[0]:minmax_lons[1]]
                createTxt(str(targetdir)+'lat_tot_'+str(suffix)+'.txt', lat_tot_zone)
                createTxt(str(targetdir)+'lon_tot_'+str(suffix)+'.txt', lon_tot_zone)
                # Store for whole globe
                lat_tot = f.variables['lat'][:]
                lon_tot = f.variables['lon'][:]
                createTxt(str(targetdir)+'lat_tot.txt', lat_tot)
                createTxt(str(targetdir)+'lon_tot.txt', lon_tot)

            # Special case of zones on both sides of 0 degree of longitude
            if lon_slice.start < lon_slice.stop:
                lon = f.variables['lon'][minmax_lons[0]:minmax_lons[1]]
            else:
                lon1 = f.variables['lon'][minmax_lons[0]:]
                lon2 = f.variables['lon'][:minmax_lons[1]]
                lon = np.concatenate((lon1, lon2))
        f.close()

        # Begin tracking
        i_minmax = (0, len(lat))
        j_minmax = (0, len(lon))
        timesteps = np.shape(all_data)[0]
        for t in xrange(timesteps):
            print 'nb_day, t', nb_day, t
            data = all_data[t]

            # Extract clusters with watershed and remove large-scale clusters
            clusters = FeatureExtractor(data, thresh_low=min_prec, thresh_high=max_prec, \
                           mask=np.flipud(cm.lArea), frac=frac_mask).getClusters(min_axis)

            # Check time connectivity between clusters
            tcc.addTime(clusters,frac_ellipse)

            # Harvest the dead tracks and write to file
            if harvestPeriod and (t + 1) % harvestPeriod == 0:
                tcc.harvestTracks(targetdir+suffix, i_minmax, j_minmax, np.flipud(cm.sArea), frac_mask,
                                  max_cells, t_life*timesteps, pickle_index, dead_only=True)
#            if t==2:
#                sys.exit()
        os.remove(newfilename)
        del all_data, data, clusters, data_unzip

        # store restart file?
        if restart_interval is not None and (nb_day + 1) % restart_interval == 0:
            # create object for dumping to file
            restart_data = {
                'tcc': tcc,
                'minmax_lons': minmax_lons,
                'minmax_lats': minmax_lats,
                'fyear': date + td(days=1),
                'list_filename': list_filename,
                'pickle_index': pickle_index,
            }

            # if the restart file already exists we temporarily rename it in case writing the new one fails
            tmpRestart = False
            if os.path.exists(restart_file):
                os.rename(restart_file, restart_file + '.tmp')
                tmpRestart = True

            # create the restart file
            print "Writing restart file:", restart_file
#            with open(restart_file, "wb") as fh:
            with gzip.GzipFile(restart_file, "wb") as fh:
                cPickle.dump(restart_data, fh)

            # remove temporary file
            if tmpRestart:
                os.remove(restart_file + '.tmp')

            # increment pickle_index, so we know which pickle files were written after the restart file
            # we will have to delete these when restarting, otherwise they will be duplicated
            pickle_index += 1

    # Final harvest (all tracks)
    print "final harvest (pickle index is %d)" % pickle_index
    tcc.harvestTracks(targetdir+suffix, i_minmax, j_minmax, np.flipud(cm.sArea), frac_mask,
                      max_cells, t_life*timesteps, pickle_index, dead_only=True)

    # Save filenames for post-processing:
    createTxt(str(targetdir)+'filenames.txt', list_filename)

    if save:
        tcc.save('cmorph.pckl_'+str(suffix))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tracking.')
    parser.add_argument('-w', action='store_true', help='Print copyright')
    parser.add_argument('-d1', dest='date1', default='2010-02-19', help='Start date YYYY-MM-DD')
    parser.add_argument('-d2', dest='date2', default='2010-02-21', help='End date YYYY-MM-DD')
    parser.add_argument('-lons', dest='lons', default='1700:2200', help='Min and max longitude \
                           indices LONMIN,LONMAX')
    parser.add_argument('-lats', dest='lats', default='200:500', help='Min and max latitude \
                           indices LATMIN,LATMAX')
    parser.add_argument('-suffix', dest='suffix', default='', help='Suffix for output')
    parser.add_argument('-harvest', dest='harvestPeriod', type=int, default=10,
                         help='Number of time steps before dead tracks area saved to disk (0 for no harvest)')
    parser.add_argument('-restart-dir', default=None, help="Directory for storing and loading restart files")
    parser.add_argument('-restart-interval', type=int, default=None, help="If set then write restart files at this interval (days)")
    args = parser.parse_args()

    # get the lat-lon box
    print np.array(args.lons.split(':'))
    try:
        minmax_lons = np.array(args.lons.split(':')).astype(np.int)
    except:
        raise RuntimeError, 'Wrong specification of longitude bound indices, use -lons LONMIN:LONMAX'
    try:
        minmax_lats = np.array(args.lats.split(':')).astype(np.int)
    except:
        raise RuntimeError, 'Wrong specification of latitude bound indices, use -lats LATMIN:LATMAX'
    try:
        fyear = datetime.strptime(args.date1,'%Y-%m-%d')
        lyear = datetime.strptime(args.date2,'%Y-%m-%d')
    except IndexError,ValueError:
        sys.stdout.write(helpstring+'\n')
        sys.exit()
    tracking(fyear, lyear, minmax_lons, minmax_lats, args.suffix, args.restart_dir,
             args.restart_interval, args.harvestPeriod)
