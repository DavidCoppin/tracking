from netCDF4 import Dataset as nc
import argparse
import numpy as np
import matplotlib.pyplot as mpl
import time
from time_connected_clusters import TimeConnectedClusters
from feature_extractor import FeatureExtractor
from cluster import Cluster
from coastal_mapping import CoastalMapping
from output_from_pickle import OutputFromPickle
from write_output_pp import createTxt, readTxt
import configparser
import sys,os,string
import bz2
from datetime import datetime,timedelta as td

def tracking(fyear, lyear, minmax_lons, minmax_lats, suffix, harvestPeriod=0):

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

    # Load class used for tracking and dates for all the files
    tcc = TimeConnectedClusters()
    delta = lyear - fyear
    dates = [fyear + td(days=i) for i in xrange(delta.days + 1)]
    list_filename=[]

    #########################################################################
    # Loop over days
    #########################################################################
    for nb_day in xrange(len(dates)):
        date=dates[nb_day]
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
        for t in xrange(np.shape(all_data)[0]):
            print 'nb_day, t', nb_day, t
            data = all_data[t]

            # Extract clusters with watershed and remove large-scale clusters
            clusters = FeatureExtractor(data, thresh_low=min_prec, thresh_high=max_prec, \
                           mask=np.flipud(cm.lArea), frac=frac_mask).getClusters(min_axis)

            # Check time connectivity between clusters
            tcc.addTime(clusters,frac_ellipse)

            # Harvest the dead tracks and write to file
            if harvestPeriod and (t + 1) % harvestPeriod == 0:
                tcc.harvestTracks(prefix=targetdir+suffix, i_minmax=i_minmax, j_minmax=j_minmax, \
                                   mask=np.flipud(cm.sArea), frac=frac_mask, dead_only=True)
        os.remove(newfilename)
        del all_data, data, clusters, data_unzip

    # Final harvest (all tracks)
    tcc.harvestTracks(prefix=targetdir+suffix,i_minmax=i_minmax, j_minmax=j_minmax, \
                       mask=np.flipud(cm.sArea), frac=frac_mask)

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
    tracking(fyear,lyear,minmax_lons, minmax_lats, args.suffix, harvestPeriod=args.harvestPeriod)
