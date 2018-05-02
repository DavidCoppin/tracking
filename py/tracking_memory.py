from netCDF4 import Dataset as nc
import argparse
import numpy as np
import matplotlib.pyplot as mpl
import time
from time_connected_clusters import TimeConnectedClusters
from feature_extractor import FeatureExtractor
from cluster import Cluster
from coastal_mapping import CoastalMapping
from output_file import OutputFile
import configparser
import sys,os,string
import bz2
from datetime import datetime,timedelta as td

@profile
def tracking(fyear, lyear, minmax_lons, minmax_lats, suffix):

    ### Need to choose if we import parameters as arguments from testCmorph or if we use
    ### configuration.py to do that for us
    config_full = configparser.ConfigParser()
    config_full.read('config.cfg')
    print config_full.sections()
    C = config_full['clusters']  # config parser expects sections
    # Read values from config.cfg
    lsm = os.path.expandvars(C.get('lsm_path'))
    print 'lsm', lsm
    reso = C.getint('reso')
    print 'reso', reso
    min_prec = C.getfloat('min_prec', 0)
    max_prec = C.getfloat('max_prec', 3)
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
    #####################################################

    lon_slice = slice(minmax_lons[0], minmax_lons[1])
    lat_slice = slice(minmax_lats[0], minmax_lats[1])
    # Get the two coastal masks
    cm = CoastalMapping(lsm, np.int(reso), lat_slice, lon_slice, np.int(szone), \
                         np.int(lzone), np.int(min_size), np.int(max_size))
#    mpl.contourf(np.flipud(cm.sArea))
#    mpl.savefig('mask_'+str(suffix)+'.png')
#    mpl.show()
    llat = minmax_lats[1] - minmax_lats[0]
    llon = minmax_lons[1] - minmax_lons[0]
    tcc = TimeConnectedClusters()
    of = OutputFile(tcc)
    delta = lyear - fyear
    dates = [fyear + td(days=i) for i in xrange(delta.days + 1)]
    precip = np.zeros((llat,llon))
    print 'dates', dates
    list_filename=[]
    for nb_day in xrange(len(dates)):
        date=dates[nb_day]
        filename=os.path.join('Data/CMORPH/Cmorph-' \
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
        all_data = f.variables["CMORPH"][:, lat_slice, lon_slice]
        all_time = f.variables["time"][:]
        # Store once and for all lat, lon, unit
        if nb_day == 0:
            lat = f.variables['lat'][minmax_lats[0]:minmax_lats[1]]
            lon = f.variables['lon'][minmax_lons[0]:minmax_lons[1]]
            unit = f.variables["time"].units
        f.close()
        for t in xrange(len(all_time)) :
            print 'nb_day, t', nb_day, t
            data = all_data[t]
            # Extract clusters with watershed and remove large-scale clusters
            clusters = FeatureExtractor(data, thresh_low=min_prec, thresh_high=max_prec, \
                           mask=np.flipud(cm.lArea), frac=frac_mask).getClusters(min_axis)
            tcc.addTime(clusters,frac_ellipse)
        of.getTime(all_time)
        os.remove(newfilename)
        del all_data, data, clusters, data_unzip
    # Remove tracks in large mask but never in small
    tcc.removeTracksByValidMask(valid_mask=np.flipud(cm.sArea), frac=frac_mask)
    # get 3D array of clusters from TimeConnectedClusters
    tracks = tcc.toArray(of.time, i_minmax=(0, len(lat)), j_minmax=(0, len(lon)))
    for nb_day in xrange(len(dates)):
        print 'write_output, nb_day', nb_day
        name = list_filename[nb_day]
        ini = nb_day*48
        end = (nb_day+1)*48
        of.writeFile(str(suffix), list_filename[nb_day], tracks[ini:end], lat, lon, ini, end, \
                      unit, lat_slice, lon_slice)
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
    parser.add_argument('-suffix', dest='suffix', default='', help='suffix for output')
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
    tracking(fyear,lyear,minmax_lons, minmax_lats, args.suffix)
