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
from configuration import Configuration
#from scipy import ndimage
#from skimage.morphology import watershed
#import cv2
import sys,os,string
import bz2
from datetime import datetime,timedelta as td

def testCmorph(lsm, fyear, lyear, minmax_lons, minmax_lats, reso, min_ellipse_axis, min_prec, \
                max_prec, frac_mask, frac_ellipse, suffix, szone, lzone, min_size, max_size, \
                save=False):

    ### Need to choose if we import parameters as arguments from testCmorph or if we use 
    ### configuration.py to do that for us
    C=Configuration('config.py')
    """
    Checking that we can create a time connected cluster from image
    """
    lon_slice = slice(minmax_lons[0], minmax_lons[1])
    lat_slice = slice(minmax_lats[0], minmax_lats[1])
    # Get the two coastal masks
    cm = CoastalMapping(lsm, np.int(reso), lat_slice, lon_slice, np.int(szone), \
                         np.int(lzone), np.int(min_size), np.int(max_size))
#    mpl.contourf(np.flipud(cm.sArea))
#    mpl.show()
    llat = minmax_lats[1] - minmax_lats[0]
    llon = minmax_lons[1] - minmax_lons[0]
    tcc = TimeConnectedClusters()
    of = OutputFile()
    delta = lyear - fyear
    dates = [fyear + td(days=i) for i in xrange(delta.days + 1)]
    precip = np.zeros((llat,llon))
    print 'dates', dates
    for nb_day in xrange(len(dates)):
        date=dates[nb_day]
        filename=os.path.join('Data/CMORPH/Cmorph-' \
               + str(date.year) + '_' + str(date.month).zfill(2) + '_'\
               + str(date.day).zfill(2) + '.nc.bz2')
        # Open the file
        filename = filename.replace('--','-').replace('__','_')
        print 'filename', filename
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
        for t in xrange(len(all_time)) :
            print 'nb_day, t', nb_day, t
            data = all_data[t]
            # Extract clusters with watershed and remove large-scale clusters
            clusters = FeatureExtractor(data, thresh_low=min_prec, thresh_high=max_prec, \
            mask=np.flipud(cm.lArea), frac=frac_mask).getClusters(min_ellipse_axis)
            tcc.addTime(clusters,frac_ellipse)
        of.getPrecip(all_data, all_time)
        os.remove(newfilename)
    # write to file
    lat = f.variables['lat'][minmax_lats[0]:minmax_lats[1]]
    lon = f.variables['lon'][minmax_lons[0]:minmax_lons[1]]
    unit = f.variables["time"].units
    f.close()
    tcc.removeTracksByValidMask(valid_mask=np.flipud(cm.sArea), frac=frac
    _mask)
    of.writeFile('cmorph.nc_'+str(suffix), unit, lat, lon, i_minmax=(0, len(lat)), \
                   j_minmax=(0, len(lon)))
    if save:
        tcc.save('cmorph.pckl_'+str(suffix))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tracking.')
    parser.add_argument('-w', action='store_true', help='Print copyright')
    parser.add_argument('-save', dest='save', action='store_true', help='Save time connected clusters \
                           object for debugging')
    parser.add_argument('-d1', dest='date1', default='2010-02-19', help='Start date YYYY-MM-DD')
    parser.add_argument('-d2', dest='date2', default='2010-02-21', help='End date YYYY-MM-DD')
    parser.add_argument('-lons', dest='lons', default='1700:2200', help='Min and max longitude \
                           indices LONMIN,LONMAX')
    parser.add_argument('-lats', dest='lats', default='200:500', help='Min and max latitude \
                           indices LATMIN,LATMAX')
    parser.add_argument('-min_axis', dest='min_axis', type=float, default=6, help='Min ellipse \
                           axis in pixels')
    parser.add_argument('-lsm', dest='lsm', default='Data/LSM/Cmorph_slm_8km.nc', help='path to \
                           land-sea data')
    parser.add_argument('-reso', dest='reso', default='8', help='resolution of the dataset in km')
    parser.add_argument('-precmin', dest='min_prec', type=float, default='0.', \
                           help='low threshold for watershed on precipitation data')
    parser.add_argument('-precmax', dest='max_prec', type=float, default='2.5', \
                           help='high threshold for watershed on precipitation data')
    parser.add_argument('-frac_mask', dest='frac_mask', type=float, default=0.8, \
                           help='threshold to keep tracks once clusters overlap with mask')
    parser.add_argument('-frac_ellipse', dest='frac_ellipse', type=float, default=0.8, \
                           help='threshold to merge overlapping ellipses')
    parser.add_argument('-suffix', dest='suffix', default='', help='suffix for output')
    parser.add_argument('-sz', dest='szone', default='6', help='small distance to coast in pixels')
    parser.add_argument('-lz', dest='lzone', default='50', help='large distance to coast in pixels')
    parser.add_argument('-smin', dest='min_size', default='300', help='area below which islands \
                           are deleted in km2')
    parser.add_argument('-smax', dest='max_size', default='800000', help='max area of filled \
                           islands in km2')
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
    testCmorph(args.lsm, fyear,lyear,minmax_lons, minmax_lats, args.reso, args.min_axis, \
                args.min_prec, args.max_prec, args.frac_mask, args.frac_ellipse, \
                args.suffix, args.szone, args.lzone, args.min_size, args.max_size, args.save)

