from netCDF4 import Dataset as nc
import argparse
import numpy
import matplotlib
import time
from time_connected_clusters import TimeConnectedClusters
from feature_extractor import FeatureExtractor
from cluster import Cluster
#from scipy import ndimage
#from skimage.morphology import watershed
#import cv2
import sys,os,string
import bz2
from datetime import datetime,timedelta as td

def testCmorph(fyear, lyear, minmax_lons, minmax_lats, min_ellipse_area):
    """
    Checking that we can create a time connected cluster from image
    """
    lon_slice = slice(minmax_lons[0], minmax_lons[1])
    lat_slice = slice(minmax_lats[0], minmax_lats[1])
    llat = minmax_lats[1] - minmax_lats[0]
    llon = minmax_lons[1] - minmax_lons[0]
    tcc = TimeConnectedClusters(min_ellipse_area=min_ellipse_area)
    delta = lyear - fyear
    dates = [fyear + td(days=i) for i in xrange(delta.days + 1)]
    print 'dates', dates
    for nb_day in xrange(len(dates)):
        date=dates[nb_day]
        filename=os.path.join('Data/CMORPH/Cmorph-' \
               + str(date.year) + '_' + str(date.month).zfill(2) + '_'\
               + str(date.day).zfill(2) + '.nc.bz2')
        #Open the file
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
        for t in xrange(48) :
            print 'nb_day, t', nb_day, t
            data = f.variables["CMORPH"][t, lat_slice, lon_slice]
            clusters = FeatureExtractor(numpy.flipud(data), thresh_min=0., thresh_max=2.5).getClusters()
#            print clusters
            tcc.addTime(clusters)

    # write to file
    lat = f.variables['lat'][minmax_lats[0]:minmax_lats[1]]
    lon = f.variables['lon'][minmax_lons[0]:minmax_lons[1]]
    f.close()
    tcc.writeFile('cmorph.nc', i_minmax=(0, len(lat)), j_minmax=(0, len(lon)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tracking.')
    parser.add_argument('-w', action='store_true', help='Print copyright')
    parser.add_argument('-d1', dest='date1', default='2010-01-12', help='Start date YYYY-MM-DD')
    parser.add_argument('-d2', dest='date2', default='2010-01-13', help='End date YYYY-MM-DD')
    parser.add_argument('-lons', dest='lons', default='1450:1700', help='Min and max longitude indices LONMIN,LONMAX')
    parser.add_argument('-lats', dest='lats', default='300:550', help='Min and max latitude indices LATMIN,LATMAX')
    parser.add_argument('-min_area', dest='min_area', default=30, help='Min ellipse area in pixels')
    args = parser.parse_args()

    # get the lat-lon box
    print numpy.array(args.lons.split(':'))
    try:
        minmax_lons = numpy.array(args.lons.split(':')).astype(numpy.int)
    except:
        raise RuntimeError, 'Wrong specification of longitude bound indices, use -lons LONMIN:LONMAX'
    try:
        minmax_lats = numpy.array(args.lats.split(':')).astype(numpy.int)
    except:
        raise RuntimeError, 'Wrong specification of latitude bound indices, use -lats LATMIN:LATMAX'
    try:
        fyear = datetime.strptime(args.date1,'%Y-%m-%d')
        lyear = datetime.strptime(args.date2,'%Y-%m-%d')
    except IndexError,ValueError:
        sys.stdout.write(helpstring+'\n')
        sys.exit()
    testCmorph(fyear,lyear,minmax_lons, minmax_lats, args.min_area)

