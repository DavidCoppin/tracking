from netCDF4 import Dataset as nc
import argparse
import numpy
import matplotlib.pyplot as mpl
import time
from time_connected_clusters import TimeConnectedClusters
from feature_extractor import FeatureExtractor
from cluster import Cluster
from coastal_mapping import CoastalMapping
#from scipy import ndimage
#from skimage.morphology import watershed
#import cv2
import sys,os,string
import bz2
from datetime import datetime,timedelta as td

def testCmorph(lsm, fyear, lyear, minmax_lons, minmax_lats, reso, min_ellipse_axis, frac_ellipse, \
                suffix, szone, lzone, min_size, max_size, save=False):
    """
    Checking that we can create a time connected cluster from image
    """
    lon_slice = slice(minmax_lons[0], minmax_lons[1])
    lat_slice = slice(minmax_lats[0], minmax_lats[1])
    # Get the two coastal masks
    cm = CoastalMapping(lsm, numpy.int(reso), lat_slice, lon_slice, numpy.int(szone), \
                         numpy.int(lzone), numpy.int(min_size), numpy.int(max_size))
    llat = minmax_lats[1] - minmax_lats[0]
    llon = minmax_lons[1] - minmax_lons[0]
    tcc = TimeConnectedClusters()
    delta = lyear - fyear
    dates = [fyear + td(days=i) for i in xrange(delta.days + 1)]
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
        for t in xrange(48) :
            print 'nb_day, t', nb_day, t
            data = f.variables["CMORPH"][t, lat_slice, lon_slice]
            # Extract clusters with watershed
            clusters = FeatureExtractor(data, thresh_low=0., thresh_high=2.5).getClusters(min_ellipse_axis)
            new_clusters = clusters
            n=0
            for cl in clusters:
                # Apply large mask to remove clusters too far away from islands
#                print 'cl', cl
#                print 'cl.cells', cl.cells, len(cl.cells)
#                print 'Cluster(cm.lArea).cells', Cluster(cm.lArea).cells
                if cl.isClusterInsideOf(Cluster(cm.sArea), 0.9):
                    pass
                else:
                    new_clusters.remove(cl)
                    n=n+1
                    print 'n', n
#                   print 'cl removed', cl
            print 'len(clusters), len(new_clusters)', len(clusters), len(new_clusters)
            test_clusters = numpy.zeros((numpy.shape(data)[0], numpy.shape(data)[1]))
            test_new_clusters = numpy.zeros((numpy.shape(data)[0], numpy.shape(data)[1]))
            test_lArea = numpy.zeros((numpy.shape(data)[0], numpy.shape(data)[1]))
            test_sArea = numpy.zeros((numpy.shape(data)[0], numpy.shape(data)[1]))
            for cl in clusters:
                test_clusters[numpy.argwhere(cl.cells)]=1
            for new_cl in new_clusters:
                test_clusters[numpy.argwhere(new_cl.cells)]=1
#            test_lArea[numpy.where(Cluster(cm.lArea).cells)]=1
#            test_sArea[numpy.where(Cluster(cm.sArea).cells)]=1
            print 'numpy.max(test_clusters), numpy.max(test_new_clusters)', numpy.max(test_clusters), numpy.max(test_new_clusters)
            mpl.subplot(1,2,1)
            mpl.contourf(numpy.flipud(test_clusters))
#            mpl.contour(numpy.flipud(test_lArea))
            mpl.subplot(1,2,2)
            mpl.contourf(numpy.flipud(test_new_clusters))
#            mpl.contour(numpy.flipud(test_lArea))
            mpl.show()
            sys.exit()
#            print clusters
            tcc.addTime(new_clusters,frac_ellipse)
        os.remove(newfilename)
    # write to file
    lat = f.variables['lat'][minmax_lats[0]:minmax_lats[1]]
    lon = f.variables['lon'][minmax_lons[0]:minmax_lons[1]]
    f.close()
    tcc.writeFile('cmorph.nc_'+str(suffix), i_minmax=(0, len(lat)), j_minmax=(0, len(lon)))
    if save:
        tcc.save('cmorph.pckl_'+str(suffix))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tracking.')
    parser.add_argument('-w', action='store_true', help='Print copyright')
    parser.add_argument('-save', dest='save', action='store_true', help='Save time connected clusters \
                           object for debugging')
    parser.add_argument('-d1', dest='date1', default='2010-01-12', help='Start date YYYY-MM-DD')
    parser.add_argument('-d2', dest='date2', default='2010-01-13', help='End date YYYY-MM-DD')
    parser.add_argument('-lons', dest='lons', default='1450:1700', help='Min and max longitude \
                           indices LONMIN,LONMAX')
    parser.add_argument('-lats', dest='lats', default='300:550', help='Min and max latitude \
                           indices LATMIN,LATMAX')
    parser.add_argument('-min_axis', dest='min_axis', type=float, default=6, help='Min ellipse \
                           axis in pixels')
    parser.add_argument('-lsm', dest='lsm', default='Data/LSM/Cmorph_slm_8km.nc', help='path to \
                           land-sea data')
    parser.add_argument('-reso', dest='reso', default='8', help='resolution of the dataset in km')
    parser.add_argument('-frac_ellipse', dest='frac_ellipse', type=float, default=0.8, \
                           help='Threshold to merge overlapping ellipses')
    parser.add_argument('-suffix', dest='suffix', default='', help='suffix for output')
    parser.add_argument('-sz', dest='szone', default='6', help='small distance to coast in pixels')
    parser.add_argument('-lz', dest='lzone', default='50', help='large distance to coast in pixels')
    parser.add_argument('-smin', dest='min_size', default='300', help='area below which islands \
                           are deleted in km2')
    parser.add_argument('-smax', dest='max_size', default='800000', help='max area of filled \
                           islands in km2')
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
    testCmorph(args.lsm, fyear,lyear,minmax_lons, minmax_lats, args.reso, args.min_axis, \
                args.frac_ellipse, args.suffix, args.szone, args.lzone, args.min_size, \
                args.max_size, args.save)

