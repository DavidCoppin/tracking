from netCDF4 import Dataset as nc
import argparse
import numpy
import matplotlib
import time
from time_connected_clusters import TimeConnectedClusters
from cluster import Cluster
from scipy import ndimage
from skimage.morphology import watershed
import cv2
import sys,os,string
import bz2
from datetime import datetime,timedelta as td

def extractClusters(data, thresh_min, thresh_max):
    """.
    Extract clusters from an image data
    @param data
    @param thresh_min
    @param thresh_max
    @return list of clusters
    """
    # remove data below minimum threshold
    ma_data = numpy.ma.masked_where(data <= thresh_min, data)
    # building black and white image with lower threshold to create borders for watershed
    tmp_data = ma_data.filled(fill_value=0)
    tmp_data[numpy.where(tmp_data !=0)] = 255.
    bw_data = tmp_data.astype(numpy.uint8)
    border = cv2.dilate(bw_data, None, iterations=5)
    border -= cv2.erode(border, None)

    # remove data below minimum threshold
    ma_conv = numpy.ma.masked_where(data <= thresh_max, data)
    # building black and white image with high threshold to serve as markers for watershed..
    tmp_conv = ma_conv.filled(fill_value=0)
    tmp_conv[numpy.where(tmp_conv !=0)] = 255.
    bw_conv = tmp_conv.astype(numpy.uint8)
    markers = ndimage.label(bw_conv, structure=numpy.ones((3, 3)))[0]
    # add border on image with high threshold to tell the watershed where it should fill in
    markers[border == 255] = 255.
    # labels each feature
    labels = watershed(-data, markers, mask=bw_data)

    # load into clusters
    res = []
    for idVal in range(1, labels.max()+1):
        iVals, jVals = numpy.where(labels == idVal)
        numVals = len(iVals)
        if numVals > 0:
            cells = [(iVals[i], jVals[i]) for i in range(len(iVals))]
            # store this cluster as a list with one element (so far). Each.
            # element will have its own ID
            res.append(Cluster(cells))
    return res


def testCmorph(fyear,lyear,minmax_lons, minmax_lats):
    """
    Checking that we can create a time connected cluster from image
    """
    lon_slice = slice(minmax_lons[0], minmax_lons[1])
    lat_slice = slice(minmax_lats[0], minmax_lats[1])
    llat = minmax_lats[1] - minmax_lats[0]
    llon = minmax_lons[1] - minmax_lons[0]
    tcc = TimeConnectedClusters()
    delta = lyear - fyear
    dates = [fyear + td(days=i) for i in xrange(delta.days + 1)]
    print 'dates', dates
    for nb_day in xrange(len(dates)):
        date=dates[nb_day]
        filename=os.path.join('/home/dcop696/clusters/py/Data/CMORPH/Cmorph-'\
               +str(date.year)+'_'+str(date.month).zfill(2)+'_'\
               +str(date.day).zfill(2)+'.nc.bz2')
        #Open the file
        filename = filename.replace('--','-').replace('__','_')
        print 'filename', filename
        zipfile = bz2.BZ2File(filename)
        data_unzip = zipfile.read()
        newfilename = filename[:-4]
        open(newfilename, 'wb').write(data_unzip)
        try:
            f=nc(newfilename)
        except RuntimeError:
            f=nc(newfilename.replace('-','_'))
    for t in xrange(48) :
        print 'nb_day, t', nb_day, t
        data = f.variables["CMORPH"][t,lat_slice, lon_slice]
        clusters = extractClusters(numpy.flipud(data), thresh_min=0., thresh_max=2.5)
        print clusters
        tcc.addTime(clusters)

    print tcc

    # write to file
    tcc.writeFile('cmorph.nc', i_minmax=(0, len(lats)), j_minmax=(0, len(lons)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tracking.')
    parser.add_argument('-w', action='store_true', help='Print copyright')
    parser.add_argument('-d1', dest='date1', default='2010-01-12', help='Start date YYYY-MM-DD')
    parser.add_argument('-d2', dest='date2', default='2010-01-13', help='End date YYYY-MM-DD')
    parser.add_argument('-lons', dest='lons', default='1450:1700', help='Min and max longitude indices LONMIN,LONMAX')
    parser.add_argument('-lats', dest='lats', default='300:550', help='Min and max latitude indices LATMIN,LATMAX')
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
    testCmorph(fyear,lyear,minmax_lons, minmax_lats)

