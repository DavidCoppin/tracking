'''
Created on March 28, 2018

@author: David Coppin
@institution: Department of Physics, University of Auckland

@description: A Class that creates the different coastal mask used in the algorithm
'''
import numpy as np
from scipy import ndimage
import cv2
from netCDF4 import Dataset as nc
import argparse
import matplotlib.pyplot as mpl

class CoastalMapping:

    def __init__(self, dataname, reso, lat_slice, lon_slice, szone, lzone, min_size, max_size):
        """
        Extract clusters from an image data
        @param dataname: land-sea mask
        @param reso
        @param lat_slice
        @param lon_slice
        @param szone: number of pixels to define area close to islands
        @param lzone: number of pixels to define large coastal area where clusters are tracked
        @param min_size: minimal size of island below which islands are removed from mask
        @param max_size: maximal size of island that should be filled
        @return mask of the two coastal areas
        """
        self.reso = reso
        self.min_size = min_size
        self.max_size = max_size
        print type(self.reso), self.reso
        if self.reso==8:
            print dataname
            slm = nc(dataname).variables['lsm'][lat_slice,lon_slice]
        else :
            print 'WTF reso', reso, self.reso
            slm_3d = nc(dataname).variables['lsmask'][:,:,:]
            len_lat = np.shape(slm_3d)[1]
            slm_short = slm_3d[:,len_lat/6:len_lat-len_lat/6,:] # remove first and last 30 degres
            slm = slm_short.squeeze()
        new_slm = np.flipud(slm)
        #Create the coastline and fill islands
        slm_fill = self.fillIslands(new_slm.copy())
        land_fill = (1-new_slm)+slm_fill
        land_fill[np.where(land_fill>=1)] = 1
        #Remove very small islands to speed up tracking
        slm_nosmall = self.eraseIslands(land_fill,new_slm)
        mask_coast = self.findCoastline(slm_nosmall,smooth_radius=2)
        mask = np.where(((slm_fill+mask_coast)/2.) >= 0.5, 1, 0)
        #Get the caostal area via inverse Box-counting
        self.lArea = self.createCoastalArea(mask_coast,lzone,1)
        self.sArea = self.createCoastalArea(mask_coast,szone,1)
        self.sArea[np.where(slm_fill>=1)] = 1


    def findCoastline(self,data,smooth_radius=2.5):
        """
        Finds coastlines in a slm-array
        smooth_radius = sigma-value for the canny algorithm
        data= the slm array
        """
        try:
            from skimage.feature import canny
        except ImportError:
            # old version
            from skimage.filter import canny
        from skimage.morphology import convex_hull_image
        return canny(data, sigma=smooth_radius).astype(np.int8)


    def fillIslands(self,slm):
        """
        Fills in islands up to a certain size
	    @param slm: land-sea mask
        @return array with big islands filled with 
        """
        slm = np.ma.masked_where(slm==0,slm).filled(255)
        slm = np.ma.masked_where(slm==1,slm).filled(0)
        tmp = slm.astype(np.uint8)
        result = cv2.findContours(tmp,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
        cnt = result[-2]
        hir = result[-1]
        hir = hir[0]
        #How many km2 representing one gridbox-pixel
        A = self.reso*self.reso*1.
        #size in grid points
        size = self.max_size/A
        for i,c in zip(hir,cnt):
            if cv2.contourArea(c) > size:
                cv2.drawContours(slm,[c],-1,0,-1)
        slm = np.ma.masked_where(slm==255,slm).filled(1)
        return slm.astype(np.int8)


    def eraseIslands(self,land,slm):
        """
        Removes small islands from land-sea mask
        @param land: array with islands filled and coastline.
        @param slm: array to get good coastline detection on side of domain
        @return new_mask: the new land-sea mask (without smaller islands)
        """
        Tmp = land.astype(np.uint8)
        result = cv2.findContours(Tmp,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)
        cnts = result[-2]
        tmpary = np.zeros(land.shape)
        new_mask = np.zeros(land.shape)
        for index,cnt in enumerate(cnts):
            cimg = np.zeros_like(land)
            cv2.drawContours(cimg, cnts, index, color=255, thickness=-1)
            ## Calculate the area covered by the contour km**2
            Area = self.reso**2*cv2.contourArea(cnt)
            if Area > self.min_size:
                cv2.drawContours(tmpary,[cnt],-1,1,-1)
        new_mask=1-tmpary
        new_mask[0,:] = slm[0,:]
        new_mask[:,0] = slm[:,0]
        new_mask[:,-1] = slm[:,-1]
        new_mask[-1,:] = slm[-1,:]
        return new_mask


    def createCoastalArea(self,data,radius,val):
        """
        Creates circles of certain radius around islands
        @param data: array of the islands
        @param zone: size of the zone around islands
        @param val: value of lands in mask
        @return new_mask: array with coastal areas
        """
        nY,nX = data.shape
        new_data = np.zeros((nY,nX))
        #Create circle mask around each coastal point
        for i in xrange(nY):
            for j in xrange(nX):
                if data[i,j] == val:
                    y,x = np.ogrid[-i:nY-i, -j:nX-j]
                    zone = x*x + y*y <= (2*radius+1)*(2*radius+1)
                    new_data[zone] = 1
        return new_data


#############################################################################################
def testCoastalArea(data, reso, minmax_lats, minmax_lons, szone, lzone, min_size, max_size):
    lon_slice = slice(minmax_lons[0], minmax_lons[1])
    lat_slice = slice(minmax_lats[0], minmax_lats[1])
    cm = CoastalMapping(data, np.int(reso), lat_slice, lon_slice, np.int(szone), np.int(lzone), 
                         np.int(min_size), np.int(max_size))
    mpl.subplot(1,2,1)
    mpl.contourf(np.flipud(cm.sArea))
    mpl.subplot(1,2,2)
    mpl.contourf(np.flipud(cm.lArea))
    mpl.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test Coastal Area')
    parser.add_argument('-w', action='store_true', help='Print copyright')
    parser.add_argument('-d', dest='data', default='Data/LSM/Cmorph_slm_8km.nc')
    parser.add_argument('-r', dest='reso', default='8', help='resolution of the land-sea mask in km')
    parser.add_argument('-lons', dest='lons', default='1200:2200', help='Min and max longitude indices LONMIN,LONMAX')
    parser.add_argument('-lats', dest='lats', default='50:650', help='Min and max latitude indices LATMIN,LATMAX')
    parser.add_argument('-sz', dest='szone', default='6', help='small distance to coast in pixels')
    parser.add_argument('-lz', dest='lzone', default='50', help='large distance to coast in pixels')
    parser.add_argument('-smin', dest='min_size', default='300', help='area below which islands are deleted in km2')
    parser.add_argument('-smax', dest='max_size', default='800000', help='max area of filled islands in km2')
    args = parser.parse_args()
    try:
        minmax_lons = np.array(args.lons.split(':')).astype(np.int)
    except:
        raise RuntimeError, 'Wrong specification of longitude bound indices, use -lons LONMIN:LONMAX'
    try:
        minmax_lats = np.array(args.lats.split(':')).astype(np.int)
    except:
        raise RuntimeError, 'Wrong specification of latitude bound indices, use -lats LATMIN:LATMAX'

    testCoastalArea(args.data, args.reso, minmax_lats, minmax_lons, args.szone, args.lzone, args.min_size, args.max_size)

