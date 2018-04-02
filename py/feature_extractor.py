'''
Created on March 28, 2018

@author: David Coppin
@institution: Department of Physics, University of Auckland

@description: A Class that extracts the clusters from the original data 
              using watershed algorithm from opencv. 
              Uses two different thresholds, the highest to estimate the convective cores,
              the lowest to keep the enveloppe of lower intensity around these cores.
              If the lower threshold is zero, we keep all the data
'''

import numpy
from scipy import ndimage
from skimage.morphology import watershed
from cluster import Cluster
import cv2

"""
A Class that extracts the clusters from the original data using watershed algorithm from opencv
"""

class FeatureExtractor:

    def __init__(self, data, thresh_low, thresh_high):
        """
        Extract clusters from an image data
        @param data
        @param thresh_low
        @param thresh_high
        @return list of clusters
        """
	    # remove data below minimum threshold
        ma_data = numpy.ma.masked_where(data <= thresh_low, data)
        # build black and white image with lower threshold to create borders for watershed
        tmp_data = ma_data.filled(fill_value=0)
        tmp_data[numpy.where(tmp_data !=0)] = 255
        bw_data = tmp_data.astype(numpy.uint8)
        border = cv2.dilate(bw_data, None, iterations=5)
        border -= cv2.erode(border, None)
    
        # remove data below minimum threshold
        ma_conv = numpy.ma.masked_where(data <= thresh_high, data)
	    # build black and white image with high threshold to serve as markers for watershed	
        tmp_conv = ma_conv.filled(fill_value=0)
        tmp_conv[numpy.where(tmp_conv !=0)] = 255
        bw_conv = tmp_conv.astype(numpy.uint8)
        markers = ndimage.label(bw_conv, structure=numpy.ones((3, 3)))[0]
	    # add border on image with high threshold to tell the watershed where it should fill in
        markers[border == 255] = 255
        # label each feature
        self.labels = watershed(-data, markers, mask=bw_data)
    

    def getClusters(self, min_ellipse_axis=1):
        """
        Get the features as clusters
        @param min_ellipse_cluster min ellipse axis size
        @return list of clusters, each cluster is a feature
        """
        res = []
        for idVal in range(1, self.labels.max()+1):
            iVals, jVals = numpy.where(self.labels == idVal)
            numVals = len(iVals)
            if numVals > 0:
            	# store as a set
                cells = {(iVals[i], jVals[i]) for i in range(len(iVals))}
                res.append(Cluster(cells, min_ellipse_axis))
        return res
    


#############################################################################################
'''
def testRectangle():
    ell = Ellipse({(i, 0) for i in range(3)}.union({(i, 1) for i in range(3)}))
    print(ell)

    pts = []
    # these points should be inside
    pt = numpy.array([1., 0.5])
    assert(ell.isPointInside(pt))
    pts.append(pt)

    pt = numpy.array([1.8, 0.5])
    assert(ell.isPointInside(pt))
    pts.append(pt)

    pt = numpy.array([1., 0.99])
    assert(ell.isPointInside(pt))
    pts.append(pt)

    # these points should be outside
    pt = numpy.array([1.82, 0.5])
    #assert(not ell.isPointInside(pt))
    pts.append(pt)

    pt = numpy.array([1., 1.01])
    #assert(not ell.isPointInside(pt))
    pts.append(pt)

    #ell.show(pts)


if __name__ == '__main__':
    testRectangle()
'''

