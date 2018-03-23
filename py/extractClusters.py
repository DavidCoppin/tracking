import numpy
from scipy import ndimage
from skimage.morphology import watershed
from cluster import Cluster
import cv2

"""
A Class that extracts the clusters from the original data using watershed algorithm from opencv
"""

class ExtractClusters:

    def __init__(self, data, thresh_min, thresh_max):
        """
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
        tmp_data[numpy.where(tmp_data !=0)] = 255
        bw_data = tmp_data.astype(numpy.uint8)
        border = cv2.dilate(bw_data, None, iterations=5)
        border -= cv2.erode(border, None)
    
        # remove data below minimum threshold
        ma_conv = numpy.ma.masked_where(data <= thresh_max, data)
	# building black and white image with high threshold to serve as markers for watershed	
        tmp_conv = ma_conv.filled(fill_value=0)
        tmp_conv[numpy.where(tmp_conv !=0)] = 255
        bw_conv = tmp_conv.astype(numpy.uint8)
        markers = ndimage.label(bw_conv, structure=numpy.ones((3, 3)))[0]
	# add border on image with high threshold to tell the watershed where it should fill in
        markers[border == 255] = 255
        # labels each feature
        self.labels = watershed(-data, markers, mask=bw_data)
    
    def extractPoints(self):
        # load into clusters
        res = []
        for idVal in range(1, self.labels.max()+1):
            iVals, jVals = numpy.where(self.labels == idVal)
            numVals = len(iVals)
            if numVals > 0:
                cells = [(iVals[i], jVals[i]) for i in range(len(iVals))]
                # store this cluster as a list with one element (so far). Each 
                # element will have its own ID
                res.append(Cluster(cells))
#		res.append(cells)    
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

