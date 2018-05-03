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

import numpy as np
from scipy import ndimage
from skimage.morphology import watershed
from cluster import Cluster
import cv2
from scipy.sparse import csr_matrix
import matplotlib.pyplot as mpl

"""
A Class that extracts the clusters from the original data using watershed algorithm from opencv
"""

class FeatureExtractor:

    def __init__(self, data, thresh_low, thresh_high, mask, frac):
        """
        Extract clusters from an image data
        @param data
        @param thresh_low
        @param thresh_high
        @param mask
        @param frac: overlap threshold for mask
        @return list of clusters
        """
        self.mask = mask
        self.frac = frac
	    # remove data below minimum threshold
        ma_data = np.ma.masked_where(data <= thresh_low, data)
        # build black and white image with lower threshold to create borders for watershed
        tmp_data = ma_data.filled(fill_value=0)
        tmp_data[np.where(tmp_data !=0)] = 255
        bw_data = tmp_data.astype(np.uint8)
        border = cv2.dilate(bw_data, None, iterations=5)
        border -= cv2.erode(border, None)
    
        # remove data below minimum threshold
        ma_conv = np.ma.masked_where(data <= thresh_high, data)
	    # build black and white image with high threshold to serve as markers for watershed	
        tmp_conv = ma_conv.filled(fill_value=0)
        tmp_conv[np.where(tmp_conv !=0)] = 255
        bw_conv = tmp_conv.astype(np.uint8)
        markers = ndimage.label(bw_conv, structure=np.ones((3, 3)))[0]
	    # add border on image with high threshold to tell the watershed where it should fill in
        markers[border == 255] = 255
        # label each feature
        self.labels = watershed(-data, markers, mask=bw_data)
        # remove clusters outside of mask
        self.labels = self.newremoveLargeScale()


    def getIndicesSparse(self, data):
        M = self.compute_M(data)
        return [np.unravel_index(row.data, data.shape) for row in M]


    def compute_M(self, data):
        cols = np.arange(data.size)
        return csr_matrix((cols, (data.ravel(), cols)),
                            shape=(data.max() + 1, data.size))


    def removeLargeScale(self):
        """
        Remove clusters outside of the mask
        """
        new_label = np.zeros_like(self.labels)
        nb = 1
        for label in np.unique(self.labels):
            # if the label is zero, we are examining the 'background'
            if label == 0:
                continue
            # get the indices of this feature
            inds = np.where(self.labels == label)
            # number of cells for this feature
            num_label_elems = len(inds[0])
            if self.mask[inds].sum() >= self.frac * num_label_elems:
                new_label[inds] = nb
                nb += 1
        return new_label

    def newremoveLargeScale(self):
        new_label = np.zeros_like(self.labels)
        nb = 1
        inds = self.getIndicesSparse(self.labels)
        # inds[0] is the background, we don't want it
#        inds = np.delete(inds2, (0), axis=0)
        for num in range(1, len(inds)):
            num_elems = len(inds[num])
#            print 'inds[nb]', inds[num]
            if self.mask[inds[num]].sum() >= self.frac * num_elems:
                new_label[inds[num]] = nb
                nb += 1
        return new_label


    def getClusters(self, min_ellipse_axis=1):
        """
        Get the features as clusters
        @param min_ellipse_cluster min ellipse axis size
        @return list of clusters, each cluster is a feature
        """
        res = []
        for idVal in range(1, self.labels.max()+1):
            iVals, jVals = np.where(self.labels == idVal)
            numVals = len(iVals)
            if numVals > 0:
            	# store as a set
                cells = {(iVals[i], jVals[i]) for i in range(len(iVals))}
                res.append(Cluster(cells, min_ellipse_axis))
        return res


    def newgetClusters(self, min_ellipse_axis=1):
        """
        Get the features as clusters
        @param min_ellipse_cluster min ellipse axis size
        @return list of clusters, each cluster is a feature
        """
        res = []
        inds = self.getIndicesSparse(self.labels)
        # inds[0] is the background, we don't want it
        for num in range(1, len(inds)):
            num_elems = len(inds[num])
            if num_elems > 0:
                # store as a set
                cells = {(inds[num][0][i], inds[num][1][i]) for i in range(len(inds[num][0]))}
                res.append(Cluster(cells, min_ellipse_axis))
        return res


#############################################################################################
'''
def testRectangle():
    ell = Ellipse({(i, 0) for i in range(3)}.union({(i, 1) for i in range(3)}))
    print(ell)

    pts = []
    # these points should be inside
    pt = np.array([1., 0.5])
    assert(ell.isPointInside(pt))
    pts.append(pt)

    pt = np.array([1.8, 0.5])
    assert(ell.isPointInside(pt))
    pts.append(pt)

    pt = np.array([1., 0.99])
    assert(ell.isPointInside(pt))
    pts.append(pt)

    # these points should be outside
    pt = np.array([1.82, 0.5])
    #assert(not ell.isPointInside(pt))
    pts.append(pt)

    pt = np.array([1., 1.01])
    #assert(not ell.isPointInside(pt))
    pts.append(pt)

    #ell.show(pts)


if __name__ == '__main__':
    testRectangle()
'''

