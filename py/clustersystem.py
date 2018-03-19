from cluster import Cluster
import numpy as np
from scipy import ndimage
from skimage.morphology import watershed
import cv2

"""
A class that works on a system of clusters, all existing at a given time index
"""

class ClusterSystem:

    def __init__(self, data, thresh_min, thresh_max):

        ma_data = np.ma.masked_where(data <= thresh_min, data).astype(np.uint8)
        
        # building threshold
        tmp_data = ma_data.filled(0)
        tmp_data[np.where(tmp_data!=0)] = 255
        bw_data = tmp_data.astype(np.uint8)
        
        # building markers and borders
        ma_conv = np.ma.masked_where(data <= thresh_max, data)
        tmp_conv = ma_conv.filled(0)
        tmp_conv[np.where(tmp_conv!=0)] = 255
        bw_conv = tmp_conv.astype(np.uint8)
        markers= ndimage.label(bw_conv, structure=np.ones((3, 3)))[0]
        border = cv2.dilate(bw_data, None, iterations=5)
        border -= cv2.erode(border, None)
        markers[border == 255] = 255

        # watershed
        ws = watershed(-data, markers, mask=bw_data)

        # load into cluster
        self.clusters = []
        for idVal in range(1, ws.max()):
            iVals, jVals = np.where(ws == idVal)
            numVals = len(iVals)
            if numVals > 0:
                cells = [(iVals[i], jVals[i]) for i in range(len(iVals))]
                self.clusters.append(Cluster(cells))

    def removeLargeScale(self):
        pass

    def mergeClusters(self):
        """
        Merge the overlapping clusters, ie clusters whose centres are inside each other's ellipse
        """
        numClusters = len(self.clusters)
        clusterIndexToRemove = []
        for i in range(numClusters):
            cli = self.clusters[i]
            for j in range(i + 1, numClusters):
                clj = self.clusters[j]
                if cli.isCentreInsideOf(clj) and clj.isCentreInsideOf(cli):
                    cli.merge(clj)
                    clusterIndexToRemove.append(j)
        clusterIndexToRemove.reverse()
        for j in clusterIndexToRemove:
            del self.clusters[j]


    def checkCoastal(self):
        pass



    def removeNonCoastal(self):
        pass



###############################################################################
def test(filename, t, lat_slice, lon_slice):
    import netCDF4
    
    f = netCDF4.Dataset(filename, 'r')
    data = np.flipud(f.variables['CMORPH'][t, lat_slice, lon_slice])
    cs = ClusterSystem(data, thresh_min=0, thresh_max=0.8)
    cs.mergeClusters()
    

if __name__ == '__main__':
    filename = 'Cmorph-2010_01_10.nc'
    test(filename, 0, slice(50, 350), slice(1200, 1700))
