fom cluster import Cluster
import numpy
import netCDF4

"""

"""

class TimeClusters:

    def __init__(self):

        self.data = {}

    def addTimeStep(self, tIndex, otherCluster):
        """
        Add time step 
        @param tIndex present time index
        @param otherCluster will be added to this if its centre is inside one of the clusters
        OR DO WE WANT TO MERGE?
        """

        for cl in self.data[tIndex - 1]:
            if cl.containsCentreOf(otherCluster.getCentre()):
                cl.merge(otherCluster)


    def trackBackwards(self):
        pass

    def writeFile(self, filename):
        """
        Write the time clusters to file
        """
        pass

