from clustersystem import ClusterSystem
import numpy as np

"""
A class that handles cluster systems in time
"""

class TimeClusterSystem:

    def __init__(self,):

    	self.time_index = -1
    	self.time_clusters = {}


    def addTime(self, data, thresh_low, thresh_high):
    	"""
    	"""
    	new_cluster_system = ClusterSystem(data, thresh_low, thresh_high)
    	new_num_clusters = new_cluster_system.getNumClusters()
    	old_cluster_system = self.time_clusters.get(self.time_index - 1, [])
 		old_num_clusters = 0
 		if old_cluster_system:
 			old_num_clusters = old_cluster_system.getNumClusters()

    	self.time_index += 1
    	self.time_clusters[self.time_index] = new_cluster_system

    	# remove the clusters that don't overlap with clusters at time t - dt
    	for i in range(new_num_clusters):
    		cli = new_cluster_system.getCluster(i)
    		for j in range(old_num_clusters):
    			clj = old_cluster_system.getCluster(j)
    			if not cli.isCentreInsideOf(clj):
    				new_cluster_system.remove(i)

    def trackBackwards(self, time_index):
    	"""
    	"""
    	current_clusters = self.time_clusters[time_index]
    	current_num_clusters = current_clusters.getNumClusters()

    	previous_clusters = self.time_clusters.get(time_index - 1, [])
    	if not previous_clusters:
    		# we're at the beginning of time
    		return
    	previous_num_clusters = previous_clusters.getNumClusters()

    	for i in range(current_num_clusters):
    		cli = current_clusters[i]
    		for j in range(previous_clusters):
    			clj = previous_clusters[j]
    			if clj.isCentreInsideOf(cli):
    				# join the two clusters
    				self.cluster_connectivity[(time_index, i)] = self.cluster_connectivity.get((time_index, i), []) + TO DO 



###############################################################################
def test(filename):
    import netCDF4
    
    f = netCDF4.Dataset(filename, 'r')
    data = np.flipud(f.variables['CMORPH'][:, lat_slice, lon_slice])
    # TO DO 
    

if __name__ == '__main__':
    filename = 'Cmorph-2010_01_10.nc'
    test(filename)
