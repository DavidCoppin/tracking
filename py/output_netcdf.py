import cPickle
import numpy as np
import netCDF4
from netCDF4 import Dataset as nc
import bz2
import os
import matplotlib.pyplot as mpl

"""
Write file with two variables: cprec (coastal precipitation) and nb (indexing of clusters)
after the clusters have been transformed into arrays
"""

class OutputNetcdf:

    def __init__(self, nb_day, lat, lon, id):#, data):
        """
        Constructor
        """
        # precipitation data
        self.precip = []

        # clusters
        self.clusters = []

        # time
        self.time = []

        # time step for the beginning
        self.ini = nb_day*48

        # time step for the beginning
        self.end = (nb_day+1)*48

        # keep time connected clusters
#        self.tcc = data

        self.filenames = []

        # Dictionary with name of pickle file as key, value is a np.array with same length as 
        # number of tracks from pickle, set at 0 by default and replace by track Id when track
        # over several output files
        self.track_id = {}

        self.lat = lat

        self.lon = lon

        self.id = id

    def getTime(self, time):
        """
        Store precipitation data and append it
        @param var data
        @param time time index
        """
        self.time = np.append(self.time, time)


    def selectPickles(self, prefix):
        """
        Gather files from the same regions
        """
        files = [i for i in os.listdir('.') if os.path.isfile(os.path.join('.',i)) \
                  and i.startswith(prefix)]
        print 'files', files, len(files)
        for nb in range(len(files)):
            num = [int(s) for s in files[nb].split('_') if s.isdigit()]
            if num[0] <= self.end :
                print 'files[nb]', files[nb]
                self.filenames.append(files[nb])
        return self.filenames


    def extractTracks(self, files):
        """
        Extract tracks that correspond to the time of the file
        @param files: all the pickles files that will considered for this day
        """
        self.clusters = np.zeros((self.end-self.ini, len(self.lat), len(self.lon)))
#        self.clusters = np.zeros((self.end-self.ini, nb_lat, nb_lon))
        test = []
        for i in files:
           print i
           f = open(i)
           tracks = cPickle.load(f)
           # Set default track Id = 0 for all tracks if first time that file is read
           if len(self.track_id) == 0:
               self.track_id = {i: np.zeros(len(tracks))}
           if i not in self.track_id :
               self.setTrackId(i, len(tracks))
           for nb in range(len(tracks)):
               keys = tracks[nb].keys()
               # Check if track has an Id different from 0
               if self.track_id.get(i)[nb] > 0 :
                   new_id = self.track_id.get(i)[nb]
               else :
                   new_id = self.id + 1
                   self.id = self.id + 1
               # Fill in clusters with new_id where the track is
               for k in keys:
                   if k < self.end:
                       for cl in tracks[nb][k]:
                           i_index, j_index, mat = cl.toArray()
                           self.clusters[k-self.ini, i_index[0]:i_index[-1]+1, j_index[0]:\
                                          j_index[-1]+1][np.where(mat==1)]= new_id
               # Replace track ID kept for next output file if track goes further than future output
               if keys[-1] >= self.end-self.ini:
                   self.track_id[i][nb] = new_id
#               print 'i, nb, self.track_id[i][nb]', i, nb, self.track_id[i][nb]
#           mpl.contourf(self.clusters[4, :, :])
#           mpl.show()


    def setTrackId(self, filename, nb_tracks):
        """
        If first time reading pickle file, add it to dictionary and set track id to 0 for 
        every track
        """
        self.track_id[filename] = np.zeros(nb_tracks)
#        self.track_id.append({filename: np.zeros(nb_tracks)})


    def deletePickles(self):
        """
        Delete pickle once they are all used
        """
        for nb in range(len(self.filenames)):
            num = [int(s) for s in self.filenames[nb].split('_') if s.isdigit()]
            if num[1] <= self.ini :
                print 'self.filenames[nb]', self.filenames[nb]
                os.remove(self.filenames[nb])
                # Delete track_id associated with pickle
                deleteTrackId(self.filenames[nb])


    def deleteTrackId(self, filename):
        """
        Delete key/value associated to a filename in track_id
        """
        self.track_id.pop(filename)


    def writeFile(self, suffix, old_filename, unit, lat_slice, \
                       lon_slice):
        """
        Write data to netcdf file
        @param new_filename new file name
        @param latitude
        @param longitude
        @param ini and end: indicates the timesteps when the data needs to be written
        """
        new_filename = 'tracking'+str(old_filename[-18:-7])+'_'+str(suffix)
        f = netCDF4.Dataset(new_filename, 'w')

        # create dimensions
        num_i = len(self.lat)
        f.createDimension('lat', size=num_i)

        num_j = len(self.lon)
        f.createDimension('lon', size=num_j)

        # infinite dimension
        time = f.createDimension('time', size=None)

        # create variables

        i_index = f.createVariable('lat', 'f', ('lat',) ,zlib=True,complevel=9,\
                     least_significant_digit=4)
        f.variables['lat'].units='degrees_north'
        f.variables['lat'].standard_name='latitude'
        f.variables['lat'].axis='Y'
        f.variables['lat'].long_name='latitude'

        j_index = f.createVariable('lon', 'f', ('lon',), zlib=True,complevel=9,\
                     least_significant_digit=4)
        f.variables['lon'].units='degrees_east'
        f.variables['lon'].standard_name='longitude'
        f.variables['lon'].axis='X'
        f.variables['lon'].long_name='longitude'

        t_index = f.createVariable('time', 'f', ('time',), zlib=True,complevel=9,\
                     least_significant_digit=4)
        f.variables['time'].axis='T'
        f.variables['time'].calendar='standard'
        f.variables['time'].units=unit

        precip = f.createVariable('cprec','f',('time','lat','lon'),zlib=True,complevel=9,\
                    least_significant_digit=4)
        f.variables['cprec'].gridtype='lonlat'
        f.variables['cprec'].code=999
        f.variables['cprec'].long_name='detected coastal precipitation'
        f.variables['cprec'].standard_name= 'detected_precipitation'
        f.variables['cprec'].short_name='cprec'
        f.variables['cprec'].units='mm/h'

        nb_var = f.createVariable('nb', 'i4', ('time', 'lat', 'lon'))
        f.variables['nb'].gridtype='lonlat'
        f.variables['nb'].code=0
        f.variables['nb'].long_name='identification number of the clusters'
        f.variables['nb'].standard_name= 'number of the clusters'
        f.variables['nb'].short_name='nb'
        f.variables['nb'].units='no unit'

        # write the output
        zipfile = bz2.BZ2File(old_filename)
        data_unzip = zipfile.read()
        filename = old_filename[:-4]
        open(filename, 'wb').write(data_unzip)
        try:
            ori = nc(filename)
        except RuntimeError:
            ori = nc(filename.replace('-','_'))
        if lon_slice.start < lon_slice.stop:
            var = ori.variables["CMORPH"][:, lat_slice, lon_slice]
        else:
            var1 = ori.variables["CMORPH"][:, lat_slice, lon_slice.start:]
            var2 = ori.variables["CMORPH"][:, lat_slice, :lon_slice.stop]
            var = np.concatenate((var1, var2), axis=2)
            del var1, var2
        os.remove(filename)

        i_index[:] = self.lat[:]
        j_index[:] = self.lon[:]
        t_index[:] = self.time[self.ini:self.end]
        # now write all the data in one go
        mask=np.zeros((self.end-self.ini, num_i, num_j))
        mask[np.where(self.clusters != 0)] = 1
        precip[:] = var * mask
        nb_var[:] = self.clusters
        del var, mask, data_unzip
        f.close()


###################
def testWrite():
        id = 0
        ini = 0
        end = 40
        on = OutputNetcdf(ini, end)
        files = on.selectPickles('png')
        print 'files', files
        on.extractTracks(files, 300, 500, id)
#        on.deletePickles()
#        on.writeFile(str(suffix), list_filename[nb_day], lat, lon, \
#                              unit, lat_slice, lon_slice)


if __name__ == '__main__':
    testWrite()
