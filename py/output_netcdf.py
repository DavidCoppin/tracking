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

    def __init__(self):#, data):
        """
        Constructor
        """
        # precipitation data
        self.precip = []

        # clusters
        self.clusters = []

        # time
        self.time = []

        # keep time connected clusters
#        self.tcc = data

        self.filenames = []

    def getTime(self, time):
        """
        Store precipitation data and append it
        @param var data
        @param time time index
        """
        self.time = np.append(self.time, time)


    def selectPickles(self, prefix, t_index):
        """
        Gather files from the same regions
        """
        files = [i for i in os.listdir('.') if os.path.isfile(os.path.join('.',i)) \
                  and i.startswith(prefix)]
        print 'files', files, len(files)
        for nb in range(len(files)):
            num = [int(s) for s in files[nb].split('_') if s.isdigit()]
            if num[0] <= t_index :
                print 'files[nb]', files[nb]
                self.filenames.append(files[nb])
        return self.filenames


    def extractTracks(self, files, lat, lon, id):
        """
        Extract tracks that correspond to the time of the file
        """
        self.id = id
#        self.clusters = np.zeros((48, len(lat), len(lon)))
        self.clusters = np.zeros((48, lat, lon))
        for i in files:
           print i
           f = open(i)
           tracks = cPickle.load(f)
           for nb in range(len(tracks)):
               keys = tracks[nb].keys()
               for k in keys:
                   for cl in tracks[nb][k]:
                       i_index, j_index, mat = cl.toArray()
                       self.clusters[k, i_index[0]:i_index[-1]+1, j_index[0]:j_index[-1]+1]\
                                     [np.where(mat==1)]= self.id+1
               self.id = self.id + 1
           mpl.contourf(self.clusters[4, :, :])
           mpl.show()


    def deletePickles(self, t_index):
        """
        Delete pickle once they are all used
        """
        for nb in range(len(self.filenames)):
            num = [int(s) for s in self.filenames[nb].split('_') if s.isdigit()]
            if num[1] <= t_index :
                print 'self.filenames[nb]', self.filenames[nb]
                os.remove(self.filenames[nb])


    def writeFile(self, suffix, old_filename, data, lat, lon, ini, end ,unit, lat_slice, \
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
        num_i = len(lat)
        f.createDimension('lat', size=num_i)

        num_j = len(lon)
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

        i_index[:] = lat[:]
        j_index[:] = lon[:]
        t_index[:] = self.time[ini:end]
        # now write all the data in one go
        mask=np.zeros((end-ini, num_i, num_j))
        mask[np.where(data != 0)] = 1
        precip[:] = var * mask
        nb_var[:] = data
        del var, mask, data_unzip
        f.close()


###################
def testWrite():
        id = 0
        on = OutputNetcdf()
        files = on.selectPickles('png', 5)
        print 'files', files
        on.extractTracks(files, 300, 500, id)
#        on.deletePickles(48)


if __name__ == '__main__':
    testWrite()
