from time_connected_clusters import TimeConnectedClusters
import numpy as np
import netCDF4


"""
Write file with two variables: cprec (coastal precipitation) and nb (indexing of clusters) 
after the clusters have been transformed into arrays
"""

class OutputFile:

    def __init__(self):
        """
        Constructor
        """
        # precipitation data
        self.precip = []

        # time
        self.time = []


    def getPrecip(self, var, time):
        """
        Store precipitation data and append it
        @param var data
        @param time time index
        """
        if np.shape(self.precip)[0]==0:
            self.precip = var
            self.time = time
        else:
            self.precip = np.vstack([self.precip, var])
            self.time = np.append(self.time, time)


    def writeFile(self, filename, unit, lat, lon, i_minmax=[], j_minmax=[]):
        """
        Write data to netcdf file
        @param filename file name
        @param i_minmax min/max lat indices
        @param j_minmax min/max lon indices
        """
        f = netCDF4.Dataset(filename, 'w')

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

        # write the data
        i_index[:] = lat[:]
        j_index[:] = lon[:]
        t_index[:] = self.time[:]

        # get 3D array of clusters from TimeConnectedClusters
        data = TimeConnectedClusters().toArray(i_minmax=[], j_minmax=[])
        # now write all the data in one go
        mask=np.zeros((self.time, num_i, num_j))
        mask[np.where(data != 0)] = 1
        precip[:] = self.precip*mask 
        nb_var[:] = data

        f.close()

