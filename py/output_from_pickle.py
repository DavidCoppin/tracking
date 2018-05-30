'''
@description: A Class that selects the pickles needed for the output, extract tracks from
              the pickles and put them with a unique in a 3D array that will be used as
              output and as a mask on precipitation, and write the final output
'''

import cPickle
import numpy as np
import netCDF4
from netCDF4 import Dataset as nc
import bz2
import os
import gzip
import matplotlib.pyplot as mpl


def createTxt(filename, list):
    """
    @param filename: file name of txt file to be written
    @param list
    @return txt file
    """
    info = open(filename, 'w')
    for item in list:
        info.write("%s\n" % item)


def readTxt(filename):
    """
    @param filename: file name of the file to be read
    @return list of each line in the txt file
    """
    name = []
    with open(filename, 'r') as file:
        for line in file:
            # remove linebreak which is the last character of the string and add to list
            name.append(line[:-1])
    return name


class OutputFromPickle:

    def __init__(self, nb_day, lat, lon, inputdir, outputdir, list_prefix, dict_pickles, \
                  track_id, id):
        """
        Constructor
        """
        # precipitation data
        self.precip = []

        # clusters
        self.clusters = []

        # time step for the beginning
        self.ini = nb_day*48

        # time step for the beginning
        self.end = (nb_day+1)*48

        # name of the pickle files to be read
        self.filenames = []

        # dictionary with name of pickle file as key, value is a np.array with same length as
        # number of tracks from pickle, set at 0 by default and replace by track Id when track
        # over several output files
        self.track_id = track_id

        # latitude range of the final output
        self.lat = lat

        # longitude range of the final output
        self.lon = lon

        # id of the cluster
        self.id = id

        # path where the pickles are stored
        self.inputdir = inputdir

        # path where the output is stored
        self.outputdir = outputdir

        # list of prefix of the different regions to include in the output
        self.list_prefix = list_prefix

        # dictionary with pickle name as key and file as value
        self.dict_pickles = dict_pickles


    def selectPickles(self):
        """
        @return list with all the pickles needed for a day
        """
        for n in range(len(self.list_prefix)):
            test = [i for i in os.listdir(self.inputdir)]
            files = [i for i in os.listdir(self.inputdir) if \
                      os.path.isfile(os.path.join(self.inputdir,i)) \
                      and i.startswith(self.list_prefix[n])]
            for nb in range(len(files)):
                num = [int(s) for s in files[nb].split('_') if s.isdigit()]
                if num[0] <= self.end :
                    self.filenames.append(self.inputdir+files[nb])
        return self.filenames


    def extractTracks(self, files):
        """
        Extract tracks that correspond to the time of the file and gives them
        a unique id
        @param files: all the pickles files for this day
        """
        self.clusters = np.zeros((self.end-self.ini, len(self.lat), len(self.lon)))
        tot_lon = len(self.lon)
        for i in files:
            lat_min, lat_max, lon_min, lon_max = self.getLatLon(i)

            # Check if file already unzipped
            if i not in self.dict_pickles:
                with gzip.GzipFile(i) as gzf:
                    tracks = cPickle.load(gzf)
                    self.dict_pickles[i]=tracks
            else:
                tracks = self.dict_pickles[i]

#            with gzip.GzipFile(i) as gzf:
#                tracks = cPickle.load(gzf)

            # Set default track Id = 0 for all tracks if first time that file is read
            if len(self.track_id) == 0:
                self.track_id = {i: np.zeros(len(tracks))}
            if i not in self.track_id :
                print 'i not in self.track_id', i
                self.setTrackId(i, len(tracks))
            print 'i', i
            for nb in range(len(tracks)):
                keys = sorted(tracks[nb].keys())

                # Check if track has an Id different from 0
                if self.track_id.get(i)[nb] > 0 and ((keys[0] >= self.ini and keys[0] <= self.end) \
                        or (keys[-1] >= self.ini and keys[-1] <= self.end)):
                    new_id = self.track_id.get(i)[nb]
                else :
                    new_id = self.id + 1
                    self.id = self.id + 1

                # Fill in clusters with new_id
                for k in keys:
                    if k >= self.ini and k < self.end:
                        for cl in tracks[nb][k]:
                            i_index, j_index, mat = cl.toArray()
                            if int(lon_min) < int(lon_max):
                                if len(self.list_prefix)==1:
                                    self.clusters[k-self.ini, i_index[0]:i_index[-1]\
                                               +1,j_index[0]:j_index[-1]\
                                               +1][np.where(mat==1)]= new_id
                                else:
                                    self.clusters[k-self.ini, i_index[0]+int(lat_min):i_index[-1]\
                                               +int(lat_min)+1,j_index[0]+int(lon_min):j_index[-1]\
                                               +int(lon_min)+1][np.where(mat==1)]= new_id

                            # Special case for region on both sides of 0 degree longitude
                            else :
                                len_west = tot_lon - int(lon_min)

                                # Only positive longitudes
                                if j_index[0] > len_west :
                                    self.clusters[k-self.ini, i_index[0]+int(lat_min):i_index[-1]\
                                                   +int(lat_min)+1,j_index[0]-len_west:j_index[-1]\
                                                   -len_west+1][np.where(mat==1)]= new_id

                                # Only negative longitudes
                                elif j_index[-1] < len_west :
                                    self.clusters[k-self.ini, i_index[0]+int(lat_min):i_index[-1]\
                                                   +int(lat_min)+1,j_index[0]+int(lon_min):j_index[-1]\
                                                   +int(lon_min)+1][np.where(mat==1)]= new_id

                                # Both
                                else :
                                    mat_west = mat[:,0:(len_west-j_index[0])]
                                    mat_east = mat[:,(len_west-j_index[0]):]
                                    self.clusters[k-self.ini, i_index[0]+int(lat_min):i_index[-1]\
                                                   +int(lat_min)+1,0:len(j_index)\
                                                   -(len_west-j_index[0])+1]\
                                                   [np.where(mat_east==1)]= new_id
                                    self.clusters[k-self.ini, i_index[0]+int(lat_min):i_index[-1]\
                                                   +int(lat_min)+1,j_index[0]+int(lon_min):]\
                                                   [np.where(mat_west==1)]= new_id

                            # Replace track ID kept for next output file if track goes further
                            # than end of this day
                            if keys[-1] >= self.end-self.ini:
                                self.track_id[i][nb] = new_id

                # if id not used, do not waste id
                if keys[-1] < self.ini or keys[0] >= self.end:
                    self.id = self.id -1
                    new_id = self.id - 1
                print 'new_id, self.id', new_id, self.id


    def getLatLon(self, file):
        """
        @param file: pickle file
        @return latitude and longitude to place area into global file
        """
        for i in self.list_prefix:
            if i in file:
                lat_min, lat_max, lon_min, lon_max = readTxt(self.inputdir+\
                          'lat-lon_'+str(i)+'.txt')
            else:
                pass
        return lat_min, lat_max, lon_min, lon_max


    def setTrackId(self, filename, nb_tracks):
        """
        If first time reading pickle file, add it to dictionary and set track id to 0 for
        every track
        @param filename: pickle file
        @param nb_tracks: number of tracks in pickle file
        """
        self.track_id[filename] = np.zeros(nb_tracks)


    def deletePickles(self):
        """
        Delete pickle once all its tracks are finished
        """
        for nb in range(len(self.filenames)):
            num = [int(s) for s in self.filenames[nb].split('_') if s.isdigit()]
            if num[1] < self.end :
                print 'should delete self.filenames[nb]', self.filenames[nb]
                os.remove(self.filenames[nb])
                self.deleteTrackId(self.filenames[nb])
                del self.dict_pickles[self.filenames[nb]]
#                print 'len(self.dict_pickles)', len(self.dict_pickles)
#                print 'self.dict_pickles.keys()', self.dict_pickles.keys()


    def deleteTrackId(self, filename):
        """
        Delete key/value associated to a filename in track_id
        """
        self.track_id.pop(filename)


    def writeFile(self, suffix, old_filename, list_lat_lon):
        """
        Write data to netcdf file
        @param suffix: suffix for output
        @param old_filename: name of cmorph file corresponding
        """
        new_filename = str(self.outputdir)+'tracking'+str(old_filename[-18:-7])+'_'\
                        +str(suffix)+'.nc'
        f = netCDF4.Dataset(new_filename, 'w')

        # create dimensions
        num_i = len(self.lat)
        f.createDimension('lat', size=num_i)
        num_j = len(self.lon)
        f.createDimension('lon', size=num_j)

        # infinite dimension
        time = f.createDimension('time', size=None)

        lat_min, lat_max, lon_min, lon_max = list_lat_lon

        # read data needed
        zipfile = bz2.BZ2File(old_filename)
        data_unzip = zipfile.read()
        filename = old_filename[:-4]
        open(filename, 'wb').write(data_unzip)
        try:
            ori = nc(filename)
        except RuntimeError:
            ori = nc(filename.replace('-','_'))
        if max(list_lat_lon)>0:
            var = ori.variables["CMORPH"][:,lat_min:lat_max,lon_min:lon_max]
        else:
            var = ori.variables["CMORPH"][:,:,:]
        tint = ori.variables["time"][:]
        unit = ori.variables["time"].units
        os.remove(filename)

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

        nb_var = f.createVariable('nb', 'i4', ('time', 'lat', 'lon'), zlib=True)
        f.variables['nb'].gridtype='lonlat'
        f.variables['nb'].code=0
        f.variables['nb'].long_name='identification number of the clusters'
        f.variables['nb'].standard_name= 'number of the clusters'
        f.variables['nb'].short_name='nb'
        f.variables['nb'].units='no unit'

        i_index[:] = self.lat[:]
        j_index[:] = self.lon[:]
        t_index[:] = tint[:]

        # Write data for whole day
        mask=np.zeros((self.end-self.ini, num_i, num_j))
        mask[np.where(self.clusters != 0)] = 1
        precip[:] = var * mask
        nb_var[:] = self.clusters
        unique = np.unique(self.clusters)
        del var, mask, data_unzip
        f.close()


