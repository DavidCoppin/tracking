Tracking algorithm for coastal precipitating systems
----------------------------------------------------

This repository contains source-code for the tracking of coastal rainfall precipitation.

For more details, see the scientific publication.
`Detection and analysis of high-resolution tropical coastal precipitation using a tracking algorithm: application to interactions between diurnal cycle and Madden-Julian Oscillation over the Maritime Continent`

It runs in two consecutive steps:
 * Tracking = tracks coastal systems and creates of pickles
 * Post-processing = generates netcdf daily outputs from the pickles created by Tracking

Note: Post-processing over the same time period as Tracking should be done once Tracking is finished. The last day of post-processing should not be after the last day of tracking.

This repository includes these modules and configuration files:
 * tracking.py - python module to run the tracking algorithm
 * config.cfg - configuration file used by tracking.py


Prerequisites
-------------

Requirements are:

 * python 2.7
 * scipy 0.19.1
 * netCDF4 1.3.1
 * scikit-image 0.12.3 (we found that 0.13.x produces different results and has a memory leak)
 * nc_time-axis 1.0.2
 * opencv 2.4.11 or opencv3
 * configparser
 * Cython

The code will currently not run with opencv 3.x.

All packages are open source and should be available via the package manager of
your OS.


Compiling
---------

The ellipse.pyx and time_connected_clusters.pyx modules must be compiled using Cython, to do this run:
```
./compile.sh
```
You **must** run this whenever you make changes in ellipse.pyx and/or time_connected_clusters.pyx.


Data and format
---------------
The code expects ``daily`` rainfall input files between 30S-30N in netcdf format. As long as it is sub-daily the temporal resolution doesn't matter. If you are not using the 8km spatial resolution or other latitudes, you will have to specify arguments lats and lons.

You will also need to supply a land-sea mask in netcdf format. This mask has to have the same spatial resolution as the rainfall estimates and is supposed to be a 2D-array (latlon).


Building
--------
Additional building should not be necessary


Usage
-----
The tracking algorithm has to be used in single threading mode.
Before running the code, you will have to edit the file config.cfg. This file provides all necessary parameters to run the code. The following variables are set:
 * data_path     = the directory where the precipitation data is stored
 * lsm_path      = the path to the land-sea mask
 * targetdir     = the directory where the pickles are stored after harvest
 * varname       = the name of the precipitation variable
 * units         = the units of the precipitation data
 * reso          = the spatial resolution of the data in km
 * min_prec      = the minimum value of precipitation used for the cluster detection
 * max_prec      = the minimum value of precipitation used to identify the convective core
 * szone         = the distance to the coast (in pixels) used for the coastal mask
 * lzone         = the distance to the coast (in pixels) used for the large-scale mask
 * frac_mask     = the minimum overlap with the large-scale mask needed to keep a cluster and track it
 * frac_ellipse  = the minimum overlap between two ellipses needed to merge them
 * min_axis      = the minimum length of the ellipse axes
 * min_size      = the minimum area (in km2) below which islands are deleted from the masks
 * max_size      = the maximum area of islands (in km2) filled by the coastal mask
 * max_cells     = the maximum size (in pixels) needed for a track to be considered large-scale and filtered (if t_life exceeds its threshold too)
 * t_life        = the minimum lifetime (in days) needed for a track to be considered large-scale and removed (if max_cells exceeds its threshold too)
 * t_life_lim    = same as t_life but for track that touch the latitudinal limits of the domain
 * frac_decrease = the decrease in size (between two time steps) over which a track is cut 
 * save          = to save time connected clusters for debugging

The module that reads the configuration is written in a way that you can add any variable you want to the config file. It simply has to have the following structure: key = value


---------------------
Running the algorithm
=====================

First part
----------
To launch the tracking (first part), you can simply run:
 `python tracking.py`

In that case, the default values defined as arguments in tracking.py are used.
To specify starting and ending dates, use -d1 and -d2. To run over specific regions, use lats and lons to define the latitude. The code runs faster over small regions. A prefix can be given to name files differently over different regions (with -suffix, used for other parts of the code too). The frequency to extract the tracks that have finished is defined by -harvest.
For these arguments, this gives:

` python tracking.py -d1 Start -d2 End -lats Lat_Start:Lat_End -lons Lon_Start:Lon_End -suffix \
                      SUFFIX -harvest Nb_time_steps `

where Start and End must have to be in YYYY-MM-DD format.

Note1: lats and lons needs pixels indices, not latitudes and longitudes
Note2: This algorithm uses CMORPH data for which latitude goes from 30S to 30N. Newest data may go from north to south.

Example :
 python tracking.py -d1 2011-10-01 -d2 2011-12-31 -lons 335:2717 -lats 0:827 -suffix io-cm \
                          -harvest 24 -restart_dir restart -restart_interval 183

Because the algorithm can take a long time to run, it is recommended to run it with -restart_dir and -restart_interval options. -restart_dir is the directory where the restart file will be kept. -restart_interval controls after how many days a restart file is saved. This way, if the algorithm stops for any reason or if you want to continue a previous tracking, the tracking algorithm will start with the tracks from the last day directly.


Second part
-----------
To run the post-processing, you can simply run:
` python write_output_pp.py `

In that case, the default values defined as arguments in write_output_pp.py are used.
Only the end date (after -d) is needed. The post-processing starts from the first day when pickles are available. The directory where the pickles needed are located is specified after -i. The post-processing can aggregate files with different prefixes. They need to be specified after -p. If the prefix of the pickles (-suffix in tracking) is png, use python write_output_pp.py -p png

Because the algorithm can take a long time to run, it is recommended to run it with -restart_dir, which is the directory where an info_pp.pkl file will be kept. This file contains the last filename used, the ongoing tracks at the end of the last day and the last ID used. This way, we can restart from the previous output file and know the tracks and ID used for the last day.

Example:
 python write_output_pp.py -d 2015_12_31 -s  -restart_dir restart_pp


Contributing
------------
We welcome all types of contributions, from blueprint designs to documentation, to testing or to deployment scripts.


Bugs
----
Bugs should be filed to **d.coppin@auckland.ac.nz**
