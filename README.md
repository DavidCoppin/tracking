# Optimisation of tracking algorithm for precipitation systems

## How to run tracking

Requirements are:

 * python 2.7
 * scipy 0.19.1
 * netCDF4 1.3.1
 * scikit-image 0.12.3
 * nc_time-axis 1.0.2
 * opencv 2.4.11 or opencv3
 * configparser
 * Cython

The code will currently not run with opencv 3.x.

### Compiling

The ellipse.pyx module must be compiled using Cython, to do this run:

```
./compile.sh
```

You **must** run this whenever you make changes in ellipse.pyx.

### On Pan

To run under Slurm:
```
sbatch slurm/run_pan.sl
```
To load dependencies (e.g. for testing on a build node):
```
module load OpenCV/3.4.0-gimkl-2017a-Python-2.7.13

### On kupe CS500 node (hafs01)

```
module load Anaconda2/5.0.1-GCC-4.8.5
source activate my_conda
python tracking.py -d1 2010-02-19 -d2 2010-02-20 -lons 1700:2200 -lats 200:500 -suffix png_3.0-8
```
Note: The starting date (-d1), the ending date (-d2) and the configuration have default values "2010-02-19", "2010-02-23" and "config23" respectively so it is no longer necessary to pass these arguments. As of March 2 2018, it is possible to specify the longitude/latitude index box with options -lons and -lats respectively. Type Tracking_clean.py --help for a full set of options. 

## How to profile

To generate a graphical representation of execution time spent in each function
```
python -m cProfile -o output.pstats tracking.py [options]
gprof2dot --colour-nodes-by-selftime -f pstats output.pstats | dot -Tpng -o output.png
```
 

