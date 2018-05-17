## Test Borneo
#python tracking_pickle.py -d1 2010-01-12 -d2 2010-01-13 -lons 1450:1700 -lats 300:550 -suffix borneo -harvest 12

## Test Philippines
#python tracking_pickle.py -d1 2013-11-05 -d2 2013-11-08 -lons 1400:1900 -lats 350:750 -suffix phil -harvest 12

## Test Papua New Guinea
python tracking_pickle.py -d1 2010-02-19 -d2 2010-02-20 -lons 1700:2200 -lats 200:500 -suffix png -harvest 12

## Test Sumatra
#python tracking_pickle.py -d1 2008-11-04 -d2 2008-11-06 -lons 1250:1550 -lats 250:550 -suffix sumatra -harvest 12

## Test Madagascar
#python tracking_pickle.py -d1 1999-01-01 -d2 1999-01-04 -lons 250:750 -lats 0:300 -suffix madag_3.0-8

## Test Maritime Continent
#python gprof2dot.py --colour-nodes-by-selftime -f pstats output.pstats | dot -Tpng -o output.png
#python -m cProfile -o output.pstats tracking.py -d1 1999-01-01 -d2 1999-01-03 -lons 1700:2200 -lats 50:650 -suffix cm -harvest 24

## Larger area
#python -m cProfile tracking_pickle.py -d1 1999-01-01 -d2 1999-01-01 -lons 328:2726 -lats 0:818 -suffix all

## Profiling
#sh5util -j "${SLURM_JOB_ID}" -o profile.h5
#python profile_plot.py profile.h5

