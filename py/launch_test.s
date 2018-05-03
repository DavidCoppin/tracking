## Test Borneo
#python tracking.py -d1 2010-01-12 -d2 2010-01-13 -lons 1450:1700 -lats 300:550 -suffix borneo_3.0-8

## Test Philippines
#python tracking.py -d1 2013-11-05 -d2 2013-11-08 -lons 1400:1900 -lats 350:750 -suffix phil_3.0-8

## Test Papua New Guinea
python tracking.py -d1 2010-02-19 -d2 2010-02-20 -lons 1700:2200 -lats 200:500 -suffix png_3.0-8

## Test Sumatra
#python tracking.py -d1 2008-11-04 -d2 2008-11-06 -lons 1250:1550 -lats 250:550 -suffix sumatra_3.0-8

## Test Madagascar
#python tracking.py -d1 1999-01-01 -d2 1999-01-04 -lons 250:750 -lats 0:300 -suffix madag_3.0-8

## Test Maritime Continent
#python tracking.py -d1 1999-01-01 -d2 1999-01-03 -lons 1700:2200 -lats 550:650 -suffix cm

## Larger area
#python -m cProfile tracking.py -d1 1999-01-01 -d2 1999-01-01 -lons 328:2726 -lats 0:818 -suffix all

## Plot ellipses
#python plotEllipses.py -n 162 -f cmorph.pckl_png -t 50:100
