## Test Borneo
#python testCmorph.py -d1 2010-01-12 -d2 2010-01-13 -lons 1450:1700 -lats 300:550 -suffix borneo

## Test Philippines
#python testCmorph.py -d1 2013-11-05 -d2 2013-11-08 -lons 1400:1900 -lats 350:750 -suffix phil

## Test Papua New Guinea 
#python testCmorph.py -d1 2010-02-19 -d2 2010-02-20 -lons 1700:2200 -lats 200:500 -suffix png -min_axis 6 

## Test Sumatra
#python testCmorph.py -d1 2008-11-04 -d2 2008-11-06 -lons 1250:1550 -lats 250:550 -suffix sumatra

## Test Madagascar
#python testCmorph.py -d1 1999-01-01 -d2 1999-01-01 -lons 250:750 -lats 0:300 -suffix madag

## Larger area
python -m cProfile testCmorph.py -d1 1999-01-01 -d2 1999-01-01 -lons 328:2726 -lats 0:818 -suffix all

## Plot ellipses
#python plotEllipses.py -n 162 -f cmorph.pckl_png -t 50:100
