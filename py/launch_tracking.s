## Test Borneo
#python tracking.py -d1 2010-01-12 -d2 2010-01-13 -lons 1450:1700 -lats 300:550 -suffix borneo -harvest 24

## Test Philippines
#python tracking.py -d1 2013-11-05 -d2 2013-11-08 -lons 1400:1900 -lats 350:750 -suffix phil -harvest 24

## Test Papua New Guinea
#python tracking.py -d1 2010-02-19 -d2 2010-02-20 -lons 1700:2200 -lats 200:500 -suffix png -harvest 24

## Test Sumatra
#python tracking.py -d1 2008-11-04 -d2 2008-11-06 -lons 1250:1550 -lats 250:550 -suffix sumatra \
#                     -harvest 48 #-restart_dir restart -restart_interval 3
#python tracking.py -d1 2008-11-04 -d2 2008-11-06 -lons 1700:2200 -lats 200:500 -suffix png -harvest 24

## Test Madagascar
#python tracking.py -d1 1999-01-01 -d2 1999-01-04 -lons 250:750 -lats 0:300 -suffix madag -harvest 24


## Four final areas
# Indian_ocean + Maritime Continent
python tracking.py -d1 1999-01-01 -d2 1999-01-03 -lons 335:2717 -lats 0:827 -suffix io-cm -harvest 24
# Pacific ocean
#python tracking.py -d1 1999-01-01 -d2 1999-01-03 -lons 2647:3216 -lats 0:827 -suffix pacific -harvest 24
# Americas
#python tracking.py -d1 1999-01-01 -d2 1999-01-03 -lons 3216:4577 -lats 0:827 -suffix america -harvest 24
# Africa
#python tracking.py -d1 1999-01-01 -d2 1999-01-03 -lons 4497:335 -lats 0:827 -suffix africa -harvest 24

#python tracking.py -d1 2011-12-16 -d2 2011-12-16 -lons 0:4948 -lats 0:827 -suffix debug \
#                    -harvest 24 #-restart_dir restart -restart_interval 1
