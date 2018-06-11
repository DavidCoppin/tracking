set -x

# run without restart
#time python tracking.py -d1 2008-11-04 -d2 2008-11-06 -lons 1250:1550 -lats 250:550 -suffix sumatra -harvest 24

# run with a restart
# first delete any current restart files
# writing out restarts every day (just while testing)
rm -rf restart.test
rm -f pickles/sumatra*
echo "Running for first 2 days..."
time python tracking.py -d1 2008-11-04 -d2 2008-11-05 -lons 1250:1550 -lats 250:550 -suffix sumatra -harvest 24 -restart-dir restart.test -restart-interval 1
echo "Running for last day..."
time python tracking.py -d1 2008-11-04 -d2 2008-11-06 -lons 1250:1550 -lats 250:550 -suffix sumatra -harvest 24 -restart-dir restart.test -restart-interval 1
