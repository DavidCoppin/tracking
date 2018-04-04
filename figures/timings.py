from matplotlib import pylab
import numpy

domain_size = numpy.array([(550-300)*(1700 -1450), (550-300)*(1700 - 0), (550-300)*4948, 827*4948])

time_kupe_clusters = numpy.array([50.7, 2*60+1.5, 13*60+45.6, 78*60+59.6])
time_kupe_tracking = numpy.array([2*60+50.3, 6*60+53.5, 64*60+50.3, 451*60+3.3])
#time_pan_tracking = [ ,  ,  , ]

pylab.plot(domain_size/1000000., time_kupe_tracking/time_kupe_clusters, 'ks')
pylab.plot(domain_size/1000000., time_kupe_tracking/time_kupe_clusters, 'b-')
pylab.title('speedup vs domain size')
pylab.xlabel('number of cells (million)')
pylab.show()



