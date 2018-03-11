from matplotlib import pylab

num_cells = [1000*500, 2000*1000, 4000*2000, 8000*4000, 16000*8000, 32000*16000,]
old_times = [ 0.0213229656219, 0.065455198288, 0.287713050842, 1.21719193459, 4.99809193611, 20.9846560955]
new_times = [0.0230829715729, 0.0230529308319, 0.0226380825043, 0.0228049755096, 0.0229768753052, 0.0227789878845]

pylab.loglog(num_cells, old_times, 'g-', num_cells, new_times, 'r-')
pylab.legend(['old', 'new'])
pylab.loglog(num_cells, old_times, 'kx', num_cells, new_times, 'k+')
pylab.xlabel('number of domain cells')
pylab.ylabel('time to compute ellipse (s)')
pylab.title('oblique rectangle 11000 cells')
pylab.show()
