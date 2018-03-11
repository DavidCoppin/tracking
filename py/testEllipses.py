import numpy
import cv2
from ellipse import Ellipse
from opencv import getContours, getEllipse
import matplotlib.pyplot as mpl
import matplotlib
import math
import time

"""
Compare ellipses
"""

class EllipseTester:
    def __init__(self, label, domain_size, cells):

        self.domain_size = domain_size
        self.label = label
        self.cells = cells

        self.data = {
            'old': {'centre':[0., 0.], 'axes':[0., 0.], 'time': None, 'x': [], 'y': []},
            'new': {'centre':[None, None], 'axes':[None, None], 'time': None, 'x': [], 'y': []},
        }

        # new
        tic = time.time()
        ell = Ellipse(cells)
        toc = time.time()
        self.data['new']['time']= toc - tic

        self.data['new']['centre'] = ell.centre
        self.data['new']['axes'] = [ell.a, ell.b]
        self.data['new']['x'], self.data['new']['y'] = ell.getPolyline()

        # old, not always working
        cellData = numpy.zeros(domain_size, numpy.int32)
        for c in cells:
            cellData[c] = 1
        try:
            tic = time.time()
            test_result, test_contours = getContours(cellData)
            ellipse_parameters = getEllipse(test_result,test_contours)
            toc = time.time()
            self.data['old']['time'] = toc - tic

            # note: reversed order
            yc, xc = ellipse_parameters[0]
            b, a = ellipse_parameters[1]
            b /= 2.0
            a /= 2.0
            angle = ellipse_parameters[2] # angle in degrees

            self.data['old']['centre'] = (xc, yc)
            self.data['old']['axes'] = (a, b)

            # constructor contour
            cos_angle = math.cos(angle*math.pi/180.)
            sin_angle = math.sin(angle*math.pi/180.)
            nt = 32
            ts = numpy.linspace(0., 2*math.pi, nt + 1)
            dxPrime = a*numpy.cos(ts)
            dyPrime = b*numpy.sin(ts)
            dx =  dxPrime * cos_angle + dyPrime * sin_angle
            dy = -dxPrime * sin_angle + dyPrime * cos_angle
            self.data['old']['x'] = xc + dx
            self.data['old']['y'] = yc + dy
        except:
            print('Failed to use old method!')

    def __repr__(self):
        area_new = math.pi*self.data['new']['axes'][0]*self.data['new']['axes'][1]
        area_old = math.pi*self.data['old']['axes'][0]*self.data['old']['axes'][1]
        res = """
EllipseTester {}
center (old new): {}    ; {}
axes   (old new): {}    ; {}
area   (old new): {}    ; {} exact: {}
times  (old new): {}s   ; {}s for domain size {}
        """.format(self.label,
                   self.data['old']['centre'], self.data['new']['centre'],
                   self.data['old']['axes'], self.data['new']['axes'],
                   area_old, area_new, len(self.cells),
                   self.data['old']['time'], self.data['new']['time'],
                   self.domain_size)
        return res

    def plot(self):

        fig, ax = mpl.subplots()
        
        # plot the cells
        patches = []
        for c in self.cells:
            i, j = c
            pts = numpy.array([[i-0.5, j-0.5], [i+0.5, j-0.5], [i+0.5, j+0.5], [i-0.5, j+0.5]])
            patch = matplotlib.patches.Polygon(pts, closed=True)
            patches.append(patch)
        p = matplotlib.collections.PatchCollection(patches, alpha=0.4)
        ax.add_collection(p)

        # plot the old ellipse
        mpl.plot(self.data['old']['x'], self.data['old']['y'], 'g-')

        # plot the new ellipse
        mpl.plot(self.data['new']['x'], self.data['new']['y'], 'r-')

        # add labels
        mpl.legend(['old', 'new'])

        # add the centres
        mpl.plot([self.data['old']['centre'][0]], [self.data['old']['centre'][1]], 'gx')
        mpl.plot([self.data['new']['centre'][0]], [self.data['new']['centre'][1]], 'r+')

        mpl.xlabel('i')
        mpl.ylabel('j')
        mpl.title(self.label)

        # show the result
        mpl.show()



#####################################################################
def testObliqueRectangle():
    domain_size = (30, 40)
    cells = set()
    for j in range(10):
        cells = cells.union({(i+j, i) for i in range(3, 14)})
    et = EllipseTester('oblique rectangle', domain_size, cells)
    print et
    et.plot()

def testConcave():
    domain_size = (30, 40)
    x, y = 10, 15
    cells = {(-2+x,3+y), (-3+x,3+y), (-3+x,2+y), (-3+x,1+y),
         (-3+x,0+y), (-2+x,0+y), (-1+x,0+y), (0+x,0+y), (1+x,0+y), (2+x,0+y), (3+x,0+y),
         (3+x,1+y), (3+x,2+y), (2+x,2+y), (1+x,2+y), (0+x,2+y),}
    et = EllipseTester('concave', domain_size, cells)
    print et
    et.plot()

def testConcave2():
    domain_size = (30, 40)
    x, y = 10, 0
    cells = {(-2+x,3+y), (-3+x,3+y), (-3+x,2+y), (-3+x,1+y),
         (-3+x,0+y), (-2+x,0+y), (-1+x,0+y), (0+x,0+y), (1+x,0+y), (2+x,0+y), (3+x,0+y),
         (3+x,1+y), (3+x,2+y), (2+x,2+y), (1+x,2+y), (0+x,2+y),}
    et = EllipseTester('concave', domain_size, cells)
    print et
    et.plot()

def testIslands():
    domain_size = (30, 40)
    cells = {(0,0), (1,0), (1,1), (1,4)}
    et = EllipseTester('islands', domain_size, cells)
    print et
    et.plot()

def testLargeDomain(domain_size):
    cells = set()
    for j in range(100):
        cells = cells.union({(i+j, i) for i in range(30, 140)})
    et = EllipseTester('large domain {}'.format(domain_size), domain_size, cells)
    print et
    #et.plot()

def testWrongEllipse():
    domain_size = (400, 1000)
    cells = {
    (241, 113), 
    (241, 114), 
    (242, 112), 
    (242, 113), 
    (242, 114), 
    (242, 115), 
    (242, 116), 
    (242, 118), 
    (242, 119), 
    (242, 120), 
    (242, 121), 
    (243, 111), 
    (243, 112), 
    (243, 113), 
    (243, 114), 
    (243, 115), 
    (243, 116), 
    (243, 117), 
    (243, 118), 
    (243, 119), 
    (243, 120), 
    (243, 121), 
    (244, 110), 
    (244, 111), 
    (244, 112), 
    (244, 113), 
    (244, 114), 
    (244, 115), 
    (244, 116), 
    (244, 117), 
    (244, 118), 
    (244, 119), 
    (244, 120), 
    (245, 109), 
    (245, 110), 
    (245, 111), 
    (245, 112), 
    (245, 113), 
    (245, 114), 
    (245, 115), 
    (245, 116), 
    (245, 117), 
    (245, 118), 
    (245, 119), 
    (245, 120), 
    (246, 109), 
    (246, 110), 
    (246, 111), 
    (246, 112), 
    (246, 113), 
    (246, 114), 
    (246, 115), 
    (246, 116), 
    (246, 117), 
    (246, 118), 
    (246, 119), 
    (247, 109), 
    (247, 110), 
    (247, 111), 
    (247, 112), 
    (247, 113), 
    (247, 114), 
    (247, 115), 
    (247, 116), 
    (247, 117), 
    (247, 118), 
    (247, 119), 
    (248, 108), 
    (248, 109), 
    (248, 110), 
    (248, 111), 
    (248, 112), 
    (248, 113), 
    (248, 114), 
    (248, 115), 
    (248, 116), 
    (248, 117), 
    (248, 118), 
    (248, 119), 
    (249, 109), 
    (249, 110), 
    (249, 111), 
    (249, 112), 
    (249, 113), 
    (249, 114), 
    (249, 115), 
    (249, 116), 
    (249, 117), 
    (249, 118), 
    (249, 119), 
    (250, 109), 
    (250, 110), 
    (250, 111), 
    (250, 112), 
    (250, 113), 
    (250, 114), 
    (250, 115), 
    (250, 116), 
    (250, 117), 
    (250, 118), 
    (250, 119), 
    (250, 120), 
    (251, 109), 
    (251, 110), 
    (251, 111), 
    (251, 112), 
    (251, 113), 
    (251, 114), 
    (251, 115), 
    (251, 116), 
    (251, 117), 
    (251, 118), 
    (251, 119), 
    (251, 120), 
    (252, 109), 
    (252, 110), 
    (252, 111), 
    (252, 112), 
    (252, 113), 
    (252, 114), 
    (252, 115), 
    (252, 116), 
    (252, 117), 
    (252, 118), 
    (252, 119), 
    (252, 120), 
    (252, 121), 
    (253, 109), 
    (253, 112), 
    (253, 113), 
    (253, 114), 
    (253, 115), 
    (253, 116), 
    (253, 117), 
    (253, 118), 
    (253, 119), 
    (253, 120), 
    (253, 121), 
    (254, 112), 
    (254, 113), 
    (254, 114), 
    (254, 115), 
    (254, 116), 
    (254, 117), 
    (254, 118), 
    (254, 119), 
    (254, 120), 
    (254, 121), 
    (255, 110), 
    (255, 112), 
    (255, 113), 
    (255, 114), 
    (255, 115), 
    (255, 116), 
    (255, 117), 
    (255, 118), 
    (255, 119), 
    (255, 120), 
    (255, 121), 
    (255, 122), 
    (256, 108), 
    (256, 110), 
    (256, 111), 
    (256, 112), 
    (256, 113), 
    (256, 114), 
    (256, 115), 
    (256, 116), 
    (256, 117), 
    (256, 118), 
    (256, 119), 
    (256, 120), 
    (256, 121), 
    (256, 122), 
    (256, 123), 
    (257, 106), 
    (257, 107), 
    (257, 108), 
    (257, 109), 
    (257, 110), 
    (257, 111), 
    (257, 112), 
    (257, 113), 
    (257, 114), 
    (257, 115), 
    (257, 116), 
    (257, 117), 
    (257, 118), 
    (257, 119), 
    (257, 120), 
    (257, 121), 
    (257, 122), 
    (257, 123), 
    (257, 124), 
    (258, 103), 
    (258, 104), 
    (258, 105), 
    (258, 106), 
    (258, 107), 
    (258, 108), 
    (258, 109), 
    (258, 110), 
    (258, 111), 
    (258, 112), 
    (258, 113), 
    (258, 114), 
    (258, 115), 
    (258, 116), 
    (258, 117), 
    (258, 118), 
    (258, 119), 
    (258, 120), 
    (258, 121), 
    (258, 122), 
    (258, 123), 
    (258, 124), 
    (259, 103), 
    (259, 104), 
    (259, 105), 
    (259, 106), 
    (259, 107), 
    (259, 108), 
    (259, 109), 
    (259, 110), 
    (259, 111), 
    (259, 112), 
    (259, 113), 
    (259, 114), 
    (259, 115), 
    (259, 116), 
    (259, 117), 
    (259, 118), 
    (259, 119), 
    (259, 120), 
    (259, 121), 
    (259, 122), 
    (259, 123), 
    (259, 124), 
    (260, 106), 
    (260, 107), 
    (260, 108), 
    (260, 109), 
    (260, 110), 
    (260, 111), 
    (260, 112), 
    (260, 113), 
    (260, 114), 
    (260, 115), 
    (260, 116), 
    (260, 117), 
    (260, 118), 
    (260, 119), 
    (260, 120), 
    (260, 121), 
    (260, 122), 
    (260, 123), 
    (260, 124), 
    (261, 107), 
    (261, 110), 
    (261, 111), 
    (261, 112), 
    (261, 113), 
    (261, 114), 
    (261, 115), 
    (261, 116), 
    (261, 117), 
    (261, 118), 
    (261, 119), 
    (261, 120), 
    (261, 121), 
    (261, 122), 
    (261, 123), 
    (261, 124), 
    (261, 125), 
    (262, 115), 
    (262, 116), 
    (262, 117), 
    (262, 118), 
    (262, 119), 
    (262, 120), 
    (262, 121), 
    (262, 122), 
    (262, 123), 
    (262, 124), 
    (262, 125), 
    (262, 126), 
    (262, 127), 
    (262, 128), 
    (263, 116), 
    (263, 117), 
    (263, 118), 
    (263, 119), 
    (263, 120), 
    (263, 121), 
    (263, 122), 
    (263, 123), 
    (263, 127), 
    (263, 128), 
    (263, 129), 
    (263, 130), 
    (263, 131), 
    (263, 132), 
    (263, 133), 
    (263, 134), 
    (263, 135), 
    (264, 117), 
    (264, 118), 
    (264, 119), 
    (264, 120), 
    (264, 121), 
    (264, 122), 
    (264, 123), 
    (264, 127), 
    (264, 128), 
    (264, 129), 
    (264, 130), 
    (264, 131), 
    (264, 132), 
    (264, 133), 
    (264, 134), 
    (264, 135), 
    (264, 136), 
    (264, 138), 
    (264, 139), 
    (265, 118), 
    (265, 119), 
    (265, 120), 
    (265, 121), 
    (265, 122), 
    (265, 123), 
    (265, 131), 
    (265, 134), 
    (265, 135), 
    (265, 136), 
    (265, 137), 
    (265, 138), 
    (265, 139), 
    (265, 140), 
    (265, 141), 
    (265, 142), 
    (265, 143), 
    (266, 120), 
    (266, 121), 
    (266, 122), 
    (266, 135), 
    (266, 136), 
    (266, 137), 
    (266, 138), 
    (266, 139), 
    (266, 140), 
    (266, 141), 
    (266, 142), 
    (266, 143), 
    (266, 144), 
    (266, 145), 
    (267, 121), 
    (267, 122), 
    (267, 142), 
    (267, 143), 
    (267, 144), 
    (267, 145), 
    (267, 146), 
    (267, 147), 
    (267, 148), 
    (267, 149), 
    (267, 150), 
    (267, 151), 
    (268, 147), 
    (268, 148), 
    (268, 149), 
    (268, 150), 
    (268, 151), 
    (268, 152), 
    (268, 153), 
    (268, 154), 
    (268, 155), 
    (268, 156), 
    (268, 157), 
    (269, 149), 
    (269, 151), 
    (269, 152), 
    (269, 153), 
    (269, 154), 
    (269, 155), 
    (269, 156), 
    (269, 157), 
    (270, 156), 
    }    
    et = EllipseTester('wrong ellipse', domain_size, cells)
    print et
    et.plot()


if __name__ == '__main__':
    testObliqueRectangle()
    testConcave()
    testConcave2()
    testIslands()
    testWrongEllipse()
    testLargeDomain([1000, 500])
    testLargeDomain([2000, 1000])
    testLargeDomain([4000, 2000])
    testLargeDomain([8000, 4000])
    testLargeDomain([16000, 8000])
    testLargeDomain([32000, 16000])
    testLargeDomain([64000, 32000])
