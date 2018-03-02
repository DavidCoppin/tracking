import numpy
from box import Box

class Blob:

	def __init__(self, ndims):
		"""
		Constructor 
		"""
        self.ndims = ndims
        self.cells = {}
        self.lo = []
        self.hi = []

    def merge(self, otherBlob):
        displ = numpy.zeros((ndims,), numpy.int32)
        for dim in range(self.ndims):
            for pm in (-1, 1):
                offsetBlob = self.shift(dim, pm)
                if offsetBlob.overlaps(otherBlob):
                    self.absorb(otherBlob)
                    return True
        return False

    def shift(self, dim, pm):
        res = Blob(self.ndims)
        displ = numpy.zeros((self.ndims,), numpy.int32)
        res.cells = {tuple(numpy.array(c) + displ) for c in self.cells}
        res.lo = self.lo + displ
        res.hi = self.hi + displ
        return res

    def overlaps(self, otherBlob):
        # quick test
        if numpy.all(self.hi < otherBlob.lo):
            return False
        if numpy.all(self.lo > otherBlob.hi):
            return False
        # there is a chance of an overlap
        s = self.cells.intersect(otherBlob.cells)
        if len(s) > 0:
            return True
        return False

    def absorb(self, otherBlob):
        cells = self.cells.union(otherBlob.cells)
        lo = numpy.minimum(self.lo, otherBlob.lo)
        hi = numpy.maximum(self.hi, otherBlob.hi)
        self.cells = cells
        self.lo = lo
        self.hi = hi









		indices = numpy.where(data[slab] > 0)

		self.ndims = len(indices)
		numCells = len(indices[0])

		# set of cell indices where data is nonzero
		self.cells = {tuple([indices[i][j] for i in range(ndims)]) for j in range(numCells)}

		# window of interest
		self.lo = numpy.array([min([cell[i] for cell in self.cells])])
		self.hi = numpy.array([max([cell[i] for cell in self.cells])])


	def build(self):
		"""
		Build the object by fusing neighbouring cells that share a common face
		"""
		while self.merge():
			pass

    def merge(self):



	def touches(self, otherBlob):
		pass

	def absorbs(self, otherBlob):
		pass