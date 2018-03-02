import numpy
from box import Box

class Blob:

	def __init__(self, ndims, cells={}):
		"""
		Constructor 
        @param ndims number of dimensions
        @param cells set of (i,j...) tuples
		"""
        self.ndims = ndims
        self.cells = cells
        self.lo = []
        self.hi = []
        if len(cells) > 0:
            list_cells = list(cells)
            self.lo = reduce(lambda x,y: numpy.minimum(x, y), list_cells)
            self.hi = reduce(lambda x,y: numpy.maximum(x, y), list_cells)

    def merge(self, otherBlob):
        """
        Merge a blob with another blob
        @param otherBlob other blob
        @return True if this and the other blob can be merged, False otherwise
        """
        # Shift this blob by one cell in every direction. If there is overlap
        # then this and the other blob can be merged. 
        displ = numpy.zeros((ndims,), numpy.int32)
        for dim in range(self.ndims):
            for pm in (-1, 1):
                offsetBlob = self.shift(dim, pm)
                if offsetBlob.overlaps(otherBlob):
                    self.absorb(otherBlob)
                    return True
        return False

    def shift(self, dim, pm):
        """
        Shift blob
        @param dim dimension (axis)
        @param pm either -1 or 1
        @return new blob that is shifted along dimension "dim" and direction "pm"
        """
        res = Blob(self.ndims)
        displ = numpy.zeros((self.ndims,), numpy.int32)
        # shift the cells up or down
        res.cells = {tuple(numpy.array(c) + displ) for c in self.cells}
        res.lo = self.lo + displ
        res.hi = self.hi + displ
        return res

    def overlaps(self, otherBlob):
        """
        Check if two blobs overlap (share at least one cell)
        @param otherBlob other blob
        @return True if there is an overlap, False otherwise
        """
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
        """
        Absorb blob
        @param otherBlob other blob
        """
        cells = self.cells.union(otherBlob.cells)
        lo = numpy.minimum(self.lo, otherBlob.lo)
        hi = numpy.maximum(self.hi, otherBlob.hi)
        self.cells = cells
        self.lo = lo
        self.hi = hi

