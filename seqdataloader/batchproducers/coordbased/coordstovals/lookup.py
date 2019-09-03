from __future__ import division, print_function, absolute_import
import numpy as np
from .core import AbstractSingleNdarrayCoordsToVals
from ..core import Coordinates


class SimpleLookup(AbstractSingleNdarrayCoordsToVals):

    def __init__(self, lookup_file, default_returnval=0.0):
        self.lookup_file = lookup_file
        self.default_returnval = default_returnval
        self.lookup = {}
        self.num_labels = None
        for line in (gzip.open(self.lookup_file) if ".gz"
                     in self.lookup_file else open(self.lookup_file)):
            (chrom, start_str, end_str, *labels) =\
              line.decode("utf-8").rstrip().split("\t")[0:3]
            coord = Coordinates(chrom=chrom,
                                start=int(start_str),
                                end=int(end_str))
            labels = [float(x) for x in labels] 
            self.lookup[coord] = labels
            if (self.num_labels is None):
                self.num_labels = len(labels)
            else:
                assert len(labels)==self.num_labels,(
                  "Unequal label lengths; "+str(len(labels), self.num_labels))
    
    def _get_ndarray(self, coors):
        to_return = []
        for coor in coors:
            if coor not in self.lookup:
                to_return.append(np.ones(self.num_labels)
                                 *self.default_returnval)
            else:
                to_return.append(self.lookup[coor])
        return np.array(to_return)
