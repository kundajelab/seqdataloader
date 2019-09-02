from __future__ import division, print_function, absolute_import
from .core import Coordinates
import numpy as np


def get_revcomp(coordinate):
    return Coordinates(chrom=coordinate.chrom,
                       start=coordinate.start, end=coordinate.end,
                       isplusstrand=(coordinate.isplusstrand==False))


class AbstractCoordBatchTransformer(object):
  
    def __call__(self, coords):
        """
        Args:
            coords (:obj:`list` of :obj:`Coordinates` objects):

        Returns:
            another :obj:`list` of :obj:`Coordinates`
        """
        raise NotImplementedError()
    
    def chain(self, coord_batch_transformer):
        return lambda coords: coord_batch_transformer(self(coords))
      
      
class ReverseComplementAugmenter(AbstractCoordBatchTransformer):
    """
        Returns a list of Coordinates twice the length of the
            original list by appending the reverse complements
            of the original coordinates at the end
    """
    def __call__(self, coords):
        return coords + [get_revcomp(x) for x in coords]
      
      
class UniformJitter(AbstractCoordBatchTransformer):
    
    def __init__(self, maxshift, seed=1234, chromsizes_file=None):
        """
          Returns a list of Coordinates jittered relative to the original
            coordinates by a shift of up to +/- maxshift. Size of the
            shift is sampled from a uniform distribution.
            
          Args:
            maxshift (:obj:`int`): maximum possible shift to sample
            chromsizes (:obj:`string`): path to a chromsizes file. If
                specified, shifts will be adjusted so as to avoid going
                over the end of the chromosome. Default is None.
        """
        self.rng = np.random.RandomState(seed)
        self.maxshift = maxshift
        self.chromsizes = (
            self._read_chromsizes(chromsizes_file=chromsizes_file)
            if chromsizes_file is not None else None)
    
    def _read_chromsizes(self, chromsizes_file):
        chrom_to_size = {}
        for row in open(chromsizes_file):
            chrom,chromlen = row.rstrip().split("\t")
            chromlen = int(chromlen)
            chrom_to_size[chrom] = chromlen
        return chrom_to_size
   
    def __call__(self, coords):
        a_list = []
        for coord in coords:
            chrom = coord.chrom
            start = coord.start
            end = coord.end
            isplusstrand = coord.isplusstrand
            shift_size = int(self.rng.uniform(low=0, high=(2*self.maxshift + 1))
                             - self.maxshift)
            shift_size = max(-start, shift_size)
            if self.chromsizes is not None:
                shift_size = min(self.chromsizes[chrom]-end, shift_size)
            start = start + shift_size
            end = end + shift_size
            a_list.append(Coordinates(chrom=chrom, start=start,
                                      end=end, isplusstrand=isplusstrand))
        return a_list
