from __future__ import absolute_import, division, print_function
import gzip
from .core import Coordinates
import numpy as np


class KerasSequenceApiCoordsBatchProducer(object): 
  
    def __getitem__(self, index):
        """
        Args:
            index (:obj:`int`): index of the batch
        
        Returns:
            :obj:`list`: the coordinates for a complete batch
        """
        raise NotImplementedError()
    
    def __len__(self):
        """
        Returns:
            The total number of batches to return
        """
        raise NotImplementedError()
   
    def on_epoch_end(self):
        """
        Things to be executed after the epoch - like shuffling the coords
        """
        raise NotImplementedError()


class SimpleCoordsBatchProducer(KerasSequenceApiCoordsBatchProducer):

    """
    Args:
        bed_file (string): file with the bed coordinates.
            Assumes coordinates are on the positive strand.
        batch_size (int): note that if you apply some kind of augmentation,
            then this value will end up being half of the actual batch size.
        coord_batch_transformer (AbstracCoordBatchTransformer): does things
            like revcomp and random jitter
        shuffle_before_epoch (boolean, optional): default False
        seed (int): default 1234; needed if shuffle=True
    """
    def __init__(self, bed_file,
                       batch_size,
                       hastitle=False,
                       coord_batch_transformer=None,
                       shuffle_before_epoch=False,
                       seed=1234):  
        print("Heads up: coordinates in bed file"
              +" are assumed to be on the positive strand;"
              +" if strand in the bed file is improtant to you, please"
              +" add that feature to SimpleCoordsBatchProducer")
        self.bed_file = bed_file
        self.batch_size = batch_size
        self.hastitle = hastitle
        self.coord_batch_transformer = coord_batch_transformer
        self.coords_list = self._read_bed_file(bed_file=self.bed_file)
        self.shuffle_before_epoch = shuffle_before_epoch
        self.seed = seed
        if (self.shuffle_before_epoch):
            self.rng = np.random.RandomState(self.seed)
            self._shuffle_coords_list()

    def _read_bed_file(self, bed_file):
        coords_list = []
        for linenum,line in enumerate((gzip.open(bed_file) if ".gz"
                             in bed_file else open(bed_file))):
            if (linenum > 0 or self.hasttitle==False):
                (chrom, start_str, end_str) =\
                  line.decode("utf-8").rstrip().split("\t")[0:3]
                coords_list.append(Coordinates(chrom=chrom,
                                              start=int(start_str),
                                              end=int(end_str)))
        return coords_list
    
    def _shuffle_coords_list(self):
        self.rng.shuffle(self.coords_list)
        
    def __getitem__(self, index):
        orig_batch = self.coords_list[index*self.batch_size:
                                      (index+1)*self.batch_size]
        if (self.coord_batch_transformer is not None):
            return self.coord_batch_transformer(orig_batch)
        else:
            return orig_batch
        
    def __len__(self):
        return int(np.ceil(len(self.coords_list)/float(self.batch_size)))
   
    def on_epoch_end(self):
        if (self.shuffle_before_epoch):
            self._shuffle_coords_list()
