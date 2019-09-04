from __future__ import absolute_import, division, print_function
import gzip
from .core import Coordinates
import numpy as np


class KerasSequenceApiCoordsBatchProducer(object): 

    """
    Args:
        batch_size (int): note that if you apply some kind of augmentation,
            then this value will end up being half of the actual batch size.
        shuffle_before_epoch (boolean, optional): default False
        seed (int): default 1234; needed if shuffle=True
    """
    def __init__(self, batch_size, shuffle_before_epoch, seed):
        self.coords_list = self._get_coordslist()
        self.batch_size = batch_size
        self.shuffle_before_epoch = shuffle_before_epoch
        self.seed = seed
        if (self.shuffle_before_epoch):
            self.rng = np.random.RandomState(self.seed)
            self._shuffle_coordslist()

    def _get_coordslist(self):
        raise NotImplementedError()
    
    def _shuffle_coordslist(self):
        self.rng.shuffle(self.coords_list)
  
    def __getitem__(self, index):
        """
        Args:
            index (:obj:`int`): index of the batch
        
        Returns:
            :obj:`list`: the coordinates for a complete batch
        """
        return self.coords_list[index*self.batch_size:
                                (index+1)*self.batch_size]
    
    def __len__(self):
        """
        Returns:
            The total number of batches to return
        """
        return int(np.ceil(len(self.coords_list)/float(self.batch_size)))
   
    def on_epoch_end(self):
        """
        Things to be executed after the epoch - like shuffling the coords
        """
        if (self.shuffle_before_epoch):
            self._shuffle_coordslist()


class BedFileObj(object):
    def __init__(self, bed_file, hastitle=False):
        print("Heads up: coordinates in bed file"
              +" are assumed to be on the positive strand;"
              +" if strand in the bed file is improtant to you, please"
              +" add that feature to SimpleCoordsBatchProducer")
        self.bed_file = bed_file 
        self.hastitle = hastitle
        self.coords_list = self._read_bed_file()

    def _read_bed_file(self):
        coords_list = []
        for linenum,line in enumerate((gzip.open(self.bed_file) if ".gz"
                                       in self.bed_file
                                       else open(self.bed_file))):
            if (linenum > 0 or self.hastitle==False):
                (chrom, start_str, end_str) =\
                  line.decode("utf-8").rstrip().split("\t")[0:3]
                coords_list.append(Coordinates(chrom=chrom,
                                              start=int(start_str),
                                              end=int(end_str)))
        return coords_list

    def __len__(self):
        return len(self.coords_list)

    def get_strided_subsample(self, offset, stride):
        return self.coords_list[offset::stride]

    def assert_sorted(self):
        prev_entry = self.coords_list[0]
        for entry in self.coords_list[1:]:
            if entry.chrom==prev_entry.chrom:
                assert entry.start >= prev_entry.start, ("Bed file "+
                        self.bed_file+" is not sorted; "+str(entry)
                        +" follows "+str(prev_entry))
            prev_entry = entry
            

class DownsampleNegativesCoordsBatchProducer(
        KerasSequenceApiCoordsBatchProducer):

    def __init__(self, pos_bed_file, neg_bed_file,
                       target_proportion_positives, **kwargs):

        print("Reading in positive bed file")
        self.pos_bedfileobj = BedFileObj(bed_file=pos_bed_file)
        print("Got",len(self.pos_bedfileobj.coords_list),
              " coords in positive bed file")
        print("Reading in negative bed file")
        self.neg_bedfileobj = BedFileObj(bed_file=neg_bed_file)
        print("Got",len(self.neg_bedfileobj.coords_list),
              " coords in negative bed file")
        self.neg_bedfileobj.assert_sorted()

        self.target_proportion_positives = target_proportion_positives
        self.subsample_factor = int(np.ceil(
            (len(self.neg_bedfileobj.coords_list)
             *(self.target_proportion_positives/
               (1-self.target_proportion_positives)) )/
            len(self.pos_bedfileobj.coords_list)))
        print("The target proportion of positives of",
              self.target_proportion_positives,"requires the negative set"
              +" to be subsampled by a factor of",self.subsample_factor,
              "which will result in a #neg of",
              int(len(self.neg_bedfileobj.coords_list)/self.subsample_factor))
        self.last_used_offset = -1
        super(DownsampleNegativesCoordsBatchProducer, self).__init__(**kwargs)

    def _shuffle_coordslist(self):
        self.rng.shuffle(self.subsampled_neg_coords)
        self.rng.shuffle(self.pos_coords)
        fracpos = len(self.pos_coords)/(
                    len(self.pos_coords) + len(self.subsampled_neg_coords))
        #interleave evenly
        pos_included = 0
        neg_included = 0
        new_coordslist = []
        for i in range(len(self.pos_coords)+len(self.subsampled_neg_coords)):
            if (pos_included < (pos_included+neg_included)*(fracpos)):
                new_coordslist.append(self.pos_coords[pos_included])
                pos_included += 1
            else:
                new_coordslist.append(self.subsampled_neg_coords[neg_included])
                neg_included += 1
        assert pos_included==len(self.pos_coords)
        assert neg_included==len(self.subsampled_neg_coords)
        self.coords_list = new_coordslist

    def _get_coordslist(self):
        self.last_used_offset += 1
        self.last_used_offset = self.last_used_offset%self.subsample_factor
        print("Using an offset of ",self.last_used_offset," before striding")
        self.last_used_offset = self.last_used_offset%self.subsample_factor
        subsampled_neg_coords = self.neg_bedfileobj.get_strided_subsample(
                                offset=self.last_used_offset,
                                stride=self.subsample_factor) 
        pos_coords = self.pos_bedfileobj.coords_list
        self.subsampled_neg_coords = subsampled_neg_coords
        self.pos_coords = pos_coords
        return pos_coords+subsampled_neg_coords
   
    def on_epoch_end(self):
        #get negative set with potentially different stride
        self.coords_list = self._get_coordslist()
        #perform shuffling as needed
        super(DownsampleNegativesCoordsBatchProducer, self).on_epoch_end()        


class SimpleCoordsBatchProducer(KerasSequenceApiCoordsBatchProducer):

    """
    Args:
        bed_file (string): file with the bed coordinates.
            Assumes coordinates are on the positive strand.
        coord_batch_transformer (AbstracCoordBatchTransformer): does things
            like revcomp and random jitter
    """
    def __init__(self, bed_file,
                       hastitle=False,
                       coord_batch_transformer=None,
                       **kwargs):  
        self.bed_file = BedFileObj(bed_file=bed_file, hastitle=hastitle)
        if (coord_batch_transformer is not None):
            raise DeprecationWarning(
             "Moving forward, coords_batch_transformer should be"
             +" specified as an argument to KerasBatchGenerator"
             +", not as an arugment to the CoordsBatchProducer."
             +" This is to allow different CoordsBatchProducer"
             +" implementations to be used with the same"
             +" coords_batch_transformer code.")
        self.coord_batch_transformer = coord_batch_transformer
        super(SimpleCoordsBatchProducer, self).__init__(**kwargs)

    def _get_coordslist(self):
        return [x for x in self.bed_file.coords_list]
        
    def __getitem__(self, index):
        orig_batch = self.coords_list[index*self.batch_size:
                                      (index+1)*self.batch_size]
        if (self.coord_batch_transformer is not None):
            return self.coord_batch_transformer(orig_batch)
        else:
            return orig_batch
