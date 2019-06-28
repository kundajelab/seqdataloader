from __future__ import division, print_function, absolute_import
import keras


class Coordinates(object):
    
    def __init__(self, chrom, start, end, isplusstrand=True):
        """
        Args:
            chrom (string)
            start (int)
            end (int)
            isplusstrand (boolean, optional): default True
            
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.isplusstrand = isplusstrand
    
    def __str__(self):
       return (self.chrom+":"+str(self.isplusstrand)+":"
               +str(self.start)+"-"+str(self.end))
    
    def __repr__(self):
       return self.__str__()
    
    def get_revcomp(self):
        return Coordinates(chrom=self.chrom, start=self.start, end=self.end,
                           isplusstrand=(self.isplusstrand==False))


class KerasBatchGenerator(keras.utils.Sequence):
  
    """
    Args:
        coordsbatch_producer (KerasSequenceApiCoordsBatchProducer)
        inputs_coordstovals (CoordsToVals)
        targets_coordstovals (CoordsToVals)
        sampleweights_coordstovals (CoordsToVals)
    """
    def __init__(self, coordsbatch_producer,
                       inputs_coordstovals,
                       targets_coordstovals,
                       sampleweights_coordstovals=None):
        self.coordsbatch_producer = coordsbatch_producer
        self.inputs_coordstovals = inputs_coordstovals
        self.targets_coordstovals = targets_coordstovals
        self.sampleweights_coordstovals = sampleweights_coordstovals
    
    def __getitem__(self, index):
        coords_batch = self.coordsbatch_producer[index]
        inputs = self.inputs_coordstovals(coords_batch)
        targets = self.targets_coordstovals(coords_batch)
        if (self.sampleweights_coordstovals is not None):
            sample_weights = self.sampleweights_coordstovals(coords_batch)
            return (inputs, targets, sample_weights)
        else:
            return (inputs, targets)
   
    def __len__(self):
        return len(self.coordsbatch_producer)
    
    def on_epoch_end(self):
        self.coordsbatch_producer.on_epoch_end()
