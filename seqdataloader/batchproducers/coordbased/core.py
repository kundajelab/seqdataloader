from __future__ import division, print_function, absolute_import
from collections import namedtuple
import keras


Coordinates = namedtuple("Coordinates",
                         ["chrom", "start", "end", "isplusstrand"])
Coordinates.__new__.__defaults__ = (True,)


def apply_mask(tomask, mask):
    if isinstance(tomask, dict):
        return dict([(key, val[mask]) for key,val in tomask.items()])
    elif isinstance(tomask, list):
        return [x[mask] for x in mask]
    else:
        return x[mask]


class KerasBatchGenerator(keras.utils.Sequence):
  
    """
    Args:
        coordsbatch_producer (KerasSequenceApiCoordsBatchProducer)
        inputs_coordstovals (CoordsToVals)
        targets_coordstovals (CoordsToVals)
        sampleweights_coordstovals (CoordsToVals)
        coordsbatch_transformer (AbstracCoordBatchTransformer)
        qc_func (callable): function that can be used to filter
            out poor-quality sequences.
        sampleweights_coordstoval: either this argument or
            sampleweights_from_inputstargets could be used to
            specify sample weights. sampleweights_coordstoval
            takes a batch of coords as inputs.
        sampleweights_from_inputstargets: either this argument or
            sampleweights_coordstoval could be used to
            specify sample weights. sampleweights_from_inputstargets
            takes the inputs and targets values to generate the weights.
    """
    def __init__(self, coordsbatch_producer,
                       inputs_coordstovals,
                       targets_coordstovals,
                       coordsbatch_transformer=None,
                       qc_func=None,
                       sampleweights_coordstovals=None,
                       sampleweights_from_inputstargets=None):
        self.coordsbatch_producer = coordsbatch_producer
        self.inputs_coordstovals = inputs_coordstovals
        self.targets_coordstovals = targets_coordstovals
        self.coordsbatch_transformer = coordsbatch_transformer
        self.sampleweights_coordstovals = sampleweights_coordstovals
        self.sampleweights_from_inputstargets =\
            sampleweights_from_inputstargets
        if sampleweights_coordstovals is not None:
            assert sampleweights_from_inputstargets is None
        if sampleweights_from_inputstargets is not None:
            assert sampleweights_coordstovals is None
        self.qc_func = qc_func
 
    def __getitem__(self, index):
        coords_batch = self.coordsbatch_producer[index]
        if (self.coordsbatch_transformer is not None):
            coords_batch = self.coordsbatch_transformer(coords_batch)
        inputs = self.inputs_coordstovals(coords_batch)
        if (self.targets_coordstovals is not None):
            targets = self.targets_coordstovals(coords_batch)
        else:
            targets=None
        if (self.qc_func is not None):
            qc_mask = self.qc_func(inputs=inputs, targets=targets)
            inputs = apply_mask(tomask=inputs, mask=qc_mask)
            if (targets is not None):
                targets = apply_mask(tomask=targets, mask=qc_mask)
        else:
            qc_mask = None
        if (self.sampleweights_coordstovals is not None):
            sample_weights = self.sampleweights_coordstovals(coords_batch)
            return (inputs, targets, sample_weights)
        elif (self.sampleweights_from_inputstargets is not None):
            sample_weights = self.sampleweights_from_inputstargets(
                                inputs=inputs, targets=targets)
            return (inputs, targets, sample_weights)
        else:
            if (self.targets_coordstovals is not None):
                return (inputs, targets)
            else:
                return inputs
   
    def __len__(self):
        return len(self.coordsbatch_producer)
    
    def on_epoch_end(self):
        self.coordsbatch_producer.on_epoch_end()
