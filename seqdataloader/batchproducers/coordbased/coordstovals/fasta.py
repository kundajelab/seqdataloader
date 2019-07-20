from __future__ import division, print_function, absolute_import
import numpy as np
from pyfaidx import Fasta
from .core import AbstractSingleNdarrayCoordsToVals


class PyfaidxCoordsToVals(AbstractSingleNdarrayCoordsToVals):

    def __init__(self, genome_fasta_path, center_size_to_use=None, **kwargs):
        """
        Args:
            genome_fasta_path (:obj:`str`): path to the genome .fa file
            **kwargs: arguments for :obj:`AbstractSingleNdarrayCoordsToVals`
        """
        super(PyfaidxCoordsToVals, self).__init__(**kwargs)
        self.center_size_to_use = center_size_to_use
        self.genome_fasta = genome_fasta_path
        self.genome_object = Fasta(genome_fasta_path)
        self.ltrdict = {
           'a':[1,0,0,0],'c':[0,1,0,0],'g':[0,0,1,0],'t':[0,0,0,1],
           'n':[0,0,0,0],'A':[1,0,0,0],'C':[0,1,0,0],'G':[0,0,1,0],
           'T':[0,0,0,1],'N':[0,0,0,0]}
        #use [0,0,0,0] as default if FASTA base is not found in the dictionary
        self.onehot_encoder = (
            lambda seq: np.array([self.ltrdict.get(x,[0,0,0,0]) for x in seq]))
    
    def _get_ndarray(self, coors):
        """
        Args:
            coors (:obj:`list` of :obj:`Coordinates): if
                center_size_to_use is not specified, all the
                coordinates must be of the same length
            
        Returns:
            numpy ndarray of dims (nexamples x width x 4)
        """
        seqs = []
        for coor in coors:
          if (self.center_size_to_use is not None):
            the_center = int((coor.start + coor.end)*0.5)
            seqs.append(self.genome_object[coor.chrom][
                the_center-int(0.5*self.center_size_to_use):
                the_center+(self.center_size_to_use
                            -int(0.5*self.center_size_to_use))])
          else:
            seqs.append(self.genome_object[coor.chrom][coor.start:coor.end])

        onehot_seqs = []
        for seq,coor in zip(seqs, coors):
            onehot = self.onehot_encoder(seq=seq.seq)
            if (coor.isplusstrand==False):
                onehot = onehot[::-1, ::-1]
            onehot_seqs.append(onehot)
        lengths = set([len(x) for x in onehot_seqs])
        assert len(lengths)==1, ("All the sequences must be of the same"
            +"lengths, but lengths are "+str(lengths))
        return np.array(onehot_seqs)
