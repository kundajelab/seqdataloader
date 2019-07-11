from __future__ import division, print_function, absolute_import
import numpy as np
import pyBigWig
from .core import CoordsToVals, get_new_coors_around_center
from ..core import Coordinates


def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def smooth_profiles(profiles, smoothing_window):
    assert len(profiles.shape)==3
    leftpadlen = int((smoothing_window-1)/2)
    rightpadlen =\
        (smoothing_window-1)-int((smoothing_window-1)/2)
    padded_profiles = np.pad(
        array=profiles,
        pad_width=((0,0),(leftpadlen, rightpadlen), (0,0)),
        mode='edge')
    smoothed_profiles = np.mean(rolling_window(
                        a=padded_profiles.transpose(0,2,1),
                        window=smoothing_window), axis=-1).transpose((0,2,1))
    return smoothed_profiles


class BigWigReader(object):

    def __init__(self, bigwig_path):
        """
        Args:
            bigwig_path (:obj:`str`): path to the .bw file
        """
        self.bigwig_path = bigwig_path
        self.bw = pyBigWig.open(bigwig_path)
        
    def read_values(self, coors):
        """
        Args:
            coords (:obj:`list` of :obj:Coordinates)
        
        Returns:
            ndarray of dims (nexamples x width). All the coordinates must be
                of the same length.
        """
        to_return = []
        for coor in coors:
            to_append = np.nan_to_num(
                          x=self.bw.values(coor.chrom, coor.start, coor.end))
            if (coor.isplusstrand==False):
                to_append = to_append[::-1]
            to_return.append(to_append)
        lengths = set([len(x) for x in to_return])
        assert len(lengths)==1, ("All the sequences must be of the same"
            +"lengths, but lengths are "+str(lengths))
        return np.array(to_return)


class LogCountsAndProfile(CoordsToVals):

    def __init__(self, bigwig_path, counts_mode_name,
                       profile_mode_name):
        self.reader = BigWigReader(bigwig_path=bigwig_path)
        self.counts_mode_name = counts_mode_name
        self.profile_mode_name = profile_mode_name
    
    def __call__(self, coors):
        profile_values = self.reader.read_values(coors=coors)
        counts = np.log(np.sum(profile_values, axis=-1)+1)
        to_return = {self.counts_mode_name: counts,
                     self.profile_mode_name: profile_values}
        return to_return


class AbstractCountAndProfileTransformer(object):

    def __call__(self, counts, profiles):
        raise NotImplementedError()

    def chain(self, count_and_profile_transformer):
        def chained_count_and_profile_transformer(counts, profiles):
            counts, profiles = self(counts=counts, profiles=profiles)
            return count_and_profile_transformer(
                    counts=counts, profiles=profiles)
        return chained_count_and_profile_transformer


class LogCountsPlusOne(AbstractCountAndProfileTransformer):

    def __call__(self, counts, profiles):
        return np.log(counts+1), profiles


class SmoothProfiles(AbstractCountAndProfileTransformer):

    def __init__(self, smoothing_windows):
        self.smoothing_windows = smoothing_windows

    def __call__(self, counts, profiles):
        profiles_to_return = np.concatenate([
           smooth_profiles(profiles=profiles, smoothing_window=x)
           for x in self.smoothing_windows], axis=-1)
        return counts, profiles_to_return


class MultiTrackCountsAndProfile(CoordsToVals):

    def __init__(self, bigwig_paths,
                       counts_and_profiles_transformer,
                       counts_mode_name,
                       profile_mode_name, center_size_to_use):
        self.bigwig_readers = [BigWigReader(bigwig_path=x)
                               for x in bigwig_paths]
        self.counts_and_profiles_transformer = counts_and_profiles_transformer
        self.counts_mode_name = counts_mode_name
        self.profile_mode_name = profile_mode_name
        self.center_size_to_use = center_size_to_use

    def _get_counts_and_vals(self, coors):
        new_coors = get_new_coors_around_center(
                            coors=coors,
                            center_size_to_use=self.center_size_to_use)
        #concatenate the results of the bigwig readers along the last axis
        profiles = np.concatenate([
                          x.read_values(coors=new_coors)[:,:,None]
                          for x in self.bigwig_readers], axis=-1)
        counts = np.sum(profiles, axis=1)
        return (counts, profiles)
    
    def __call__(self, coors):
        counts, profiles = self._get_counts_and_vals(coors=coors)
        counts_transformed, profile_transformed =\
            self.counts_and_profiles_transformer(
                  counts=counts,
                  profiles=profiles)
        return {self.counts_mode_name: counts_transformed,
                self.profile_mode_name: profile_transformed}
 

class AbstractPosAndNegStrandCountsAndProfile(CoordsToVals):

    def __init__(self, pos_strand_bigwig_path, neg_strand_bigwig_path,
                       counts_mode_name, profile_mode_name,
                       center_size_to_use):
        self.pos_strand_reader =\
          BigWigReader(bigwig_path=pos_strand_bigwig_path)
        self.neg_strand_reader =\
          BigWigReader(bigwig_path=neg_strand_bigwig_path)
        self.counts_mode_name = counts_mode_name
        self.profile_mode_name = profile_mode_name
        self.center_size_to_use = center_size_to_use
        
    def _get_pos_and_neg_counts_and_vals(self, coors):
        new_coors = get_new_coors_around_center(
                            coors=coors,
                            center_size_to_use=self.center_size_to_use)
        first_strand_profile_values = self.pos_strand_reader.read_values(
                                  coors=new_coors)
        second_strand_profile_values = np.abs(
            self.neg_strand_reader.read_values(coors=new_coors))
        pos_profile_values = []
        neg_profile_values = []
        #need to swap the pos and neg strands if the strand is negative
        for (first_strand, second_strand, coor) in zip(first_strand_profile_values,
                                                       second_strand_profile_values,
                                                       coors):
            if (coor.isplusstrand==True):
                pos_profile_values.append(first_strand)
                neg_profile_values.append(second_strand)
            else:
                pos_profile_values.append(second_strand)
                neg_profile_values.append(first_strand)
        pos_profile_values = np.array(pos_profile_values)
        neg_profile_values = np.array(neg_profile_values)        
        pos_counts = np.sum(pos_profile_values, axis=-1)
        neg_counts = np.sum(neg_profile_values, axis=-1)
        return (pos_counts, neg_counts, pos_profile_values, neg_profile_values)
    
    """
    Returns:
        ndarray: combined/transformed counts
        ndarray: combined/transformed profile
    """
    def combine_pos_and_neg_counts_and_vals(self,
      pos_counts, neg_counts, pos_profile_values, neg_profile_values):
        raise NotImplementedError()
    
    def __call__(self, coors):
        pos_counts, neg_counts, pos_profile_values, neg_profile_values =(
            self._get_pos_and_neg_counts_and_vals(coors=coors))
        counts_ndarray, profile_ndarray =\
          self.combine_pos_and_neg_counts_and_vals( pos_counts=pos_counts, neg_counts=neg_counts,
            pos_profile_values=pos_profile_values,
            neg_profile_values=neg_profile_values)
        return {self.counts_mode_name: counts_ndarray,
                self.profile_mode_name: profile_ndarray}


class PosAndNegSeparateLogCounts(AbstractPosAndNegStrandCountsAndProfile):
  
    def __init__(self, **kwargs):
        super(PosAndNegSeparateLogCounts,self).__init__(**kwargs)
  
    def combine_pos_and_neg_counts_and_vals(self,
      pos_counts, neg_counts, pos_profile_values, neg_profile_values):
      
        return (np.concatenate([np.log(pos_counts+1)[:,None], 
                                np.log(neg_counts+1)[:,None]], axis=1),
                np.concatenate(
                    [pos_profile_values[:,:,None],
                     neg_profile_values[:,:,None]], axis=2))
    

class PosAndNegSmoothWindowCollapsedLogCounts(
        AbstractPosAndNegStrandCountsAndProfile):
  
    def __init__(self, smoothing_windows, **kwargs):
        super(PosAndNegSmoothWindowCollapsedLogCounts, self).__init__(**kwargs)
        self.smoothing_windows = smoothing_windows
  
    def combine_pos_and_neg_counts_and_vals(self, pos_counts, neg_counts,
        pos_profile_values, neg_profile_values):
        
        profile_sum = (
            pos_profile_values[:,:]+
            neg_profile_values[:,:])
        
        smoothed_profiles = []
        for smoothing_window in self.smoothing_windows:
            padded_profile = smooth_profiles(profiles=profile_sum[:,:,None],
                                             smoothing_window=smoothing_window) 
            smoothed_profiles.append(padded_profile)
        
        smoothed_profiles = np.concatenate(smoothed_profiles, axis=2)
        
        return (np.log(pos_counts+neg_counts+1), smoothed_profiles)  
