from __future__ import division, print_function, absolute_import
import numpy as np
import pyBigWig
from .core import CoordsToVals


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
        pos_profile_values = self.pos_strand_reader.read_values(coors=coors)
        neg_profile_values = np.abs(
            self.neg_strand_reader.read_values(coors=coors))
        orig_len = pos_profile_values.shape[1]
        left_start = int(orig_len/2) - int(self.center_size_to_use/2)
        right_end = int(orig_len/2) + (self.center_size_to_use
                                       - int(self.center_size_to_use/2))
        pos_profile_values = pos_profile_values[:, left_start:right_end]
        neg_profile_values = neg_profile_values[:, left_start:right_end]
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


def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
      

def smooth(profile_array, smoothing_windows):
    smoothed_profiles = []
    for smoothing_window in smoothing_windows:
      strided_profile = rolling_window(a=profile_array,
                                             window=smoothing_window)
      smoothed_profile_nopad = np.mean(strided_profile, axis=-1)
      leftpadlen = int((smoothing_window-1)/2)
      rightpadlen =\
                (smoothing_window-1)-int((smoothing_window-1)/2)
      padded_profile = np.pad(
                array=smoothed_profile_nopad,
                pad_width=((0,0),(leftpadlen, rightpadlen)),
                mode='constant')
      smoothed_profiles.append(padded_profile[:,:,None])
  return smoothed_profiles
    

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
            strided_profile = rolling_window(a=profile_sum,
                                             window=smoothing_window)
            smoothed_profile_nopad = np.mean(strided_profile, axis=-1)
            leftpadlen = int((smoothing_window-1)/2)
            rightpadlen =\
                (smoothing_window-1)-int((smoothing_window-1)/2)
            padded_profile = np.pad(
                array=smoothed_profile_nopad,
                pad_width=((0,0),(leftpadlen, rightpadlen)),
                mode='constant')
            
            #print(padded_profile.shape)
            
            smoothed_profiles.append(padded_profile[:,:,None])
        
        smooth_profiles = np.concatenate(smoothed_profiles, axis=2)
        
        return (np.log(pos_counts+neg_counts+1),
                smooth_profiles)  