from math import floor,ceil
import pdb 
def label_ambiguous_bins(chrom,task_name,non_zero_bins,min_bin_start,max_bin_start,args):
    left_ambiguous_start= min_bin_start-args.stride
    left_ambiguous_end=left_ambiguous_start+args.bin_size
    left_ambiguous_bin=(chrom,left_ambiguous_start,left_ambiguous_end)

    right_ambiguous_start=max_bin_start+args.stride
    right_ambiguous_end=right_ambiguous_start+args.bin_size    
    right_ambiguous_bin=(chrom,right_ambiguous_start,right_ambiguous_end)
    
    if left_ambiguous_bin not in non_zero_bins:
        non_zero_bins[left_ambiguous_bin]=dict()
    non_zero_bins[left_ambiguous_bin][task_name]='-1'
    
    if right_ambiguous_bin not in non_zero_bins:
        non_zero_bins[right_ambiguous_bin]=dict()
    non_zero_bins[right_ambiguous_bin][task_name]='-1'
    return non_zero_bins

def peak_summit_near_bin_center(task_name,task_bed,args):
    '''
    For each peak, the summit position is determined. 

    The minimum bin with a positive label is args.bin_center_size upstream of the summit;
    The max bin with a positive label is args.bin_center_size downstream of the summit 

    Within this range, bin centers are shifted by args.stride 

    If specified in args.allow_ambiguous, then adjacent bins to the two extremes are marked as 
    ambiguous 
    '''
    non_zero_bins=dict() 
    for entry in task_bed:
        chrom=entry[0]
        peak_start=int(entry[1])
        summit=peak_start+int(entry[-1])
        min_bin_center_start=ceil((summit-args.bin_center_size)/args.stride)*args.stride
        max_bin_center_start=floor(summit/args.stride)*args.stride 

        #padding to the left and right of the bin_center_size region needed to generated the desired bin_size 
        topad=int(round((args.bin_size-args.bin_center_size)/2))
        #label all bins with center regions between min_bin_center_start and max_bin_center_start (inclusive) with 1 
        for bin_center_start in range(min_bin_center_start,max_bin_center_start+1,args.stride):
            cur_bin_start=bin_center_start-topad
            cur_bin_end=cur_bin_start+args.bin_size 
            cur_bin=(chrom,cur_bin_start,cur_bin_end)
            non_zero_bins[cur_bin]=dict()
            non_zero_bins[cur_bin][task_name]='1'

        #if user specified use of ambiguous bins,
        #label the adjacent bins to min_bin_center_start and max_bin_center_start
        #as ambiguous
        if (args.allow_ambiguous==True):
            min_bin_start=min_bin_center_start-topad
            max_bin_start=max_bin_center_start-topad
            non_zero_bins=label_ambiguous_bins(chrom,task_name,non_zero_bins,min_bin_start,max_bin_start,args)
    return non_zero_bins

def peak_percent_overlap_with_bin(task_name,task_bed,args):
    non_zero_bins=dict()
    for entry in task_bed:
        chrom=entry[0]
        peak_start=int(entry[1])
        peak_end=int(entry[2])

        peak_length=peak_end-peak_start+1
        min_overlap=int(round(min([peak_length,args.bin_size])*args.overlap_thresh))

        overlap_bin_first_start=int(ceil(peak_start-min_overlap)/args.stride)*args.stride
        overlap_bin_last_start=int(floor(peak_end-min_overlap)/args.stride)*args.stride
        for bin_start in range(overlap_bin_first_start,overlap_bin_last_start+1,args.stride ):
            cur_bin=(chrom,bin_start,bin_start+args.bin_size)
            non_zero_bins[cur_bin]=dict()
            non_zero_bins[cur_bin][task_name]='1'

        if(args.allow_ambiguous==True): 
            non_zero_bins=label_ambiguous_bins(chrom,task_name,non_zero_bins,overlap_bin_first_start,overlap_bin_last_start,args)
    return non_zero_bins

