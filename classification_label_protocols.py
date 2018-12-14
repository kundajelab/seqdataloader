from math import floor,ceil
import pdb 
def label_ambiguous_bins(chrom,task_name,non_zero_bins,first_overlap_seq_start,last_overlap_seq_start,seq_size,args):
    left_ambiguous_start= first_overlap_seq_start-args.bin_stride
    left_ambiguous_end=left_ambiguous_start+seq_size
    left_ambiguous_seq=(chrom,left_ambiguous_start,left_ambiguous_end)

    right_ambiguous_start=last_overlap_seq_start+args.bin_stride
    right_ambiguous_end=right_ambiguous_start+seq_size
    right_ambiguous_seq=(chrom,right_ambiguous_start,right_ambiguous_end)
    
    if left_ambiguous_seq not in non_zero_bins:
        non_zero_bins[left_ambiguous_seq]=dict()
    non_zero_bins[left_ambiguous_seq][task_name]='-1'
    
    if right_ambiguous_seq not in non_zero_bins:
        non_zero_bins[right_ambiguous_seq]=dict()
    non_zero_bins[right_ambiguous_seq][task_name]='-1'
    return non_zero_bins

def peak_summit_in_bin(task_name,task_bed,args):
    '''
    For each peak, the summit position is determined. 

    The minimum bin with a positive label is args.binsize upstream of the summit;
    The max bin with a positive label is args.binsize downstream of the summit 

    Within this range, bin centers are shifted by args.bin_stride 

    If specified in args.allow_ambiguous, then adjacent bins to the two extremes are marked as 
    ambiguous 
    '''
    non_zero_bins=dict()
    seq_size=args.bin_size+args.left_flank+args.right_flank
    
    
    for entry in task_bed:
        chrom=entry[0]
        peak_start=int(entry[1])
        summit=peak_start+int(entry[-1])        
        min_bin_start=ceil((summit-args.bin_size)/args.bin_stride)*args.bin_stride
        max_bin_start=floor(summit/args.bin_stride)*args.bin_stride 

        #padding to the left and right of the binsize region needed to generated the desired bin_size 
        #label all bins with center regions between min_bin_start and max_bin_start (inclusive) with 1 
        for bin_start in range(min_bin_start,max_bin_start+1,args.bin_stride):
            cur_seq_start=bin_start-args.left_flank
            cur_seq_end=cur_seq_start+seq_size
            cur_seq=(chrom,cur_seq_start,cur_seq_end)
            non_zero_bins[cur_seq]=dict()
            non_zero_bins[cur_seq][task_name]='1'

        #if user specified use of ambiguous bins,
        #label the adjacent bins to min_bin_start and max_bin_start
        #as ambiguous
        if (args.allow_ambiguous==True):
            min_seq_start=min_bin_start-args.left_flank
            max_seq_start=max_bin_start-args.left_flank
            non_zero_bins=label_ambiguous_bins(chrom,task_name,non_zero_bins,min_seq_start,max_seq_start,seq_size,args)
    return non_zero_bins

def peak_percent_overlap_with_bin(task_name,task_bed,args):
    '''
    50% of the central 200bp region in a 1kb bin must overlap with the peak 
    '''
    non_zero_bins=dict()
    seq_size=args.bin_size+args.left_flank+args.right_flank
    for entry in task_bed:
        chrom=entry[0]
        peak_start=int(entry[1])
        peak_end=int(entry[2])
        min_overlap=args.overlap_thresh*args.bin_size
        min_bin_start=peak_start-min_overlap
        max_bin_start=peak_end-min_overlap
        first_overlap_seq_start=int(ceil((min_bin_start-args.left_flank)/args.bin_stride))*args.bin_stride
        last_overlap_seq_start=int(floor((max_bin_start-args.left_flank)/args.bin_stride))*args.bin_stride 
        peak_length=peak_end-peak_start+1
        for seq_start in range(first_overlap_seq_start,last_overlap_seq_start+1,args.bin_stride ):
            cur_seq=(chrom,seq_start,seq_start+seq_size)
            non_zero_bins[cur_seq]=dict()
            non_zero_bins[cur_seq][task_name]='1'

        if(args.allow_ambiguous==True): 
            non_zero_bins=label_ambiguous_bins(chrom,task_name,non_zero_bins,first_overlap_seq_start,last_overlap_seq_start,seq_size,args)
    return non_zero_bins

