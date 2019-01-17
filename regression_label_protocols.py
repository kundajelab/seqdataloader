from math import floor,ceil
import pandas as pd
from multiprocessing.pool import ThreadPool
from utils import merge_dictionaries

import pdb 
def label_ambiguous_bins_regression(chrom,task_name,non_zero_bins,min_bin_start,max_bin_start,seq_size,task_bigwig,args):
    left_ambiguous_bin_start=min_bin_start-args.bin_stride
    left_ambiguous_bin_end=left_ambiguous_bin_start+args.bin_size
    left_ambiguous_seq_start=left_ambiguous_bin_start-args.left_flank
    left_ambiguous_seq_end=left_ambiguous_seq_start+seq_size 
    left_ambiguous_seq=(chrom,left_ambiguous_seq_start,left_ambiguous_seq_end)

    try:
        left_ambiguous_coverage=round(task_bigwig.stats(chrom,left_ambiguous_bin_start,left_ambiguous_bin_end)[0],3)
        if left_ambiguous_seq not in non_zero_bins:
            non_zero_bins[left_ambiguous_seq]=dict()
        non_zero_bins[left_ambiguous_seq][task_name]=left_ambiguous_coverage
    except:
        print("skipping left flank:"+str(left_ambiguous_seq))

    right_ambiguous_bin_start=max_bin_start+args.bin_stride
    right_ambiguous_bin_end=right_ambiguous_bin_start+args.bin_size
    right_ambiguous_seq_start=right_ambiguous_bin_start-args.right_flank
    right_ambiguous_seq_end=right_ambiguous_seq_start+seq_size 
    right_ambiguous_seq=(chrom,right_ambiguous_seq_start,right_ambiguous_seq_end)

    try:
        right_ambiguous_coverage=round(task_bigwig.stats(chrom,right_ambiguous_bin_start,right_ambiguous_bin_end)[0],3)
        if right_ambiguous_seq not in non_zero_bins:
            non_zero_bins[right_ambiguous_seq]=dict()
        non_zero_bins[right_ambiguous_seq][task_name]=right_ambiguous_coverage 
    except:
        print("skipping right flank:"+str(right_ambiguous_seq))
    return non_zero_bins

def peak_summit_in_bin_regression(task_name,task_bed,task_bigwig,args):
    '''
    For each peak, the summit position is determined. 

    The minimum bin with bedtools coverage is args.binsize upstream of the summit;
    The max bin with bedtools coverage is args.binsize downstream of the summit 

    Within this range, bin centers are shifted by args.bin_stride 

    If specified in args.allow_ambiguous, then coverage is also computed in adjacent bins to the two extremes are marked as 
    ambiguous 
    '''
    print(task_name)
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
            #get bigWig coverage in the current bin (ignoring the flanks)
            try:
                non_zero_bins[cur_seq][task_name]=round(task_bigwig.stats(chrom,bin_start,bin_start+args.bin_size)[0],3)
            except:
                print("skipping:"+str(cur_seq))
        #if user specified use of ambiguous bins,
        #label the adjacent bins to min_bin_start and max_bin_start
        #as ambiguous
        if (args.allow_ambiguous==True):
            non_zero_bins=label_ambiguous_bins_regression(chrom,task_name,non_zero_bins,min_bin_start,max_bin_start,seq_size,task_bigwig,args)
    return non_zero_bins

def peak_percent_overlap_with_bin_regression(task_name,task_bed,task_bigwig,args):
    '''
    50% of the central 200bp region in a 1kb bin must overlap with the peak for coverage to be computed in the provided bigWig 
    '''
    print(task_name)
    non_zero_bins=dict()
    seq_size=args.bin_size+args.left_flank+args.right_flank
    for entry in task_bed:
        chrom=entry[0]
        peak_start=int(entry[1])
        peak_end=int(entry[2])
        min_overlap=int(round(args.overlap_thresh*args.bin_size))
        min_bin_start=peak_start-min_overlap
        max_bin_start=peak_end-min_overlap
        first_overlap_seq_start=int(ceil((min_bin_start-args.left_flank)/args.bin_stride))*args.bin_stride
        last_overlap_seq_start=int(floor((max_bin_start-args.left_flank)/args.bin_stride))*args.bin_stride 
        peak_length=peak_end-peak_start+1
        for seq_start in range(first_overlap_seq_start,last_overlap_seq_start+1,args.bin_stride ):
            cur_seq=(chrom,seq_start,seq_start+seq_size)
            non_zero_bins[cur_seq]=dict()
            try:
                non_zero_bins[cur_seq][task_name]=round(task_bigwig.stats(chrom,seq_start+args.left_flank,seq_start+seq_size-args.right_flank)[0],3)
            except:
                print("skipping:"+str(cur_seq))
        if(args.allow_ambiguous==True): 
            non_zero_bins=label_ambiguous_bins_regression(chrom,task_name,non_zero_bins,min_bin_start,max_bin_start,seq_size,task_bigwig,args)
    return non_zero_bins

def one_chrom_genome_bins_regression(task_name,task_bigwig,chrom,chrom_size,seq_size,args):
    #make sure we don't run off the edges of the chromosome when we add left/right flanks to the bins
    print("getting coverage from chrom:"+str(chrom))
    non_zero_bins=dict()
    counter=0
    for bin_start in range(seq_size,chrom_size-seq_size,args.bin_stride):
        counter+=1
        if counter % 10000 == 0:
            print(chrom+":"+str(bin_start))
        bin_end=bin_start+args.bin_size
        #add left and right flanks to get the sequence start and end coordinates 
        seq_start=bin_start-args.left_flank
        seq_end=bin_end+args.right_flank
        cur_seq=(chrom,seq_start,seq_end)

        #compute bigWig coverage for the bin
        try:
            bin_coverage=round(task_bigwig.stats(chrom,bin_start,bin_end)[0],3)                
            #store bin coverage entry in the dictionary 
            non_zero_bins[cur_seq]=dict()
            non_zero_bins[cur_seq][task_name]=bin_coverage
        except:
            print("skipping:"+str(cur_seq))
    return non_zero_bins

def all_genome_bins_regression(task_name,task_bed,task_bigwig,args):
    '''
    compute bigWig coverage for all bins in the genome, regardless of whether a called peak overlaps the bin
    Warning: This will take a long time! 
    '''
    print(task_name)
    pool=ThreadPool(args.subthreads)
    non_zero_bins_list=[] 
    chrom_sizes=pd.read_table(args.chrom_sizes,header=None,sep='\t')
    seq_size=args.bin_size+args.left_flank+args.right_flank
    #iterate through each bin in the genome    
    for index,row in chrom_sizes.iterrows():
        chrom=row[0]
        chrom_size=int(row[1])
        non_zero_bins_list.append(pool.apply_async(one_chrom_genome_bins_regression,args=(task_name,task_bigwig,chrom,chrom_size,seq_size,args)))
    pool.close()
    pool.join()
    print("merging across chromosomes")
    non_zero_bins=non_zero_bins_list[0].get()
    for non_zero_bins_subdict in non_zero_bins_list[1::]:
        non_zero_bins.update(non_zero_bins_subdict.get())
    return non_zero_bins

