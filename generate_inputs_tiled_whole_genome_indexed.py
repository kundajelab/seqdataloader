import argparse
from pybedtools import BedTool
import pyBigWig 
import pandas as pd
import numpy as np 
import pdb
import csv
from classification_label_protocols import *
from regression_label_protocols import * 
from multiprocessing.pool import ThreadPool
from utils import merge_dictionaries
import gzip 

#Approaches to determining classification labels
#Others can be added here (imported from classification_label_protocols) 
labeling_approaches={
    "peak_summit_in_bin_classification":peak_summit_in_bin_classification,
    "peak_percent_overlap_with_bin_classification":peak_percent_overlap_with_bin_classification,
    "peak_summit_in_bin_regression":peak_summit_in_bin_regression,
    "peak_percent_overlap_with_bin_regression":peak_percent_overlap_with_bin_regression,
    "all_genome_bins_regression":all_genome_bins_regression
    }

def parse_args():
    parser=argparse.ArgumentParser(description="Generate genome-wide labeled bins for a set of narrowPeak task files ")
    parser.add_argument("--task_list",help="this is a tab-separated file with the name of the task in the first column, the path to the corresponding narrowPeak(.gz) file in the second column (optionally), and the path to the corresponding bigWig file in the third column (optionally, for regression)")
    parser.add_argument("--out_bed",help="output filename that labeled bed file will be saved to.")
    parser.add_argument("--chrom_sizes",help="chromsizes file for the reference genome. First column is chrom name; second column is chrom size")
    parser.add_argument("--bin_stride",type=int,default=50,help="bin_stride to shift adjacent bins by")
    parser.add_argument("--left_flank",type=int,default=400,help="left flank")
    parser.add_argument("--right_flank",type=int,default=400,help="right flank")
    parser.add_argument("--bin_size",type=int,default=200,help="flank around bin center where peak summit falls in a positivei bin")
    parser.add_argument("--threads",type=int,default=1)
    parser.add_argument("--subthreads",type=int,default=20,help="This is only useful for regression labels for each bin in the genome. Each task-thread will generate 20 subthreads to allow for parallel processing of chromosomes. Reduce number to use fewer threads")
    parser.add_argument("--overlap_thresh",type=float,default=0.5,help="minimum percent of bin that must overlap a peak for a positive label")
    parser.add_argument("--allow_ambiguous",default=False,action="store_true")
    parser.add_argument("--labeling_approach",choices=["peak_summit_in_bin_classification",
                                                       "peak_percent_overlap_with_bin_classification",
                                                       "peak_summit_in_bin_regression",
                                                       "peak_percent_overlap_with_bin_regression",
                                                       "all_genome_bins_regression"])    
    return parser.parse_args()

def get_labels_one_task(task_name,task_bed,task_bigwig,args):
    #determine the appropriate labeling approach
    return labeling_approaches[args.labeling_approach](task_name,task_bed,task_bigwig,args)
    

def write_output_bed(args,task_names,non_zero_bins):
    '''
    Generate output file 
    args - arguments that were passed to the main script 
    non_zero_bins - dictionary of the form tuple(chrom,start_bin,end_bin)->task_name->label 
    task_names - list of unique task names 
    '''
    #open the output file and write the header 
    outf=gzip.open(args.out_bed,'wt')
    header='\t'.join(['\t'.join(['Chrom','Start','End']),
                      '\t'.join(task_names)])
    outf.write(header+'\n')

    #load the chromosome sizes 
    chrom_sizes=pd.read_table(args.chrom_sizes,header=None,sep='\t')
    print("loaded chrom_sizes")
    
    num_tasks=len(task_names)
    seq_size=args.bin_size+args.left_flank+args.right_flank 
    #iterate through each bin in the genome 
    for index,row in chrom_sizes.iterrows():
        chrom=row[0]
        chrom_size=int(row[1])
        
        print("Writing output file entries for chrom:"+str(chrom))
        for seq_start in range(seq_size,chrom_size-seq_size,args.bin_stride):
            
            #store the current sequence as a tuple
            seq_end=seq_start+seq_size
            cur_seq=tuple([chrom,seq_start,seq_end])
            if cur_seq not in non_zero_bins:
                
                #all tasks have 0 label
                outf.write('\t'.join(['\t'.join([str(i) for i in cur_seq]),'\t'.join(['0']*num_tasks)])+'\n')
            else:
                outf.write('\t'.join([str(i) for i in cur_seq]))
                #iterate through tasks to determine appropriate labels
                for task_name in task_names:
                    if task_name in non_zero_bins[cur_seq]:
                        outf.write('\t'+str(non_zero_bins[cur_seq][task_name]))
                    else:
                        outf.write('\t0')
                outf.write('\n')
    outf.close()
    
def get_nonzero_bins(args,tasks):
    #parallelized bin labeling
    pool=ThreadPool(args.threads)
    if args.allow_ambiguous==True:
        print(' '.join(["determining positive and ambiguous bins with", str(args.threads), "threads"]))
    else:
        print(' '.join(["determining positive bins with",str(args.threads),"threads"]))
              
    #create list of dictionaries to store non-zero bin labels of the form:
    # tuple(chrom,bin_start,bin_end)-> task->label
    non_zero_bins_list=[]

    #keep track of all task names
    task_names=[]
    for index,row in tasks.iterrows():
        task_name=row[0]
        task_names.append(task_name)
        #try to get the peak file associated with the task (if it's provided)
        try:
            task_bed=BedTool(row[1])
        except:
            task_bed=None
            print("Warning! No peak file was provided for task:"+task_name+"; Is this intentional?")
        #try to get the bigWig file associated with the task (if it's provided)
        try:
            task_bigwig=pyBigWig.open(row[2])
        except:
            print("Warning! No bigWig file was provided for task:"+task_name+"; Is this intentional?")
            task_bigwig=None
            
        #get non-zero bin labels for the current task 
        non_zero_bins_list.append(pool.apply_async(get_labels_one_task,args=(task_name,task_bed,task_bigwig,args)))  
    pool.close()
    pool.join()
    print("finished parsing and labeling task bed files")
    
    #collapse the non_zero_bins list of dictionaries into a single dictionary
    print("merging non-zero label dictionaries")
    non_zero_bins_dict=non_zero_bins_list[0].get()
    for non_zero_bins_subdict in non_zero_bins_list[1::]: 
        non_zero_bins_dict=merge_dictionaries(non_zero_bins_dict,non_zero_bins_subdict.get())
    return task_names,non_zero_bins_dict

def main():
    
    #parse the input arguments
    args=parse_args()

    #read in the metadata file with:
    #task names in column 1,
    #path to peak file in column 2,
    #path to bigWig file in column 3
    tasks=pd.read_csv(args.task_list,sep='\t',header=None)

    #multi-threaded identification of non-zero bin labels for each task 
    task_names,non_zero_bins_dict=get_nonzero_bins(args,tasks) 
    
    #write the output file
    write_output_bed(args,task_names,non_zero_bins_dict)
    

if __name__=="__main__":
    main()
    
    
