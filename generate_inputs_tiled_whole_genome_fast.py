import argparse
import pybedtools
from pybedtools import BedTool
import pandas as pd
import numpy as np 
import pdb
import csv
from multiprocessing.pool import ThreadPool

def parse_args():
    parser=argparse.ArgumentParser("generated task-specific negatives for a multi-tasked input matrix")
    parser.add_argument("--task_list",help="this is a tab-separated file with the name of the task in the first column and the path to the corresponding bed file in the second column")
    parser.add_argument("--outf")
    parser.add_argument("--bins")
    parser.add_argument("--training_bin_size",type=int,default=1000)
    parser.add_argument("--threads",type=int,default=1)
    parser.add_argument("--overlap_thresh",type=float,default=0.5,help="minimum percent of bin that must overlap a peak for a positive label")
    return parser.parse_args()


def get_labels_one_task(task_bedfile,master_bed,training_bin_size,overlap_thresh,task,outfname):
    labels=[]
    processed=dict()
    print('thread!')
    try:
        intersection=master_bed.intersect(task_bedfile,wao=True)
        print("Got intersection!") 
    except:
        return None
    overlap_check=training_bin_size*overlap_thresh
    outf=open(outfname,'w')
    outf.write(task+'\n')        
    for entry in intersection:
                    
        master_peak='\t'.join(entry[0:3])
        if master_peak in processed:
            continue

        overlap=entry[-1]
        if (overlap > overlap_check ):
            outf.write('1\n')
        elif overlap>0:
            outf.write('-1\n')
        else:
            outf.write('0\n')
        processed[master_peak]=1
    #write labels to output file
    print("done writing  outputs:"+task+".labels")
    return
    
def get_labels(task_bedfiles,master_bed,training_bin_size,threads,overlap_thresh,outf):
    pool=ThreadPool(threads)    
    print("getting labels with " +str(threads)+ " threads" )
    label_dict=dict()
    for task in task_bedfiles:
        print(task)
        task_bedfile=task_bedfiles[task]
        outfname=outf+'.task'
        pool.apply_async(get_labels_one_task,args=(task_bedfile,master_bed,training_bin_size,overlap_thresh,task,outfname))
    pool.close()
    pool.join()                                           
    print("got labels") 
    


def load_bed_files_for_tasks(tasks): 
    task_bedfiles=dict()
    for index,row in tasks.iterrows():
        #pdb.set_trace() 
        taskname=row['Task']
        print(taskname)
        cur_bed=BedTool(row['File'])
        task_bedfiles[taskname]=cur_bed
    print("loaded bed files for all tasks")
    return task_bedfiles

def main():
    args=parse_args()
    tasks=pd.read_csv(args.task_list,sep='\t',header=0)

    #load the genome-wide intervals
    master_bed=BedTool(args.bins)
    print("loaded genome-wide bins")
    
    #load the task bed files and generate a master bed file 
    task_bedfiles=load_bed_files_for_tasks(tasks)
    
    #generate labels for each task
    label_dict=get_labels(task_bedfiles,master_bed,args.training_bin_size,args.threads,args.overlap_thresh,args.outf)
        
    

if __name__=="__main__":
    main()
    
    
