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
    parser.add_argument("--training_bin_size")
    parser.add_argument("--threads",type=int,default=1) 
    return parser.parse_args()

def write_outputs(args,label_dict,master_bed):
    #faster to write in a single operation as a pandas df
    print("writing output file")
    tasks=list(label_dict.keys())
    num_tasks=len(tasks)
    num_entries=len(master_bed)
    data_for_np=np.zeros((num_entries,num_tasks))
    for i in range(num_tasks):
        cur_task=tasks[i]
        data_for_np[:,i]=label_dict[cur_task].get()
    data=pd.DataFrame(data=data_for_np,
                      index=[str(i).strip('\n') for i in master_bed],
                      columns=tasks)
    data.to_csv(args.outf,sep='\t',index_label="Chrom\tStart\tEnd")

def get_labels_one_task(task_bedfile,master_bed,training_bin_size):
    labels=[]
    processed=dict()
    print('thread!')
    try:
        intersection=master_bed.intersect(task_bedfile,wao=True)
        print("Got intersection!") 
    except:
        return None 
    for entry in intersection:
        overlap=int(entry[-1])
        master_peak='\t'.join(entry[0:3])
        if master_peak in processed:
            continue
        processed[master_peak]=1
        #if there is no overlap with the task bed file, the label is 0 
        if overlap==0:
            labels.append(0)
        else:
            summit_offset=int(entry[-2])
            peak_start=int(entry[4])
            peak_end=int(entry[5])
            if summit_offset!=-1:
                #check to see if the summit falls in the master_peak
                summit=int(entry[4])+summit_offset
                master_peak_start=int(entry[1])
                master_peak_end=int(entry[2])
                if (master_peak_start <= summit <= master_peak_end):
                    labels.append(1)
                else:
                    labels.append(-1)
            else:
                task_peak_size=peak_end-peak_start
                tocheck=min([task_peak_size,training_bin_size])
                if overlap/tocheck > 0.5:
                    #it's a positive, more than 1/2 the task peak overlaps with master bed region (or vice versa)
                    labels.append(1)
                else:
                    #ambiguous, exclude from analysis
                    labels.append(-1)
    return labels
    
def get_labels(task_bedfiles,master_bed,training_bin_size,threads):
    pool=ThreadPool(threads)    
    print("getting labels with " +str(threads)+ " threads" )
    label_dict=dict()
    for task in task_bedfiles:
        print(task)
        task_bedfile=task_bedfiles[task]
        label_dict[task]=pool.apply_async(get_labels_one_task,args=(task_bedfile,master_bed,training_bin_size))
    pool.close()
    pool.join()                                           
    print("got labels") 
    return label_dict


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
    label_dict=get_labels(task_bedfiles,master_bed,args.training_bin_size,args.threads)

    #write the labels to output file 
    write_outputs(args, label_dict,master_bed)
    
        
    

if __name__=="__main__":
    main()
    
    
