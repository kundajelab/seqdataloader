import pandas as pd
import numpy as np
import pyBigWig
from pybedtools import BedTool
from itertools import islice
from collections import OrderedDict

def open_bigwig_for_parsing(fname,parallel=False):
    if not parallel:
        return pyBigWig.open(fname)
    else:
        #pybigwig objects cannot be pickled
        return fname 
def open_csv_for_parsing(fname,parallel=False):
    #if not parallel:
    return BedTool(fname)
    #else:
        #
    #    return fname

def parse_bigwig_chrom_vals(entry):
    bigwig_object=entry[0]
    if type(bigwig_object)==str:
        bigwig_object=pyBigWig.open(bigwig_object)
    chrom=entry[1]
    start=entry[2]
    end=entry[3]
    cur_attribute_info=entry[4]
    #note: pybigwig uses NA in place of 0 where there are no reads, replace with 0.
    bw_chroms=bigwig_object.chroms().keys()
    if chrom not in bw_chroms:
        print("WARNING: chromosome:"+str(chrom)+ " was not found in the bigwig file:"+str(bigwig_object))
        size=(end-start+1)
        signal_data=np.full(size,np.nan)
    else: 
        #check to see if chromosome in bigwig, if not, return all NA's & warning that chromosome is not present in the dataset
        try:
            signal_data=np.nan_to_num(bigwig_object.values(chrom,start,end))
        except Exception as e:
            print(chrom+"\t"+str(start)+"\t"+str(end)+str(cur_attribute_info))
            raise e
    return start, end, signal_data


def parse_narrowPeak_chrom_vals(entry):
    task_bed=entry[0]
    chrom=entry[1]
    start=entry[2]
    end=entry[3]
    num_entries=end-start
    chrom_coords=chrom+'\t'+str(start)+'\t'+str(end)
    chrom_bed=BedTool(chrom_coords,from_string=True)
    cur_bed=task_bed.intersect(chrom_bed)
    cur_attribute_info=entry[4]
    store_summits=None
    summit_indicator=None
    summit_from_peak_center=None
    if 'store_summits' in cur_attribute_info:
        store_summits=cur_attribute_info['store_summits']
        if store_summits is True:
            summit_from_peak_center=cur_attribute_info['summit_from_peak_center'] 
            summit_indicator=cur_attribute_info['summit_indicator']
    signal_data = np.zeros(num_entries, dtype=np.int)
    warned=False
    summits=[]
    for entry in cur_bed:
        #offset relative to start position of the interval
        entry_start=int(entry[1])-start
        entry_end=int(entry[2])-start
        signal_data[entry_start:entry_end]=1
        #add in summits in a separate step to avoid overwriting them with "1's" for overlaping peak coordinates;
        #The overwriting issue is particularly relevant for pseudobulk data. 
        if store_summits is True:
            if summit_from_peak_center is True:
                summit_pos=int(entry_start+(entry_end-entry_start)*0.5)
            else:
                try:
                    summit_pos=entry_start+int(entry[-1])
                except:
                    print("WARNING: could not add summit position from last column of narrowPeak file, falling back to peak center"+str(entry))
                    summit_pos=int(entry_start+(entry_end-entry_start)*0.5)                
            if (summit_pos < entry_end) and (summit_pos > entry_start):
                summits.append(summit_pos)
            else:
                print("WARNING: summit position outside peak region position,skipping:"+str(entry))                    
    if store_summits is True:
        signal_data[summits]=summit_indicator
    return start, end, signal_data 
    
def chunkify(iterable,chunk):
    it=iter(iterable)
    while True:
        piece=list(islice(it,chunk))
        if piece:
            yield piece
        else:
            return
        
def transform_indices_to_chrom_coords(start_chunk_index,end_chunk_index,chrom_indices):
    #print("transforming:"+str(start_chunk_index)+'-'+str(end_chunk_index))
    #print(str(chrom_indices))
    chrom_coords=[]
    while True:
        found=False
        for chrom in chrom_indices:
            cur_chrom_start_index=chrom_indices[chrom][0]
            cur_chrom_end_index=chrom_indices[chrom][1]
            cur_chrom_size=chrom_indices[chrom][2]
            if start_chunk_index >=cur_chrom_start_index:
                if start_chunk_index < cur_chrom_end_index:
                    found=True
                    #start position is on this chromosome
                    chrom_coord_start=start_chunk_index-cur_chrom_start_index
                    #check if end coordinate falls on same chromosome
                    if end_chunk_index < cur_chrom_end_index:
                        #on one chrom
                        chrom_coord_end=end_chunk_index-cur_chrom_start_index
                        chrom_coords.append((chrom,chrom_coord_start,chrom_coord_end,start_chunk_index,end_chunk_index))
                        return chrom_coords
                    else:
                        chrom_coord_end=cur_chrom_size
                        chrom_coords.append((chrom,chrom_coord_start,chrom_coord_end,start_chunk_index,cur_chrom_end_index))
                        #update start_chunk_index
                        start_chunk_index=cur_chrom_end_index
        if found is False:
            if len(chrom_coords)==0:
                raise Exception("failed to tranform indices:"+str(start_chunk_index)+"-"+str(end_chunk_index)+ " to chrom coords;"+str(chrom_indices))
            else:
                return chrom_coords
        
def transform_chrom_size_to_indices(chrom_sizes):
    '''
    chrom_sizes is a dataframe 
    get 0-based tdb coordinates for start & end of each chromosome 
    '''
    start_coord=0
    chrom_indices=OrderedDict() 
    for index,row in chrom_sizes.iterrows():
        chrom=row[0]
        size=row[1]
        chrom_indices[chrom]=[start_coord,start_coord+size,size]
        start_coord=start_coord+size
    return chrom_indices,start_coord


