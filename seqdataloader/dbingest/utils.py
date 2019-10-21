import pandas as pd
import numpy as np
import pyBigWig
from itertools import islice

def open_bigwig_for_parsing(fname):
    return fname

def open_csv_for_parsing(fname):
    return fname


def parse_bigwig_chrom_vals(bigwig_name,chrom,size,store_summits,summit_indicator):
    bigwig_object=pyBigWig.open(bigwig_name)
    #note: pybigwig uses NA in place of 0 where there are no reads, replace with 0.
    signal_data=np.nan_to_num(bigwig_object.values(chrom,0,size))
    return signal_data

def parse_narrowPeak_chrom_vals(narrowPeak_fname,chrom,size,store_summits,summit_indicator):
    narrowPeak_df=pd.read_csv(narrowPeak_fname,header=None,sep='\t') 
    signal_data = np.zeros(size, dtype=np.int)
    chrom_subset_df=narrowPeak_df[narrowPeak_df[0]==chrom]
    del narrowPeak_df
    nitems_in_row=chrom_subset_df.shape[1]
    if (store_summits==True) and (chrom_subset_df.shape[1]<10):
        print("warning! dataset does not have 10 columns of the standard narrowPeak format, falling back to peak centers")
    summits=[]
    for index,row in chrom_subset_df.iterrows():
        signal_data[row[1]:row[2]]=1
        #add in summits in a separate step to avoid overwriting them with "1's" for overlaping peak coordinates;
        #The overwriting issue is particularly relevant for pseudobulk data. 
        if store_summits is True:
            if len(row)<10:
                summit_pos=int(round(row[1]+0.5*(row[2]-row[1])))
            else: 
                summit_pos=row[1]+row[nitems_in_row-1]
            summits.append(summit_pos)
    if store_summits is True:
        signal_data[summits]=summit_indicator
    assert min(signal_data)==0 #sanity check 
    return signal_data 

    
def chunkify(iterable,chunk):
    it=iter(iterable)
    while True:
        piece=list(islice(it,chunk))
        if piece:
            yield piece
        else:
            return
        
