import pandas as pd
import numpy as np
import pyBigWig

def open_bigwig_for_parsing(fname):
    #return pyBigWig.open(fname)
    return fname

def open_csv_for_parsing(fname):
    return pd.read_csv(fname,header=None,sep='\t') 

def parse_bigwig_chrom_vals(bigwig_name,chrom,size,store_summits,summit_indicator):
    bigwig_object=pyBigWig.open(bigwig_name)
    signal_data = np.zeros(size, dtype=np.float32)
    signal_data[:]=bigwig_object.values(chrom,0,size) 
    return signal_data

def parse_narrowPeak_chrom_vals(narrowPeak_df,chrom,size,store_summits,summit_indicator):
    signal_data = np.zeros(size, dtype=np.int)
    chrom_subset_df=narrowPeak_df[narrowPeak_df[0]==chrom]
    nitems_in_row=chrom_subset_df.shape[1] 
    for index,row in chrom_subset_df.iterrows():
        signal_data[row[1]:row[2]]=1
        if store_summits is True:
            #we assume standard narrowPeak format where summit is stored in final column
            summit_pos=row[1]+row[nitems_in_row-1]
            signal_data[summit_pos]=summit_indicator            
    return signal_data 

    
