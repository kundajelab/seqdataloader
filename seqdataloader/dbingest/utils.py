import pandas as pd
import numpy as np
import pyBigWig

def open_bigwig_for_parsing(fname):
    #return pyBigWig.open(fname)
    return fname

def open_csv_for_parsing(fname):
    return pd.read_csv(fname,header=None,sep='\t') 

def parse_bigwig_chrom_vals(bigwig_name,chrom,size):
    bigwig_object=pyBigWig.open(bigwig_name)
    signal_data = np.zeros(size, dtype=np.float32)
    signal_data[:]=bigwig_object.values(chrom,0,size) 
    return signal_data

def parse_narrowPeak_chrom_vals(narrowPeak_df,chrom,size):
    signal_data = np.zeros(size, dtype=np.int)
    chrom_subset_df=narrowPeak_df[narrowPeak_df[0]==chrom]
    for index,row in chrom_subset_df.iterrows():
        signal_data[row[1]:row[2]]=1
    return signal_data 

    
