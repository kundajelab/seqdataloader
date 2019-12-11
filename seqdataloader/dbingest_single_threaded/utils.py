import pandas as pd
import numpy as np
import pyBigWig
import pdb 
from itertools import islice

def open_bigwig_for_parsing(fname):
    return fname

def open_csv_for_parsing(fname):
    return fname


def parse_bigwig_chrom_vals(bigwig_name,chrom,start,end,cur_attribute_info):
    store_summits=None
    summit_indicator=None
    if 'store_summits' in cur_attribute_info:
        store_summits=cur_attribute_info['store_summits']
    if 'summit_indicator' in cur_attribute_info:
        summit_indicator=cur_attribute_info['summit_indicator']
    bigwig_object=pyBigWig.open(bigwig_name)
    #note: pybigwig uses NA in place of 0 where there are no reads, replace with 0.
    #check to see if chromosome in bigwig, if not, return all NA's & warning that chromosome is not present in the dataset
    bw_chroms=bigwig_object.chroms().keys()
    if chrom not in bw_chroms:
        print("WARNING: chromosome:"+str(chrom) +" was not found in the bigwig file:"+str(bigwig_name))
        size=(end-start+1)
        signal_data=np.full(size,np.nan)
    else: 
        signal_data=np.nan_to_num(bigwig_object.values(chrom,start,end))
    return signal_data


def parse_narrowPeak_chrom_vals(narrowPeak_fname,chrom,start,end,cur_attribute_info):
    store_summits=None
    summit_indicator=None
    if 'store_summits' in cur_attribute_info:
        store_summits=cur_attribute_info['store_summits']
    if 'summit_indicator' in cur_attribute_info:
        summit_indicator=cur_attribute_info['summit_indicator']

    narrowPeak_df=pd.read_csv(narrowPeak_fname,header=None,sep='\t')
    signal_data = np.zeros(end, dtype=np.int)
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
    return signal_data 
