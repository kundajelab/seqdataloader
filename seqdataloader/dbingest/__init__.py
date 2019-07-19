## helper functions to ingest bigwig and narrowPeak data files into a tileDB instance.
## tileDB instances are indexed by coordinate, and a task(i.e. a unique identifier) 
## attributes include cell type, assay, project, filetype (i.e. fc bigwig, pval bigwig, 
import tiledb
import argparse
import pandas as pd
import numpy as np
from attrib_config import *
from utils import *
import multiprocessing as mp 

def parse_args():
    parser=argparse.ArgumentParser(description="ingest data into tileDB")
    parser.add_argument("--tiledb_metadata",help="fileds are: Dataset, FC_bigwig, PVAL_bigwig, COUNT_bigwig, IDR_peaks, OVERLAP_peaks, AMBIG_peaks")
    parser.add_argument("--tiledb_group")
    parser.add_argument("--overwrite",default=False,action="store_true") 
    parser.add_argument("--chrom_sizes",help="2 column tsv-separated file. Column 1 = chromsome name; Column 2 = chromosome size")
    parser.add_argument("--threads") 
    return parser.parse_args()

def create_new_array(size,
                     array_out_name,
                     default_tile_size=9000,
                     compressor='gzip',
                     compression_level=-1):
    '''
    Creates an empty tileDB array
    '''
    tile_size=min(size,default_tile_size)    
    tiledb_dim = tiledb.Dim(
        name='genome_coordinate',
        domain=(0, size - 1),
        tile=tile_size,
        dtype='uint32')
    tiledb_dom = tiledb.Domain(tiledb_dim)
    
    #generate the attribute information
    attribute_info=get_attribute_info()
    attribs=[]
    for key in attribute_info:
        attribs.append(tiledb.Attr(
            name=key,
            filters=tiledb.FilterList([tiledb.GzipFilter()]),
            dtype=attrib_info[key]['dtype']))
    tiledb_schema = tiledb.ArraySchema(
        domain=tiledb_dom,
        attrs=tuple(attribs),
        cell_order='row-major',
        tile_order='row-major')

    tiledb.DenseArray.create(array_out_name, tiledb_schema)
    print("created empty array on disk") 
    return

def write_array_to_tiledb(size,
                          dict_to_write,
                          array_out_name,
                          default_tile_size=9000,
                          compressor='gzip',
                          compression_level=-1,
                          updating=False):
    if updating is True:
        #we are only updating some attributes in the array
        with tiledb.DenseArray(array_out_name,mode='r') as cur_array:
            cur_vals=cur_array[:]
        for key in dict_to_write:
            cur_vals[key]=dict_to_write[key]
        dict_to_write=cur_vals
        
    else:
        #we are writing for the first time, make sure all attributes are provided, if some are not, use a nan array
        required_attrib=list(get_attrib_info().keys())
        for attrib in required_attrib:
            if attrib not in dict_to_write:
                signal_data=np.zeros(size)
                signal_data[:]=np.nan
                dict_to_write[attrib]=signal_data
                
    with tiledb.DenseArray(array_out_name, mode='w') as out_array:
        out_array[:] = dict_to_write

def extract_metadata_field(row,field):
    dataset=row['dataset'] 
    try:
        return row[field]
    except:
        print("tiledb_metadata has no column "+field+" for dataset:"+str(dataset))
        return None

def open_data_for_parsing(row,attribute_info):
    data_dict={}
    cols=list(row.index)
    if 'dataset' in cols:
        cols.remove('dataset')
    for col in cols:
        cur_fname=extract_metadata_field(row,col)
        if cur_fname is not None:
            data_dict[col]=attribute_info[col]['opener'](cur_fname)            
    return data_dict

def create_tiledb_array(row,args,chrom_sizes,attribute_info):
    '''
    create new tileDB array for a single dataset, overwrite if array exists and user sets --overwrite flag
    '''

    #get the current dataset 
    dataset=row['dataset']    
    #read in filenames for bigwigs
    data_dict=open_data_for_parsing(row,attribute_info)
    
    array_outf_prefix="/".join([args.tiledb_group,dataset])
    for index, row in chrom_sizes.iterrows():
        chrom=row[0]
        size=row[1]
        array_out_name='.'.join([array_outf_prefix,chrom])
        updating=False
        if tiledb.object_type(array_out_name) == "array":
            if args.overwrite==False:
                raise Exception("array:"+str(array_out_name) + "already exists; use the --overwrite flag to overwrite it. Exiting")
            else:
                print("warning: the array: "+str(array_out_name)+" already exists. You provided the --overwrite flag, so it will be updated/overwritten")
                updating=True
        else:
            #create the array:
            create_new_array(size=size,
                             array_out_name=array_out_name)
            print("created new array:"+str(array_out_name))
                  
        dict_to_write=dict()
        for attribute in data_dict:
            dict_to_write[attribute]=attribute_info[attribute]['parser'](data_dict[attribute],chrom,size)
            print("got:"+str(attribute)+" for chrom:"+str(chrom))
        write_array_to_tiledb(size=size,
                              dict_to_write=dict_to_write,
                              array_out_name=array_out_name,
                              updating=updating)
        print("wrote array to disk for dataset:"+str(dataset))         

def main():
    args=parse_args()
    attribute_info=get_attribute_info() 
    tiledb_metadata=pd.read_csv(args.tiledb_metadata,header=0,sep='\t')
    print("loaded tiledb metadata")
    chrom_sizes=pd.read_csv(args.chrom_sizes,header=None,sep='\t')
    print("loaded chrom sizes")
    
    #check if the tiledb_group exists, and if not, create it
    if tiledb.object_type(args.tiledb_group) is not 'group':        
        group_uri=tiledb.group_create(args.tiledb_group)
        print("created tiledb group") 
    else:
        group_uri=args.tiledb_group
        print("tiledb group already exists")
        
    for index,row in tiledb_metadata.iterrows():
        create_tiledb_array(row,args,chrom_sizes,attribute_info)
    
if __name__=="__main__":
    main() 
    
    
