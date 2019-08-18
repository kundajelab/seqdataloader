## helper functions to ingest bigwig and narrowPeak data files into a tileDB instance.
## tileDB instances are indexed by coordinate
import tiledb
import argparse
import pandas as pd
import numpy as np
from .attrib_config import *
from .utils import *
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool


def parse_args():
    parser=argparse.ArgumentParser(description="ingest data into tileDB")
    parser.add_argument("--tiledb_metadata",help="fileds are: Dataset, fc_bigwig, pval_bigwig, count_bigwig_plus_5p, count_bigwig_minus_5p, idr_peaks, overlap_peaks, ambig_peaks")
    parser.add_argument("--tiledb_group")
    parser.add_argument("--overwrite",default=False,action="store_true") 
    parser.add_argument("--chrom_sizes",help="2 column tsv-separated file. Column 1 = chromsome name; Column 2 = chromosome size")
    parser.add_argument("--chrom_threads",type=int,default=1,help="inner thread pool, launched by task_threads")
    parser.add_argument("--task_threads",type=int,default=1,help="outer thread pool,launched by main thread")
    parser.add_argument("--store_summits",action="store_true")
    parser.add_argument("--summit_indicator",type=int,default=2,help="integer value to use as indicator of whether a base pair is a peak summit")
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
            dtype=attribute_info[key]['dtype']))
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
    print("starting to write output") 
    if updating is True:
        #we are only updating some attributes in the array
        with tiledb.DenseArray(array_out_name,mode='r') as cur_array:
            cur_vals=cur_array[:]
        for key in dict_to_write:
            cur_vals[key]=dict_to_write[key]
        dict_to_write=cur_vals
        print("updated data dict for writing") 
    else:
        #we are writing for the first time, make sure all attributes are provided, if some are not, use a nan array
        required_attrib=list(get_attribute_info().keys())
        for attrib in required_attrib:
            if attrib not in dict_to_write:
                dict_to_write[attrib]=np.full(size,np.nan)
    print("finalizing the write")
    with tiledb.DenseArray(array_out_name, mode='w') as out_array:
            out_array[:]=dict_to_write
    print("done writing")
    return

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

def process_chrom(inputs):
    chrom=inputs[0]
    size=inputs[1]
    array_out_name=inputs[2]
    data_dict=inputs[3]
    attribute_info=inputs[4]
    overwrite=inputs[5]
    store_summits=inputs[6]
    summit_indicator=inputs[7] 
    updating=False
    if tiledb.object_type(array_out_name) == "array":
        if overwrite==False:
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
        dict_to_write[attribute]=attribute_info[attribute]['parser'](data_dict[attribute],chrom,size,store_summits,summit_indicator)
        print("got:"+str(attribute)+" for chrom:"+str(chrom))

    write_array_to_tiledb(size=size,
                          dict_to_write=dict_to_write,
                          array_out_name=array_out_name,
                          updating=updating)
    print("wrote array to disk for dataset:"+str(array_out_name))         
    return 'done'

def create_tiledb_array(inputs):
    '''
    create new tileDB array for a single dataset, overwrite if array exists and user sets --overwrite flag
    '''
    row=inputs[0]
    args=inputs[1]
    chrom_sizes=inputs[2]
    attribute_info=inputs[3] 
    #get the current dataset 
    dataset=row['dataset']    
    #read in filenames for bigwigs
    data_dict=open_data_for_parsing(row,attribute_info)
    pool=ThreadPool(args.chrom_threads)
    pool_inputs=[] 
    array_outf_prefix="/".join([args.tiledb_group,dataset])
    for index, row in chrom_sizes.iterrows():
        chrom=row[0]
        size=row[1]
        array_out_name='.'.join([array_outf_prefix,chrom])
        pool_inputs.append((chrom,size,array_out_name,data_dict,attribute_info,args.overwrite,args.store_summits, args.summit_indicator))
    pool_outputs=pool.map(process_chrom,pool_inputs)
    pool.close()
    pool.join()
    return "done"

def args_object_from_args_dict(args_dict):
    #create an argparse.Namespace from the dictionary of inputs
    args_object=argparse.Namespace()
    #set the defaults
    vars(args_object)['chrom_threahds']=1
    vars(args_object)['task_threads']=1
    vars(args_object)['overwrite']=False
    vars(args_object)['store_summits']=True
    vars(args_object)['summit_indicator']=2
    for key in args_dict:
        vars(args_object)[key]=args_dict[key]
    #set any defaults that are unset 
    args=args_object    
    return args 
    
def ingest(args):
    if type(args)==type({}):
        args=args_object_from_args_dict(args)
        
    attribute_info=get_attribute_info() 
    tiledb_metadata=pd.read_csv(args.tiledb_metadata,header=0,sep='\t')
    print("loaded tiledb metadata")
    chrom_sizes=pd.read_csv(args.chrom_sizes,header=None,sep='\t')
    print("loaded chrom sizes")
    pool=ThreadPool(args.task_threads)
    pool_inputs=[] 
    #check if the tiledb_group exists, and if not, create it
    if tiledb.object_type(args.tiledb_group) is not 'group':        
        group_uri=tiledb.group_create(args.tiledb_group)
        print("created tiledb group") 
    else:
        group_uri=args.tiledb_group
        print("tiledb group already exists")
        
    for index,row in tiledb_metadata.iterrows():
        pool_inputs.append([row,args,chrom_sizes,attribute_info])        
    task_returns=pool.map(create_tiledb_array,pool_inputs)
    pool.close()
    pool.join()

def main():
    args=parse_args()
    ingest(args)
    
    
if __name__=="__main__":
    main() 
    
    
