## helper functions to ingest bigwig and narrowPeak data files into a tileDB instance.
## tileDB instances are indexed by coordinate
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import tiledb
import argparse
import pandas as pd
import numpy as np
from collections import OrderedDict
from .attrib_config import *
from .utils import *
from multiprocessing.pool import Pool, ThreadPool 
import pdb
import gc

#graceful shutdown
import psutil
import signal 
import os


#config
tdb_Config=tiledb.Config({"sm.check_coord_dups":"false",
                          "sm.check_coord_oob":"false",
                          "sm.check_global_order":"false",
                          "sm.num_writer_threads":"50",
                          "sm.num_reader_threads":"50",
                          "sm.num_async_threads":"50",
                          "vfs.num_threads":"50"})    
def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def kill_child_processes(parent_pid, sig=signal.SIGTERM):
    try:
        parent = psutil.Process(parent_pid)
    except psutil.NoSuchProcess:
        return
    children = parent.children(recursive=True)
    for process in children:
        process.send_signal(sig)
        
def args_object_from_args_dict(args_dict):
    #create an argparse.Namespace from the dictionary of inputs
    args_object=argparse.Namespace()
    #set the defaults
    vars(args_object)['chrom_threads']=1
    vars(args_object)['overwrite']=False
    vars(args_object)['batch_size']=10000000
    vars(args_object)['tile_size']=90000
    vars(args_object)['attribute_config']='encode_pipeline'
    for key in args_dict:
        vars(args_object)[key]=args_dict[key]
    #set any defaults that are unset 
    args=args_object    
    return args 
    
def parse_args():
    parser=argparse.ArgumentParser(description="ingest data into tileDB")
    parser.add_argument("--tiledb_metadata",help="fields are: dataset, fc_bigwig, pval_bigwig, count_bigwig_plus_5p, count_bigwig_minus_5p, idr_peak, overlap_peak, ambig_peak")
    parser.add_argument("--tiledb_group")
    parser.add_argument("--overwrite",default=False,action="store_true") 
    parser.add_argument("--chrom_sizes",help="2 column tsv-separated file. Column 1 = chromsome name; Column 2 = chromosome size")
    parser.add_argument("--chrom_threads",type=int,default=1,help="inner thread pool, launched by task_threads")
    parser.add_argument("--batch_size",type=int,default=10000000,help="num entries to write at once")
    parser.add_argument("--tile_size",type=int,default=90000,help="tile size")
    parser.add_argument("--attribute_config",default='encode_pipeline',help="the following are supported: encode_pipeline, generic_bigwig")
    return parser.parse_args()

def create_new_array(size,
                     array_out_name,
                     tile_size,
                     attribute_config,
                     compressor='gzip',
                     compression_level=-1):
    '''
    Creates an empty tileDB array
    '''
    
    tile_size=min(size,tile_size)    
    tiledb_dim = tiledb.Dim(
        name='genome_coordinate',
        domain=(0, size - 1),
        tile=tile_size,
        dtype='uint32')
    tiledb_dom = tiledb.Domain(tiledb_dim,ctx=tiledb.Ctx(config=tdb_Config) )

    #generate the attribute information
    attribute_info=get_attribute_info(attribute_config)
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
    

def extract_metadata_field(row,field):
    dataset=row['dataset'] 
    try:
        return row[field]
    except:
        print("tiledb_metadata has no column "+field+" for dataset:"+str(dataset))
        return None

def open_data_for_parsing(row,attribute_info):
    try:
        data_dict={}
        cols=list(row.index)
        if 'dataset' in cols:
            cols.remove('dataset')
        for col in cols:
            cur_fname=extract_metadata_field(row,col)
            if cur_fname is not None:
                data_dict[col]=attribute_info[col]['opener'](cur_fname)
        return data_dict
    except Exception as e:
        print(repr(e))
        raise e
    
def parse_input_chunks(chrom,size,parser_chunk,data_dict,attribute_info,attribute,cur_parser,args):
    try:
        vals=np.empty((size,))
        pool_inputs=[] 
        for chrom_pos in range(0,size,parser_chunk):
            cur_start=chrom_pos
            cur_end=min([cur_start+parser_chunk,size])
            pool_inputs.append((data_dict[attribute],chrom,cur_start,cur_end,attribute_info[attribute]))
        #iterate through the tasks 
        #with Pool(args.chrom_threads,initializer=init_worker) as pool:
        with ThreadPool(args.chrom_threads) as pool:
            results=pool.map(cur_parser,pool_inputs)
        pool.close()
        pool.join()
        for elt in results:
            elt_start=elt[0]
            elt_end=elt[1]
            elt_val=elt[2]
            vals[elt_start:elt_end]=elt_val
        return vals
    except KeyboardInterrupt:
        print("keyboard interrupt detected")
        #shutdown the pool
        pool.terminate()
        # Kill remaining child processes
        kill_child_processes(os.getpid())
        raise 
    except Exception as e:
        print(repr(e))
        #shutdown the pool
        pool.terminate() 
        # Kill remaining child processes
        kill_child_processes(os.getpid())
        raise e                
    
def ingest(args):
    try:
        if type(args)==type({}):
            args=args_object_from_args_dict(args)
        overwrite=args.overwrite
        chrom_threads=args.chrom_threads
        batch_size=args.batch_size
        tile_size=args.tile_size
        attribute_config=args.attribute_config
        updating=False

        attribute_info=get_attribute_info(args.attribute_config) 
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
        for task_index,task_row in tiledb_metadata.iterrows():            
            dataset=task_row['dataset']    
            #read in filenames for bigwigs
            data_dict=open_data_for_parsing(task_row,attribute_info)
            array_outf_prefix="/".join([args.tiledb_group,dataset])
            pool_inputs=[]
            for chrom_index, chrom_row in chrom_sizes.iterrows():
                chrom=chrom_row[0]
                size=chrom_row[1]
                array_out_name='.'.join([array_outf_prefix,chrom])                
                if tiledb.object_type(array_out_name) == "array":
                    if overwrite==False:
                        raise Exception("array:"+str(array_out_name) + "already exists; use the --overwrite flag to overwrite it. Exiting")
                    else:
                        print("warning: the array: "+str(array_out_name)+" already exists. You provided the --overwrite flag, so it will be updated/overwritten")
                        updating=True
                else:
                    #create the array:
                    create_new_array(size=size,
                                     attribute_config=attribute_config,
                                     array_out_name=array_out_name,
                                     tile_size=tile_size)
                    print("created new array:"+str(array_out_name))
                pool_inputs.append((data_dict,attribute_info,chrom,size,array_out_name,updating,args))
            with Pool(chrom_threads,initializer=init_worker) as pool:
            #with ThreadPool(chrom_threads) as pool:
                print("made pool")
                res=pool.map(process_chrom,pool_inputs)
            pool.close()
            pool.join()
            print("wrote chrom array for task:"+str(dataset))
    except KeyboardInterrupt:
        print('detected keyboard interrupt')
        #shutdown the pool
        pool.terminate() 
        # Kill remaining child processes
        kill_child_processes(os.getpid())
        raise 
    except Exception as e:
        print(repr(e))
        #shutdown the pool
        pool.terminate() 
        # Kill remaining child processes
        kill_child_processes(os.getpid())
        raise e

def process_chrom(inputs):
    try:
        data_dict=inputs[0]
        attribute_info=inputs[1]
        chrom=inputs[2]
        size=inputs[3]
        array_out_name=inputs[4]
        updating=inputs[5]
        args=inputs[6]
        attribute_config=args.attribute_config
        batch_size=args.batch_size
        dict_to_write=OrderedDict() 
        for attribute in data_dict:
            cur_parser=attribute_info[attribute]['parser']
            print(cur_parser)
            if cur_parser is parse_bigwig_chrom_vals:
               dict_to_write[attribute]=parse_input_chunks(chrom,size,batch_size,data_dict,attribute_info,attribute,cur_parser,args)
            else:
                dict_to_write[attribute]=cur_parser(data_dict[attribute],chrom,0,size,attribute_info[attribute])
            print("got:"+str(attribute)+" for chrom:"+str(chrom))
            compressor='gzip'
            compression_level=-1

        if updating is True:
            #we are only updating some attributes in the array
            print(array_out_name)
            with tiledb.DenseArray(array_out_name,mode='r') as cur_array:
                cur_vals=cur_array[:]            
            print('got cur vals') 
            for key in dict_to_write:
                cur_vals[key]=dict_to_write[key]
            dict_to_write=cur_vals
            print("updated data dict for writing") 
        else:
            #we are writing for the first time, make sure all attributes are provided, if some are not, use a nan array
            required_attrib=list(get_attribute_info(attribute_config).keys())
            for attrib in required_attrib:
                if attrib not in dict_to_write:
                    dict_to_write[attrib]=np.full(size,np.nan)
        with tiledb.DenseArray(array_out_name,ctx=tiledb.Ctx(config=tdb_Config) ,mode='w') as out_array:
            out_array[:]=dict_to_write
            print("wrote to disk:"+array_out_name)
    except KeyboardInterrupt:
        print('detected keyboard interrupt')
        # Kill remaining child processes
        kill_child_processes(os.getpid())
        raise 
    except Exception as e:
        print(repr(e))
        # Kill remaining child processes
        kill_child_processes(os.getpid())
        raise e
     

def main():
    args=parse_args()
    ingest(args)
        
if __name__=="__main__":
    main() 
    
    
