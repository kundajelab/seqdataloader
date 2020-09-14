## helper functions to ingest bigwig and narrowPeak data files into a tileDB instance.
## tileDB instances are indexed by coordinate
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import math
import psutil
from multiprocessing import Pool, Process, Queue
import os
import signal
import tiledb
import pickle
import argparse
import pandas as pd
import numpy as np
from collections import OrderedDict
from ..attrib_config import *
from ..queue_config import * 
from ..utils import *
from ..tdb_config import * 
import gc
import time
import sys

def args_object_from_args_dict(args_dict):
    #create an argparse.Namespace from the dictionary of inputs
    args_object=argparse.Namespace()
    #set the defaults
    vars(args_object)['overwrite']=False
    vars(args_object)['coord_tile_size']=10000
    vars(args_boject)['task_tile_size']=1
    vars(args_object)['attribute_config']=None
    vars(args_object)['attribute_config_file']=None
    vars(args_object)['write_chunk']=30000000
    vars(args_object)['threads']=1
    vars(args_object)['max_queue_size']=30
    vars(args_object)['max_mem_g']=100
    for key in args_dict:
        vars(args_object)[key]=args_dict[key]
    #set any defaults that are unset 
    args=args_object    
    return args 
    
def parse_args():
    parser=argparse.ArgumentParser(description="ingest data into tileDB")
    parser.add_argument("--tiledb_metadata",help="each row is a dataset, each column corresponds to an attribute")
    parser.add_argument("--array_name")
    parser.add_argument("--overwrite",default=False,action="store_true") 
    parser.add_argument("--chrom_sizes",help="2 column tsv-separated file. Column 1 = chromsome name; Column 2 = chromosome size")
    parser.add_argument("--coord_tile_size",type=int,default=10000,help="coordinate axis tile size")
    parser.add_argument("--task_tile_size",type=int,default=1,help="task axis tile size")
    parser.add_argument("--attribute_config",default=None,help="the following are supported: encode_pipeline, encode_pipeline_with_controls, generic_bigwig")
    parser.add_argument("--attribute_config_file",default=None,help="file with 2 columns; first column indicates attribute name; 2nd column indicates attribute type, which is one of bigwig, bed_no_summit, bed_summit_from_peak_center, bed_summit_from_last_col")
    parser.add_argument("--write_chunk",type=int,default=30000000,help="number of bases to write to disk in one tileDB DenseArray write operation")
    parser.add_argument("--threads",type=int,default=1,help="number of chunks to process in parallel")
    parser.add_argument("--max_queue_size",type=int,default=30)
    parser.add_argument("--max_mem_g",type=int,default=100,help="maximum memory usage in Gigabytes")
    return parser.parse_args()
    
    
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

                
def create_new_array(tdb_Context,
                     size,
                     array_out_name,
                     coord_tile_size,
                     task_tile_size,
                     attribute_config,
                     attribute_config_file,
                     compressor='gzip',
                     compression_level=-1,
                     var=False):
    '''
    Creates an empty tileDB array
    size= tuple(num_indices,num_tasks)
    '''
    coord_tile_size=min(size[0],coord_tile_size)
    task_tile_size=max([1,min(size[1],task_tile_size)])
    tiledb_dim_coords = tiledb.Dim(
        name='genome_coordinate',
        domain=(0, size[0]),
        tile=coord_tile_size,
        dtype='uint32')
    tiledb_dim_tasks=tiledb.Dim(
        name='task',
        domain=(0,size[1]),#max([1,size[1]])),
        tile=task_tile_size,
        dtype='uint32')
    tiledb_dom = tiledb.Domain(tiledb_dim_coords,tiledb_dim_tasks,ctx=tdb_Context)

    #generate the attribute information
    attribute_info=get_attribute_info(attribute_config,attribute_config_file)
    attribs=[]
    for key in attribute_info:
        attribs.append(tiledb.Attr(
            name=key,
            var=var,
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
            if isinstance(cur_fname,str):
                assert os.path.exists(cur_fname), "The path:"+str(cur_fname)+" does not exist. If you meant to skip this column, leave it empty in the metadata sheet." 
            elif math.isnan(float(cur_fname)):
                continue
            elif cur_fname is  None:
                continue
            data_dict[col]=attribute_info[col]['opener'](cur_fname,parallel=True)
        return data_dict
    except Exception as e:
        print(repr(e))
        kill_child_processes(os.getpid())
        raise
    
def ingest(args):
    if type(args)==type({}):
        args=args_object_from_args_dict(args)
    if args.write_chunk > max_write_chunk:
        print("WARNING: You have specified a write_chunk size of"+str(args.write_chunk)+" but the maximum supported with python serialization is:"+str(max_write_chunk)+". It will be reset to "+str(max_write_chunk))
        args.write_chunk=max_write_chunk

    #create a queue to write the array
    global write_queue
    write_queue=Queue(maxsize=args.max_queue_size)

    #config
    tdb_Config=tiledb.Config(tdb_config_params)
    tdb_write_Context=tiledb.Ctx(config=tdb_Config)   
    tdb_read_Context=tiledb.Ctx(config=tdb_Config)
    
    overwrite=args.overwrite
    coord_tile_size=args.coord_tile_size
    task_tile_size=args.task_tile_size
    attribute_config=args.attribute_config
    attribute_config_file=args.attribute_config_file
    updating=False

    attribute_info=get_attribute_info(args.attribute_config,args.attribute_config_file)
    tiledb_metadata=pd.read_csv(args.tiledb_metadata,header=0,sep='\t')
    num_tasks=tiledb_metadata.shape[0]
    print("num_tasks:"+str(num_tasks))
    
    print("loaded tiledb metadata")
    chrom_sizes=pd.read_csv(args.chrom_sizes,header=None,sep='\t')
    print("loaded chrom sizes")
    chrom_indices,num_indices=transform_chrom_size_to_indices(chrom_sizes)
    print("num_indices:"+str(num_indices))
    array_out_name=args.array_name
    if tiledb.object_type(array_out_name) == "array":
        if overwrite==False:
            raise Exception("array:"+str(array_out_name) + "already exists; use the --overwrite flag to overwrite it. Exiting")
        else:
            print("warning: the array: "+str(array_out_name)+" already exists. You provided the --overwrite flag, so it will be updated/overwritten")
            updating=True
    else:
        #create the array:
        create_new_array(tdb_Context=tdb_write_Context,
                         size=(num_indices,num_tasks-1),
                         attribute_config=attribute_config,
                         attribute_config_file=attribute_config_file,
                         array_out_name=array_out_name,
                         coord_tile_size=coord_tile_size,
                         task_tile_size=task_tile_size,
                         var=False)
        print("created new array:"+str(array_out_name))
        #create metadata array
        metadata_dict={}
        metadata_dict['tasks']=[i for i in tiledb_metadata['dataset']]
        metadata_dict['chroms']=[i for i in chrom_indices.keys()]
        metadata_dict['sizes']=[i[2] for i in list(chrom_indices.values())]
        metadata_dict['offsets']=[i[0] for i in list(chrom_indices.values())]
        num_tasks=tiledb_metadata['dataset'].shape[0]
        
        num_chroms=len(chrom_indices.keys())
        with tiledb.DenseArray(array_out_name,ctx=tdb_write_Context,mode='w') as cur_array:
            cur_array.meta['num_tasks']=num_tasks
            cur_array.meta['num_chroms']=num_chroms
            for task_index in range(num_tasks):
                cur_array.meta['_'.join(['task',str(task_index)])]=metadata_dict['tasks'][task_index]
            for chrom_index in range(num_chroms):
                cur_array.meta['_'.join(['chrom',str(chrom_index)])]=metadata_dict['chroms'][chrom_index]
                cur_array.meta['_'.join(['size',str(chrom_index)])]=metadata_dict['sizes'][chrom_index]
                cur_array.meta['_'.join(['offset',str(chrom_index)])]=metadata_dict['offsets'][chrom_index]                                
        print("created tiledb metadata")
    pool=Pool(processes=args.threads,initializer=init_worker)
    print("made pool") 
    pool_inputs=[] 
    for task_index,task_row in tiledb_metadata.iterrows():
        dataset=task_row['dataset']
        #read in filenames for bigwigs
        data_dict=open_data_for_parsing(task_row,attribute_info)
        for start_chunk_index in range(0,num_indices,args.write_chunk):
            end_chunk_index=start_chunk_index+min([num_indices,args.write_chunk])
            #convert global indices to chrom+pos indices
            chunk_chrom_coords=transform_indices_to_chrom_coords(start_chunk_index,end_chunk_index,chrom_indices)
            if chunk_chrom_coords is None:
                raise Exception("failed to tranform indices:"+str(start_chunk_index)+"-"+str(end_chunk_index)+ " to chrom coords;"+str(chrom_indices))
            for coord_set in chunk_chrom_coords:
                pool_inputs.append((task_index,data_dict,attribute_info,coord_set,args))
    pool_feed_chunk_start=0
    pool_feed_chunk_max=len(pool_inputs)
    chunks_to_process=len(pool_inputs)
    array_writer=Process(target=write_array,args=([args,updating,chunks_to_process]))
    try:
        array_writer.start()
    except Exception as e:
        raise e

    try:
        while pool_feed_chunk_start < pool_feed_chunk_max:
            pool_feed_chunk_end=min([pool_feed_chunk_start+queue_feed_chunk_size,pool_feed_chunk_max])
            #only do mapping if queue size is not exceeded & total memory consumption is not exceeded
            write_queue_size=write_queue.qsize()
            mem_used=psutil.virtual_memory().used / (10**9)
            print("mapping to pool, queue size:"+str(write_queue_size))
            print("mapping to pool, mem used:"+str(mem_used))
            while (write_queue_size >=args.max_queue_size) or (mem_used >=args.max_mem_g):
                time.sleep(10)
            print("sending to pool:"+str(pool_feed_chunk_start)+"-"+str(pool_feed_chunk_end)+"/"+str(chunks_to_process))
            pool.map(process_chunk,pool_inputs[pool_feed_chunk_start:pool_feed_chunk_end])
            pool_feed_chunk_start+=queue_feed_chunk_size
            time.sleep(60)
        pool.close()
    except KeyboardInterrupt:
        kill_child_processes(os.getpid())
        pool.terminate()
        raise
    except Exception as e:
        print(e)
        kill_child_processes(os.getpid())
        raise 
        
    #wait until we're done writing to the tiledb array
    array_writer.join()
    print("array_writer.join() is complete")
    print("shutting down pool")
    pool.join()
    print('done!') 

def process_chunk(inputs):
    try:
        task_index=inputs[0]
        data_dict=inputs[1]
        attribute_info=inputs[2]
        coord_set=inputs[3]
        args=inputs[4]

        attribute_config=args.attribute_config
        dict_to_write=OrderedDict()
        chrom=coord_set[0]
        start_pos=coord_set[1]
        end_pos=coord_set[2]
        start_index=coord_set[3]
        end_index=coord_set[4] 
        for attribute in data_dict:
            cur_parser=attribute_info[attribute]['parser']
            cur_vals=cur_parser([data_dict[attribute],chrom,start_pos,end_pos,attribute_info[attribute]])
            dict_to_write[attribute]=cur_vals[-1] #the last entry in the tuple is the actual numpy array of values; the first entries store start and end blocks
        payload=pickle.dumps([task_index,start_index,end_index,dict_to_write],pickle.HIGHEST_PROTOCOL)
        write_queue.put(payload)
        gc.collect()
    except: 
        kill_child_processes(os.getpid())
        raise

def write_array(args, updating, chunks_to_process):    
    try:
        #config
        tdb_Config=tiledb.Config(tdb_config_params)
        tdb_write_Context=tiledb.Ctx(config=tdb_Config)   
        
        if updating is True:
            tdb_read_Context=tiledb.Ctx(config=tdb_Config)
            cur_array_toread=tiledb.DenseArray(args.array_name,ctx=tdb_read_Context,mode='r')
        cur_array_towrite=tiledb.DenseArray(args.array_name,ctx=tdb_write_Context,mode='w')
        chunks_processed=0
        while chunks_processed < chunks_to_process:
            while write_queue.empty() is True:
                time.sleep(10)
            processed_chunk=write_queue.get()
            processed_chunk_unpickled=pickle.loads(processed_chunk)
            task_index=processed_chunk_unpickled[0]
            start_index=processed_chunk_unpickled[1]
            end_index=processed_chunk_unpickled[2]
            dict_to_write=processed_chunk_unpickled[3]
            if updating is True:
                #we are only updating some attributes in the array
                cur_vals=cur_array_toread[start_index:end_index,task_index]            
                #print("got cur vals for task "+str(task_index)+" for "+str(start_index)+":"+str(end_index))
                for key in dict_to_write:
                    cur_vals[key]=dict_to_write[key]
                dict_to_write=cur_vals
                print("updated data dict for writing:"+args.array_name) 
            else:
                #we are writing for the first time, make sure all attributes are provided, if some are not, use a nan array
                required_attrib=list(get_attribute_info(args.attribute_config,args.attribute_config_file).keys())
                #print(str(required_attrib))
                for attrib in required_attrib:
                    if attrib not in dict_to_write:
                        print("augmenting")
                        dict_to_write[attrib]=np.full(end_index-start_index,np.nan)
            #write in chunks
            cur_array_towrite[start_index:end_index,task_index]=dict_to_write
            print('Gigs:', round(psutil.virtual_memory().used / (10**9), 2))
            gc.collect()
            chunks_processed+=1
            print("wrote to disk "+str(task_index)+" for "+str(start_index)+":"+str(end_index)+";"+str(chunks_processed)+"/"+str(chunks_to_process))
        assert chunks_processed >=chunks_to_process
        print("closing arrays")
        if updating is True:
            cur_array_toread.close()
        cur_array_towrite.close()
        return 

    except KeyboardInterrupt:
        kill_child_processes(os.getpid())
        #try to delete all tmp files
        raise
    except Exception as e:
        print(e)
        kill_child_processes(os.getpid())
        raise Exception(e.message) 

    
def main():
    args=parse_args()
    ingest(args)
        
if __name__=="__main__":
    main() 
    
    
