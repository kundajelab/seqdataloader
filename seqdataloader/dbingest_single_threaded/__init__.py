## helper functions to ingest bigwig and narrowPeak data files into a tileDB instance.
## tileDB instances are indexed by coordinate
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import psutil
#import multiprocessing as mp
#mpx = mp.get_context('spawn')
import tiledb
import pdb
import argparse
import pandas as pd
import numpy as np
from collections import OrderedDict
from ..attrib_config import *
from ..utils import *
from ..tdb_config import * 
import gc

    
def args_object_from_args_dict(args_dict):
    #create an argparse.Namespace from the dictionary of inputs
    args_object=argparse.Namespace()
    #set the defaults
    vars(args_object)['overwrite']=False
    vars(args_object)['coord_tile_size']=10000
    vars(args_boject)['task_tile_size']=1
    vars(args_object)['attribute_config']='encode_pipeline'
    vars(args_object)['write_chunk']=None
    for key in args_dict:
        vars(args_object)[key]=args_dict[key]
    #set any defaults that are unset 
    args=args_object    
    return args 
    
def parse_args():
    parser=argparse.ArgumentParser(description="ingest data into tileDB")
    parser.add_argument("--tiledb_metadata",help="fields are: dataset, fc_bigwig, pval_bigwig, count_bigwig_plus_5p, count_bigwig_minus_5p, count_bigwig_unstranded_5p, idr_peak, overlap_peak, ambig_peak")
    parser.add_argument("--tiledb_group")
    parser.add_argument("--overwrite",default=False,action="store_true") 
    parser.add_argument("--chrom_sizes",help="2 column tsv-separated file. Column 1 = chromsome name; Column 2 = chromosome size")
    parser.add_argument("--coord_tile_size",type=int,default=10000,help="coordinate axis tile size")
    parser.add_argument("--task_tile_size",type=int,default=1,help="task axis tile size")
    parser.add_argument("--attribute_config",default='encode_pipeline',help="the following are supported: encode_pipeline, generic_bigwig")
    parser.add_argument("--write_chunk",type=int,default=None,help="number of bases to write to disk in one tileDB DenseArray write operation") 
    return parser.parse_args()

def create_new_array(tdb_Context,
                     size,
                     array_out_name,
                     coord_tile_size,
                     task_tile_size,
                     attribute_config,
                     compressor='gzip',
                     compression_level=-1,
                     var=False):
    '''
    Creates an empty tileDB array
    size= tuple(num_indices,num_tasks)
    '''
    coord_tile_size=min(size[0],coord_tile_size)
    task_tile_size=min(size[1],task_tile_size)
    tiledb_dim_coords = tiledb.Dim(
        name='genome_coordinate',
        domain=(0, size[0]),
        tile=coord_tile_size,
        dtype='uint32')
    tiledb_dim_tasks=tiledb.Dim(
        name='task',
        domain=(0,size[1]),
        tile=task_tile_size,
        dtype='uint32')
    tiledb_dom = tiledb.Domain(tiledb_dim_coords,tiledb_dim_tasks,ctx=tdb_Context)

    #generate the attribute information
    attribute_info=get_attribute_info(attribute_config)
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
            if cur_fname is not None:
                data_dict[col]=attribute_info[col]['opener'](cur_fname)
        return data_dict
    except Exception as e:
        print(repr(e))
        raise e
    
def get_subdict(full_dict,start,end):
    subdict=dict()
    for key in full_dict:
        subdict[key]=full_dict[key][start:end]
    print(subdict.keys())
    return subdict
    
def ingest_single_threaded(args):
    if type(args)==type({}):
        args=args_object_from_args_dict(args)

    #config
    tdb_Config=tiledb.Config(tdb_config_params)
    tdb_write_Context=tiledb.Ctx(config=tdb_Config)   
    tdb_read_Context=tiledb.Ctx(config=tdb_Config)
    
    overwrite=args.overwrite
    coord_tile_size=args.coord_tile_size
    task_tile_size=args.task_tile_size
    attribute_config=args.attribute_config
    updating=False

    attribute_info=get_attribute_info(args.attribute_config) 
    tiledb_metadata=pd.read_csv(args.tiledb_metadata,header=0,sep='\t')
    num_tasks=tiledb_metadata.shape[0]
    
    print("loaded tiledb metadata")
    chrom_sizes=pd.read_csv(args.chrom_sizes,header=None,sep='\t')
    print("loaded chrom sizes")
    chrom_indices,num_indices=transform_chrom_size_to_indices(chrom_sizes)
    print("num_indices:"+str(num_indices))
    array_out_name=args.tiledb_group
    if tiledb.object_type(array_out_name) == "array":
        if overwrite==False:
            raise Exception("array:"+str(array_out_name) + "already exists; use the --overwrite flag to overwrite it. Exiting")
        else:
            print("warning: the array: "+str(array_out_name)+" already exists. You provided the --overwrite flag, so it will be updated/overwritten")
            updating=True
    else:
        #create the array:
        create_new_array(tdb_Context=tdb_write_Context,
                         size=(num_indices,num_tasks),
                         attribute_config=attribute_config,
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
    if updating is True:
        cur_array_toread=tiledb.DenseArray(array_out_name,ctx=tdb_read_Context,mode='r')
    else:
        cur_array_toread=None
    cur_array_towrite=tiledb.DenseArray(array_out_name,ctx=tdb_write_Context,mode='w')
    for task_index,task_row in tiledb_metadata.iterrows():
        dataset=task_row['dataset']
        print(dataset) 
        #read in filenames for bigwigs
        data_dict=open_data_for_parsing(task_row,attribute_info)
        for start_chunk_index in range(0,num_indices,args.write_chunk):
            print(str(start_chunk_index)+'/'+str(num_indices)) 
            end_chunk_index=start_chunk_index+min([num_indices,start_chunk_index+args.write_chunk])
            print("end chunk index:"+str(end_chunk_index))
            #convert global indices to chrom+pos indices
            chunk_chrom_coords=transform_indices_to_chrom_coords(start_chunk_index,end_chunk_index,chrom_indices)
            print("processing:"+str(chunk_chrom_coords))
            for coord_set in chunk_chrom_coords:
                print("\t"+"coord_set:"+str(coord_set))
                process_chunk(task_index,data_dict,attribute_info,coord_set,updating,args,cur_array_toread,cur_array_towrite)
                print('Gigs:', round(psutil.virtual_memory().used / (10**9), 2))            
                print("wrote chrom array for task:"+str(dataset)+"for index:"+str(start_chunk_index))
    print("closing arrays")
    if cur_array_to_read is not None:
        cur_array_toread.close()
    cur_array_towrite.close()
    print('done!') 

def process_chunk(task_index, data_dict, attribute_info, coord_set, updating, args, cur_array_toread, cur_array_towrite):
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
        print("got:"+str(attribute)+" for task "+str(task_index)+" for "+str(chrom)+":"+str(start_pos)+"-"+str(end_pos))

    if updating is True:
        #we are only updating some attributes in the array
        cur_vals=cur_array_toread[start_index:end_index,task_index]            
        print("got cur vals for task "+str(task_index)+" for "+str(chrom)+":"+str(start_pos)+"-"+str(end_pos))
        for key in dict_to_write:
            cur_vals[key]=dict_to_write[key]
        dict_to_write=cur_vals
        print("updated data dict for writing:"+array_out_name) 
    else:
        #we are writing for the first time, make sure all attributes are provided, if some are not, use a nan array
        required_attrib=list(get_attribute_info(attribute_config).keys())
        for attrib in required_attrib:
            if attrib not in dict_to_write:
                dict_to_write[attrib]=np.full(end_pos-start_pos,np.nan)
            
    #write in chunks
    cur_array_towrite[start_index:end_index,task_index]=dict_to_write
    print("wrote to disk "+str(task_index)+" for "+str(chrom)+":"+str(start_pos)+"-"+str(end_pos))
    gc.collect() 
    
def main():
    args=parse_args()
    ingest_single_threaded(args)
        
if __name__=="__main__":
    main() 
    
    
