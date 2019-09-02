## helper functions to ingest bigwig and narrowPeak data files into a tileDB instance.
## tileDB instances are indexed by coordinate
import tiledb
import argparse
import pandas as pd
import numpy as np
from .attrib_config import *
from .utils import *
from concurrent.futures import ProcessPoolExecutor
import pdb

def parse_args():
    parser=argparse.ArgumentParser(description="ingest data into tileDB")
    parser.add_argument("--tiledb_metadata",help="fileds are: Dataset, fc_bigwig, pval_bigwig, count_bigwig_plus_5p, count_bigwig_minus_5p, idr_peaks, overlap_peaks, ambig_peaks")
    parser.add_argument("--tiledb_group")
    parser.add_argument("--overwrite",default=False,action="store_true") 
    parser.add_argument("--chrom_sizes",help="2 column tsv-separated file. Column 1 = chromsome name; Column 2 = chromosome size")
    parser.add_argument("--chrom_threads",type=int,default=1,help="inner thread pool, launched by task_threads")
    parser.add_argument("--task_threads",type=int,default=1,help="outer thread pool,launched by main thread")
    parser.add_argument("--write_threads",type=int,default=1)
    parser.add_argument("--batch_size",type=int,default=1000000,help="num entries to write at once")
    parser.add_argument("--tile_size",type=int,default=10000,help="tile size") 
    return parser.parse_args()

def create_new_array(size,
                     array_out_name,
                     tile_size,
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
    #config
    tdb_Config=tiledb.Config({"sm.check_coord_dups":"false",
                              "sm.check_coord_oob":"false",
                              "sm.check_global_order":"false",
                              "sm.num_writer_threads":50,
                              "sm.num_reader_threads":50})
    tdb_Context=tiledb.Ctx(config=tdb_Config)
    tiledb_dom = tiledb.Domain(tiledb_dim,ctx=tdb_Context)
    
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

def write_chunk(inputs):
    array_out_name=inputs[0]
    start=inputs[1]
    end=inputs[2]
    sub_dict=inputs[3]
    batch_size=inputs[4]
    print("start:"+str(start)+", end:"+str(end))
    ctx = tiledb.Ctx()
    with tiledb.DenseArray(array_out_name,ctx=ctx,mode='w') as out_array: 
        #sub_df=df.iloc[start:end]
        #sub_dict=sub_df.to_dict(orient='list')
        out_array[start:end]=sub_dict
        print("done with chunk start:"+str(start)+", end:"+str(end)) 
    return "done"
    
def write_array_to_tiledb(size,
                          dict_to_write,
                          array_out_name,
                          batch_size=10000,
                          compressor='gzip',
                          compression_level=-1,
                          updating=False,
                          write_threads=1):
    print("starting to write output") 
    if updating is True:
        #we are only updating some attributes in the array
        with tiledb.DenseArray(array_out_name,mode='r') as cur_array:
            cur_vals=cur_array[:]
        print('got cur vals') 
        for key in dict_to_write:
            print(key) 
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
    df=pd.DataFrame.from_dict(dict_to_write)
    num_entries=df.shape[0]
    pool_inputs=[]
    for i in range(0,num_entries,batch_size):
        pool_inputs.append((array_out_name,i,i+batch_size,df.iloc[i:i+batch_size].to_dict(orient="list"),batch_size))
    pool_inputs.append((array_out_name,i+batch_size,num_entries,df.iloc[i:i+batch_size].to_dict(orient="list"),batch_size))
    
    with ProcessPoolExecutor(max_workers=write_threads) as pool:
        print("made pool")
        futures=pool.map(write_chunk,pool_inputs)
    print("done writing")

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
    args=inputs[5]
    overwrite=args.overwrite
    write_threads=args.write_threads
    batch_size=args.batch_size
    tile_size=args.tile_size 
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
                         array_out_name=array_out_name,
                         tile_size=tile_size)
        print("created new array:"+str(array_out_name))

    dict_to_write=dict()
    for attribute in data_dict:
        store_summits=False
        summit_indicator=None
        if 'store_summits' in attribute_info[attribute]:
            store_summits=attribute_info[attribute]['store_summits']
        print("store_summits:"+str(store_summits))
        if 'summit_indicator' in attribute_info[attribute]: 
            summit_indicator=attribute_info[attribute]['summit_indicator']
        print("summit_indicator:"+str(summit_indicator))
        dict_to_write[attribute]=attribute_info[attribute]['parser'](data_dict[attribute],chrom,size,store_summits,summit_indicator)
        print("got:"+str(attribute)+" for chrom:"+str(chrom))
    
    write_array_to_tiledb(size=size,
                          dict_to_write=dict_to_write,
                          array_out_name=array_out_name,
                          updating=updating,
                          batch_size=batch_size,
                          write_threads=write_threads)
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
    pool_inputs=[] 
    array_outf_prefix="/".join([args.tiledb_group,dataset])
    print("parsed pool inputs") 
    for index, row in chrom_sizes.iterrows():
        chrom=row[0]
        size=row[1]
        array_out_name='.'.join([array_outf_prefix,chrom])
        pool_inputs.append((chrom,size,array_out_name,data_dict,attribute_info,args))
    try:
        with ProcessPoolExecutor(max_workers=args.chrom_threads) as pool:
            cur_futures=pool.map(process_chrom,pool_inputs)
        for entry in pool_inputs:
            result=process_chrom(entry)
    except Exception as e:
        raise e 
    return "done"

def args_object_from_args_dict(args_dict):
    #create an argparse.Namespace from the dictionary of inputs
    args_object=argparse.Namespace()
    #set the defaults
    vars(args_object)['chrom_threahds']=1
    vars(args_object)['task_threads']=1
    vars(args_object)['overwrite']=False
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
    try:
        with ProcessPoolExecutor(max_workers=args.task_threads) as pool:
            results=pool.map(create_tiledb_array,pool_inputs)
        #the lines below are kept in case we need to get rid of the pooled multiprocessing
        # at any point 
        #for entry in pool_inputs:
        #    result=create_tiledb_array(entry)
    except Exception as e:
        raise e
    return "done"

def main():
    args=parse_args()
    ingest(args)
    
    
if __name__=="__main__":
    main() 
    
    
