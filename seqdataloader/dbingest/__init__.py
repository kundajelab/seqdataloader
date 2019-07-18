## helper functions to ingest bigwig and narrowPeak data files into a tileDB instance.
## tileDB instances are indexed by coordinate, and a task(i.e. a unique identifier) 
## attributes include cell type, assay, project, filetype (i.e. fc bigwig, pval bigwig, 
import tiledb
import argparse
import pandas as pdb

def parse_args():
    parser=argparse.ArgumentParser(description="ingest data into tileDB")
    parser.add_argument("--input_files",help="2 columns file, with unique identifiers in column 1 and file names in column 2")
    parser.add_argument("--tile_db_name")
    parser.add_argument("--append",default=False,action="store_true") 
    parser.add_argument("--chrom_sizes")
    parser.add_argument("--dense",default=False,action="store_true")
    return parser.parse_args()

def create_tiledb_array(chrom_sizes,tasks) 
    '''
    create new tileDB array 
    '''
    chroms=pd.read_tsv(chrom_sizes,header=None,sep='\t')
    tasks=pd.read_tsv(tasks,header=0,sep='\t')
    task_names=tasks[0]
    ctx=tiledb.Ctx()
    doms=[] 
    for index,row in chroms.iterrows():
        chrom_name=row[0]
        chrom_size=row[1]
        new_dom=tiledb.Domain(ctx,
                              tiledb.Dim(ctx,name="coord",domain=(1,chrom_size),tile=chrom_size,dtype=int.32))
        #schema will have an attribute for each task
        for task in task_names:
            schema=tiledb
                            

    
    pass

def ingest_bigwig():
    pass

def ingest_narrowPeak():
    pass 

def create_task():
    pass

def db_ingest(args): 
    pass


def main():
    args=parse_args()
    db_ingest(args)
    
if __name__=="__main__":
    main() 
    
    
