## helper functions to ingest bigwig and narrowPeak data files into a tileDB instance.
## tileDB instances are indexed by coordinate, and a task(i.e. a unique identifier) 
## attributes include cell type, assay, project, filetype (i.e. fc bigwig, pval bigwig, 
import tiledb
import argparse
import pandas as pdb
import pyBigWig

def parse_args():
    parser=argparse.ArgumentParser(description="ingest data into tileDB")
    parser.add_argument("--tiledb_metadata",help="fileds are: Dataset, FC_bigwig, PVAL_bigwig, COUNT_bigwig, IDR_peaks, OVERLAP_peaks, AMBIG_peaks")
    parser.add_argument("--tiledb_group")
    parser.add_argument("--overwrite",default=False,action="store_true") 
    parser.add_argument("--chrom_sizes",help="2 column tsv-separated file. Column 1 = chromsome name; Column 2 = chromosome size")
    return parser.parse_args()


def write_array_to_tiledb(fc_bigwig_chrom_array=fc_bigwig_chrom_array,
                              pval_bigwig_chrom_array=pval_bigwig_chrom_array,
                              count_bigwig_chrom_array=count_bigwig_chrom_array,
                              idr_peak_chrom_array=idr_peak_chrom_array,
                              overlap_peak_chrom_array=overlap_peak_chrom_array,
                              ambig_peak_chrom_array=ambig_peak_chrom_array,
                              array_out_name)

def write_array_to_tiledb(array, url, default_tile_size=9000,compressor='gzip', compression_level=-1):
    size = array.shape[0]
    tile_size = min(size, default_tile_size)
    print("tile_size:"+str(tile_size))
    tiledb_dim = tiledb.Dim(
        ctx=ctx,
        name='genome_coordinate',
        domain=(0, size - 1),
        tile=tile_size,
        dtype='uint32')
    print("got dim info") 
    tiledb_dom = tiledb.Domain(tiledb_dim,ctx=ctx)
    print("made dom")
    tiledb_attr = tiledb.Attr(
        ctx=ctx,
        name='signal_value',
        filters=tiledb.FilterList([tiledb.GzipFilter()]),
        dtype='float32')
    print("made attr") 
    tiledb_schema = tiledb.ArraySchema(
        ctx=ctx,
        domain=tiledb_dom,
        attrs=(tiledb_attr,),
        cell_order='row-major',
        tile_order='row-major')
    print("made schema") 
    tiledb.DenseArray.create(url, tiledb_schema)
    print("created empty array on disk") 
    with tiledb.DenseArray(url, ctx=ctx, mode='w') as final_array:
        final_array[:] = array.astype(np.float32)

def extract_metadata_field(dataset,row,field):        
    try:
        return row[field]
    except:
        print("tiledb_metadata has no column "+field+" for dataset:"+str(dataset))
        return None

def create_tiledb_array(row,args,chrom_sizes) 
    '''
    create new tileDB array for a single dataset, overwrite if array exists and user sets --overwrite flag
    '''
    dataset=row['dataset']
    fc_bigwig=extract_metadata_field(dataset,row,'fc_bigwig')
    pval_bigwig=extract_metadata_field('pval_bigwig')
    count_bigwig_=extract_metadata_field('count_bigwig')
    #open the bigwigs if they are not none
    if fc_bigwig_ is not None:
        fc_bigwig=pyBigWig.open(fc_bigwig)
    if pval_bigwig is not None:
        pval_bigwig=pyBigWig.open(pval_bigwig)
    if count_bigwig is not None:
        count_bigwig=pyBigWig.open(count_bigwig)
        
    #load the bed files (0-indexed by default) if they are not None 
    idr_peaks=extract_metadata_field('idr_peaks')
    if idr_peaks is not None:
        idr_peaks=pd.read_csv(idr_peaks,header=None,sep='\t') 
    overlap_peaks=extract_metadata_field('overlap_peaks')
    if overlap_peaks is not None:
        overlap_peaks=pd.read_csv(overlap_peaks,header=None,sep='\t') 
    ambig_peaks=extract_metadata_field('ambig_peaks')
    if ambig_peaks is not None:
        ambig_peaks=pd.read_csv(ambig_peaks,header=None,sep='\t')
        

    array_outf_prefix="/".join([args.tiledb_group,dataset])
    for index, row in chrom_sizes.iterrows():
        chrom=row[0]
        size=row[1]
        array_out_name='.'.join([array_outf_prefix,chrom])
        if tiledb.object_type(array_out_name) == "array":
            if args.overwrite==False:
                raise Exception("array:"+str(array_out_name) + "already exists; use the --overwrite flag to overwrite it. Exiting")
            else:
                print("warning: the array: "+str(array_out_name)+" already exists. You provided the --overwrite flag, so it will be overwritten")
        fc_bigwig_chrom_array=parse_bigwig_chrom_vals(fc_bigwig,chrom,size)
        print("parsed fc bigwig for chrom:"+str(chrom))
        pval_bigwig_chrom_array=parse_bigwig_chrom_vals(pval_bigwig,chrom,size)
        print("parsed pval bigwig for chrom:"+str(chrom))
        count_bigwig_chrom_array=parse_bigwig_chrom_vals(count_bigwig,chrom,size)
        print("parsed count bigwig for chrom:"+str(chrom))
        idr_peak_chrom_array=parse_narrowPeak_chrom_vals(idr_peaks,chrom,size)
        print("parsed idr peaks for chrom:"+str(chrom))
        overlap_peak_chrom_array=parse_narrowPeak_chrom_vals(overlap_peaks,chrom,size)
        print("parsed overlap peaks for chrom:"+str(chrom))
        ambig_peak_chrom_array=parse_narrowPeak_chrom_vals(ambig_peaks,chrom,size)
        print("parsed ambig peaks for chrom:"+str(chrom))

        write_array_to_tiledb(fc_bigwig_chrom_array=fc_bigwig_chrom_array,
                              pval_bigwig_chrom_array=pval_bigwig_chrom_array,
                              count_bigwig_chrom_array=count_bigwig_chrom_array,
                              idr_peak_chrom_array=idr_peak_chrom_array,
                              overlap_peak_chrom_array=overlap_peak_chrom_array,
                              ambig_peak_chrom_array=ambig_peak_chrom_array,
                              array_out_name)
        print("wrote array to disk for dataset:"+str(dataset))         


def parse_bigwig_chrom_vals(bigwig_object,chrom,size):
    if bigwig_object is None:
        return bigwig_object
    signal_data = np.zeros(size, dtype=np.float32)
    signal_data[:]=biwgig_object.values(chrom,0,size) 

def parse_narrowPeak_chrom_vals(narrowPeak_df,chrom,size):
    if narrowPeak_df is None:
        return narrowPeak_df
    signal_data = np.zeros(size, dtype=np.float32)
    chrom_subset_df=narrowPeak_df[narrowPeak_df[0]==chrom]
    for index,row in chrom_subset_df.iterrows():
        signal_data[row[1]:row[2]]=1
    return signal_data 

def main():
    args=parse_args()
    
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
        create_tiledb_array(row,args,chrom_sizes)
    
if __name__=="__main__":
    main() 
    
    
