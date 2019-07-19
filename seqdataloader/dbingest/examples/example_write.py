import pyBigWig
import tiledb
import numpy as np


bw = pyBigWig.open('/oak/stanford/groups/akundaje/projects/alzheimers_parkinsons/merged_tagAligns_outputs/PD_CTRL_PTMN/cromwell-executions/atac/643566f0-0543-4664-b6e9-efa0c01de273/call-macs2_signal_track/shard-0/execution/glob-7ab0340dfeb10ca109917cbdcc568548/grouped.t.PD_CTRL_PTMN.gz.pval.signal.bigwig')
bw_base_name = 'ENCFF001CTB'
db_base_name = 'encode_dnase_hg38'
db_entry='/'.join([db_base_name,bw_base_name])
chrom_info = bw.chroms()
print("opened bigwig and got chrom info")

# Create the context 
ctx = tiledb.Ctx()

def write_array_to_tiledb(array, url, ctx, default_tile_size=9000,compressor='gzip', compression_level=-1):
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

for chrom_name, chrom_size in chrom_info.items():
    if chrom_name !="chr21":
        continue
    print("chrom name"+str(chrom_name)+", chrom size:"+str(chrom_size)) 
    signal_data = np.zeros(chrom_size, dtype=np.float32)
    print("made signal data") 
    signal_data[:] = bw.values(chrom_name, 0, chrom_size)
    print("populated chrom_size") 
    chrom_s3_address = db_entry+'.' + chrom_name
    print(chrom_s3_address) 
    write_array_to_tiledb(signal_data, chrom_s3_address, ctx)
    print("wrote array to disk") 
