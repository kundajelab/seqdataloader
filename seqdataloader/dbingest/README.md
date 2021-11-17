
## Example command: 

To ingest a microglia dataset with an unstranded bigwig, idr peaks, optimal peaks, 
a blacklist of genomic regions to avoid, and negative (Non-peak) regions that are gc-matched to the idr peak set 
we would run dbingest as follows: 

contents of metadata.tsv: 

```
dataset     idr_peak                            overlap_peak                            ambig_peak              negatives_peak  count_bigwig_unstranded_5p
microglia   microglia.idr.optimal.narrowPeak    microglia.overlap.optimal.narrowPeak     blacklist/GRch38/GRch38_unified_blacklist.bed     microglia.gc.matched.negatives.bed ../data/microglia.unstranded.bw
```

contents of attributes.txt: 

```
idr_peak        bed_summit_from_last_col
overlap_peak    bed_summit_from_last_col
ambig_peak      bed_no_summit
negatives_peak  bed_no_summit
count_bigwig_unstranded_5p      bigwig
```

The attributes file indicates how each column in the metadata.tsv file should be parsed. Supported values are  

* bed_summit_from_last_col -- this assumes the input file is in narrowPeak (or similar) format, where the summit offset from the start coordinate is in the last column of the file.  File is stored as an array of 0 (no peak) 1 (peak) and 2 (summit). 

* bed_no_summit -- this assumes that the input file is a bed file without summit information -- peak intervals are centered on (start+end)/2   

* bigwig -- treat the input file as a bigwig   

* bed_no_summit -- do not calculate summits for the provided bed file (i.e. store the bed file as an array of 0 (no peak) and 1 (peak) but don't store a value of 2 to indicate summit). 

The command to run to ingest the metadata.csv file to tiledb is: 

```

db_ingest --tiledb_metadata metadata.tsv \
          --array_name microglia_db \
          --overwrite \
          --chrom_sizes hg38.chrom.sizes \
          --attribute_config_file attribs.txt \
          --coord_tile_size 10000 \
          --task_tile_size 1 \
          --write_chunk 30000000 \
          --threads 40 \
          --max_queue_size 50 \
          --max_mem_g 200
```


