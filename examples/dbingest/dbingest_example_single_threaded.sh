#!/bin/bash
db_ingest_single_threaded --tiledb_metadata tasks.dbingest.local.tsv \
    --tiledb_group hepg2_dnase_encode_single_threaded \
    --overwrite \
    --chrom_sizes hg38.chrom.sizes \
    --attribute_config encode_pipeline \
    --tile_size 9000

