#!/bin/bash
db_ingest --tiledb_metadata tasks.dbingest.local.tsv \
    --tiledb_group hepg2_dnase_encode \
    --overwrite \
    --chrom_sizes hg38.chrom.sizes \
    --chrom_threads 25 \
    --attribute_config encode_pipeline \
    --tile_size 9000 \
    --batch_size 1000000

