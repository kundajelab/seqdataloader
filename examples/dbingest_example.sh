#!/bin/bash
db_ingest --tiledb_metadata tasks.dbingest.tsv \
    --tiledb_group hepg2_dnase_encode \
    --overwrite \
    --chrom_sizes hg38.chrom21.sizes \
    --chrom_threads 25 \
    --task_threads 1 \
    --store_summits \
    --write_threads 1 \
    --summit_indicator 2

