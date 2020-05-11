db_ingest_single_threaded --tiledb_metadata tier1.encode.dnase.tasks.tsv \
	--tiledb_group db/dnase \
	--overwrite \
	--chrom_sizes hg38.chrom.sizes \
	--tile_size 10000 \
	--write_chunk 10000000
