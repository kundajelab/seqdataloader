base=/users/amtseng/test/tiledb_bench/SPI1_test

#db_ingest_single_threaded --tiledb_metadata $base/tasklist.tsv \
#	--tiledb_group amtseng_database \
#	--overwrite \
#	--chrom_sizes $base/hg38.chr10.size \
#	--tile_size 10000 \
#	--write_chunk 1000000

db_ingest_single_threaded --tiledb_metadata $base/tasklist.tsv \
	--tiledb_group amtseng_database_10million_write_chunk \
	--overwrite \
	--chrom_sizes $base/hg38.chr10.size \
	--tile_size 10000 \
	--write_chunk 10000000
