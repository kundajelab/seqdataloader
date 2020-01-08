db_ingest --tiledb_metadata tasklist.tsv \
	--tiledb_group highmem_50million_write_chunk \
	--overwrite \
	--chrom_sizes hg38.chr10.size \
	--task_threads 2 \
	--chrom_threads 25 \
	--write_threads 1 \
	--attribute_config encode_pipeline
