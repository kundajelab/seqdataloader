#save output as tsv.gz
genomewide_labels --task_list tasks.labelgen.tsv \
		  --outf classificationlabels.SummitWithin200bpCenter.tsv.gz \
		  --output_type gzip \
		  --chrom_sizes hg38.chrom21.sizes \
		  --bin_stride 50 \
		  --left_flank 400 \
		  --right_flank 400 \
		  --bin_size 200 \
		  --chrom_threads 10 \
		  --task_threads 4 \
		  --allow_ambiguous \
		  --labeling_approach peak_summit_in_bin_classification \
		  --save_label_source
#save output as hdf5
genomewide_labels --task_list tasks.labelgen.tsv \
		  --outf classificationlabels.SummitWithin200bpCenter.hdf5 \
		  --output_type hdf5 \
		  --chrom_sizes hg38.chrom21.sizes \
		  --bin_stride 50 \
		  --left_flank 400 \
		  --right_flank 400 \
		  --bin_size 200 \
		  --chrom_threads 10 \
		  --task_threads 4 \
		  --allow_ambiguous \
		  --labeling_approach peak_summit_in_bin_classification \
		  --save_label_source

