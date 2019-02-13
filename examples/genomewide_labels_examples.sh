#Classification Approach 1: Summit Must Lie Within 200 BP Bin
## Timing 
## writing to gzip 
##real 11m46.403s
##user 18m40.788s
##sys 6m18.136s

## Writing to bz2:
## real 14m44.037s
## user 21m49.384s
## sys 6m21.000s

#genomewide_labels --task_list tasks.tsv \
#		  --outf classificationlabels.SummitWithin200bpCenter.tsv.gz \
#		  --output_type gzip \
#		  --chrom_sizes hg38.chrom.sizes \
#		  --bin_stride 50 \
#		  --left_flank 400 \
#		  --right_flank 400 \
#		  --bin_size 200 \
#		  --threads 10 \
#		  --subthreads 4 \
#		  --allow_ambiguous \
#		  --labeling_approach peak_summit_in_bin_classification

#Example of restricting analysis to a single chromosome with --chroms_to_keep flag 
#genomewide_labels --task_list tasks.tsv \
#		  --outf classificationlabels.SummitWithin200bpCenter.tsv.gz \
#		  --output_type gzip \
#		  --chrom_sizes hg38.chrom.sizes \
#		  --chroms_to_keep chr21 \
#		  --bin_stride 50 \
#		  --left_flank 400 \
#		  --right_flank 400 \
#		  --bin_size 200 \
#		  --threads 10 \
#		  --subthreads 4 \
#		  --allow_ambiguous \
#		  --labeling_approach peak_summit_in_bin_classification


#Example with only positives stored
genomewide_labels --task_list tasks.tsv \
		  --outf classificationlabels.SummitWithin200bpCenter.tsv.gz \
		  --output_type gzip \
		  --chrom_sizes hg38.chrom.sizes \
		  --chroms_to_keep chr21 \
		  --bin_stride 50 \
		  --left_flank 400 \
		  --right_flank 400 \
		  --bin_size 200 \
		  --threads 10 \
		  --subthreads 4 \
		  --allow_ambiguous \
		  --store_positives_only \
		  --labeling_approach peak_summit_in_bin_classification



##Classification Approach 2: 50% Overlap Between Peak and 200 BP Bin (50% of the Smaller of the Two)
## Timing 
## real18m56.337s
## user25m23.004s
## sys7m58.104s

#genomewide_labels --task_list tasks.tsv \
#       --outf classificationlabels.50PercentOverlap.tsv.gz \
#       --output_type gzip \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 10 \
#       --subthreads 4 \
#       --allow_ambiguous \
#       --overlap_thresh 0.5 \
#       --labeling_approach peak_percent_overlap_with_bin_classification

##Regression Approach 1:Summit Must Lie Within 200 BP Bin
## Timing:
## real 18m15.728s
## user 24m25.028s
## sys 7m58.244s

#genomewide_labels --task_list tasks.tsv \
#       --outf regressionlabels.SummitWithin200bpCenter.tsv.gz \
#       --output_type gzip \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --bin_size 200 \
#       --threads 10 \
#       --subthreads 4 \
#       --allow_ambiguous \
#       --labeling_approach peak_summit_in_bin_regression
#

#Regression Approach 2: 50% Overlap Between Peak and 200 BP Bin (50% of the Smaller of the Two)
## real18m56.337s
## user25m23.004s
## sys7m58.104s

#genomewide_labels --task_list tasks.tsv \
#       --outf regressionlabels.50PercentOverlap.tsv.gz \
#       --output_type gzip \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 10 \
#       --subthreads 4 \
#       --allow_ambiguous \
#       --overlap_thresh 0.5 \
#       --labeling_approach peak_percent_overlap_with_bin_regression


##Regression Approach 3: Provide bedtools coverage in the bigWig for every bin in the genome

##Timing for hdf5 save 
## real 8m51.275s
## user 17m38.576s
## sys 6m14.768s
#genomewide_labels --task_list tasks.tsv \
#       --outf regressionlabels.allbins.hg38.hdf5 \
#       --output_type hdf5 \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 24 \
#       --subthreads 2 \
#       --labeling_approach all_genome_bins_regression

## Timing (pkl) 
## real    23m10.448s
## user    31m55.056s
## sys     5m39.880s
#genomewide_labels --task_list tasks.tsv \
#       --outf regressionlabels.allbins.hg38.pkl \
#       --output_type pkl \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 24 \
#       --subthreads 2 \
#       --labeling_approach all_genome_bins_regression


## Timing for full data frame (gzip) 
##  real29m50.597s
##  user38m2.020s
##  sys6m34.064s

## Timing for chromosome-specific dataframes (gzip) 
## real 21m35.525s
## user 51m55.496s
## sys 7m49.140s
#genomewide_labels --task_list tasks.tsv \
#       --outf regressionlabels.allbins.hg38.tsv.gz \
#       --output_type gzip \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 24 \
#       --subthreads 2 \
#       --split_output_by_chrom \
#       --labeling_approach all_genome_bins_regression
