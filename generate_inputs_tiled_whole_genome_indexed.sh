#Classification Approach 1: Summit Must Lie Within 200 BP Bin 
python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
       --outf classificationlabels.SummitWithin200bpCenter.tsv.gz \
       --output_type bed.gz \
       --chrom_sizes hg38.chrom.sizes \
       --bin_stride 50 \
       --left_flank 400 \
       --right_flank 400 \
       --bin_size 200 \
       --threads 10 \
       --subthreads 4 \
       --allow_ambiguous \
       --labeling_approach peak_summit_in_bin_classification

##Classification Approach 2: 50% Overlap Between Peak and 200 BP Bin (50% of the Smaller of the Two) 
#python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
#       --outf classificationlabels.50PercentOverlap.tsv.gz \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 1 \
#       --allow_ambiguous \
#       --overlap_thresh 0.5 \
#       --labeling_approach peak_percent_overlap_with_bin_classification
#
##Regression Approach 1:Summit Must Lie Within 200 BP Bin
#python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
#       --outf regressionlabels.SummitWithin200bpCenter.tsv.gz \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --bin_size 200 \
#       --threads 4 \
#       --allow_ambiguous \
#       --labeling_approach peak_summit_in_bin_regression
#

##Regression Approach 2: 50% Overlap Between Peak and 200 BP Bin (50% of the Smaller of the Two) 
#python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
#       --outf regressionlabels.50PercentOverlap.tsv.gz \
#       --output_type bed.gz \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 4 \
#       --subthreads 2 \
#       --allow_ambiguous \
#       --overlap_thresh 0.5 \
#       --labeling_approach peak_percent_overlap_with_bin_regression
#

##Regression Approach 3: Provide bedtools coverage in the bigWig for every bin in the genome
### Timing:
## real    23m10.448s
## user    31m55.056s
## sys     5m39.880s
#python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
#       --outf regressionlabels.allbins.hg38.pkl \
#       --output_type pkl \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 24 \
#       --subthreads 2 \
#       --labeling_approach all_genome_bins_regression
#


## Timing 
##  real29m50.597s
##  user38m2.020s
##  sys6m34.064s
## Chromosome-specific io:
## real21m35.525s
## user51m55.496s
## sys7m49.140s
#python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
#       --outf regressionlabels.allbins.hg38.tsv.gz \
#       --output_type bed.gz \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 24 \
#       --subthreads 2 \
#       --split_output_by_chrom \
#       --labeling_approach all_genome_bins_regression

##TIMING:
## real 8m51.275s
## user 17m38.576s
## sys 6m14.768s
#python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
#       --outf regressionlabels.allbins.hg38.hdf5 \
#       --output_type hdf5 \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 24 \
#       --subthreads 2 \
#       --labeling_approach all_genome_bins_regression
#

