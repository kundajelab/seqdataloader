#Classification Approach 1: Summit Must Lie Within 200 BP Bin 
#python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
#       --out_bed classificationlabels.SummitWithin200bpCenter.tsv.gz \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --bin_size 200 \
#       --threads 1 \
#       --allow_ambiguous \
#       --labeling_approach peak_summit_in_bin_classification
#
##Classification Approach 2: 50% Overlap Between Peak and 200 BP Bin (50% of the Smaller of the Two) 
#python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
#       --out_bed classificationlabels.50PercentOverlap.tsv.gz \
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
#       --out_bed regressionlabels.SummitWithin200bpCenter.tsv.gz \
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
#       --out_bed regressionlabels.50PercentOverlap.tsv.gz \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 4 \
#       --allow_ambiguous \
#       --overlap_thresh 0.5 \
#       --labeling_approach peak_percent_overlap_with_bin_regression

##Regression Approach 3: Provide bedtools coverage in the bigWig for every bin in the genome 
#python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
#       --out_bed regressionlabels.allbins.hg38.tsv.gz \
#       --chrom_sizes hg38.chrom.sizes \
#       --bin_stride 50 \
#       --left_flank 400 \
#       --right_flank 400 \
#       --threads 4 \
#       --subthreads 10 \
#       --allow_ambiguous \
#       --overlap_thresh 0.5 \
#       --labeling_approach all_genome_bins_regression
#
#
