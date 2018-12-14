#Approach 1: Summit Must Lie Within 200 BP Bin 
python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
       --out_bed labels.SummitWithin200bpCenter.tsv.gz \
       --chrom_sizes hg38.chrom.sizes \
       --bin_stride 50 \
       --left_flank 400 \
       --right_flank 400 \
       --bin_size 200 \
       --threads 1 \
       --allow_ambiguous \
       --labeling_approach peak_summit_in_bin

#Approach 2: 50% Overlap Between Peak and 200 BP Bin (50% of the Smaller of the Two) 
python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
       --out_bed labels.50PercentOverlap.tsv.gz \
       --chrom_sizes hg38.chrom.sizes \
       --bin_stride 50 \
       --left_flank 400 \
       --right_flank 400 \
       --threads 1 \
       --allow_ambiguous \
       --overlap_thresh 0.5 \
       --labeling_approach peak_percent_overlap_with_bin

