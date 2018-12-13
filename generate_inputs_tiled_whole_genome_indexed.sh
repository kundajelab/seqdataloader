#Approach 1: Summit Must Lie Within 200 BP Around Bin Center 
python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
       --out_bed labels.SummitWithin200bpCenter.tsv.gz \
       --chrom_sizes hg38.chrom.sizes \
       --stride 50 \
       --bin_size 1000 \
       --bin_center_size 200 \
       --threads 1 \
       --allow_ambiguous \
       --labeling_approach peak_summit_near_bin_center

#Approach 2: 50% Overlap Between Peak and Bin (50% of the Smaller of the Two) 
python generate_inputs_tiled_whole_genome_indexed.py --task_list tasks.tsv \
       --out_bed labels.50PercentOverlap.tsv.gz \
       --chrom_sizes hg38.chrom.sizes \
       --stride 50 \
       --bin_size 1000 \
       --threads 1 \
       --allow_ambiguous \
       --overlap_thresh 0.5 \
       --labeling_approach peak_percent_overlap_with_bin

