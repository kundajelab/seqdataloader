# seqdataloader
Sequence data label generation and ingestion into deep learning models


Sample datasets are included in the folder `peak_files_from_encode_for_label_comparison`

Execute the script:
`generate_inputs_tiled_whole_genome_indexed.sh` for examples on hot to generate classification labels on sample datasets.
The script generates binary classification labels (1,0,-1 for ambiguous) in bed file format:labels.SummitWithin200bpCenter.tsv.gz and labels.50PercentOverlap.tsv.gz

