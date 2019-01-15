# seqdataloader
Sequence data label generation and ingestion into deep learning models


Sample datasets are included in the folder `peak_files_from_encode_for_label_comparison`

Execute the script:

`generate_inputs_tiled_whole_genome_indexed.sh` for examples on how to generate classification labels on sample datasets.
The script generates binary classification labels (1,0,-1 for ambiguous) in bed file format:

http://mitra.stanford.edu/kundaje/seqdataloader/labels.50PercentOverlap.tsv.gz

http://mitra.stanford.edu/kundaje/seqdataloader/labels.SummitWithin200bpCenter.tsv.gz

Corresponding WashU Browser Tracks with optimal narrowPeak and associated bin labels are here:

http://epigenomegateway.wustl.edu/legacy/?genome=hg38&session=57a1zxqQW7&statusId=1100007190





