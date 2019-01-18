# seqdataloader
Sequence data label generation and ingestion into deep learning models


Sample datasets are included in the folder `peak_files_from_encode_for_label_comparison` and `bigwig_files_from_encode_for_label_comparison`

Execute the script:

`generate_inputs_tiled_whole_genome_indexed.sh` for examples on how to generate classification and regression labels on sample datasets.
The script generates binary classification labels (1,0,-1 for ambiguous) or continuous regression labels reflective of bigWig coverage in a bin  in bed file format:

http://mitra.stanford.edu/kundaje/seqdataloader/classificationlabels.50PercentOverlap.tsv.gz

http://mitra.stanford.edu/kundaje/seqdataloader/classificationlabels.SummitWithin200bpCenter.tsv.gz

http://mitra.stanford.edu/kundaje/seqdataloader/regressionlabels.50PercentOverlap.tsv.gz

http://mitra.stanford.edu/kundaje/seqdataloader/regressionlabels.SummitWithin200bpCenter.tsv.gz

Corresponding WashU Browser Tracks with optimal narrowPeak and associated bin labels are here:
http://epigenomegateway.wustl.edu/legacy/?genome=hg38&session=GDB2BTMGnB&statusId=1154897038


## Dependencies

Please make sure the following dependencies are installed on your system to use SeqDataLoader:
* pybedtools
* pyBigWig 
* pandas
* numpy
* csv
* multiprocessing
* gzip 



